rule pairs_to_cool:
    input:
        pairs=f"{OUTDIR}/pairs/filtered/{{sample}}.filtered.pairs.gz"
    output:
        cool=f"{OUTDIR}/matrices/{{sample}}.cool"
    params:
        binsize=lambda wc: sample_base_resolution(wc.sample),
        chromsizes=lambda wc: REF["chrom_sizes"],
        assembly=ASSEMBLY
    threads: int(THREADS.get("cooler", 8))
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/cooler/load/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.cool}) $(dirname {log})
        assembly_arg=""
        if [ -n "{params.assembly}" ]; then
          assembly_arg="--assembly {params.assembly}"
        fi
        cooler cload pairs $assembly_arg -c1 2 -p1 3 -c2 4 -p2 5 {params.chromsizes}:{params.binsize} {input.pairs} {output.cool} \
          > {log} 2>&1
        """


rule zoomify_mcool:
    input:
        cool=f"{OUTDIR}/matrices/{{sample}}.cool"
    output:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    params:
        resolutions=lambda wc: sample_resolutions_csv(wc.sample),
        threads=ZOOMIFY_THREADS
    threads: ZOOMIFY_THREADS
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/cooler/zoomify/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.mcool}) $(dirname {log})
        # NOTE: no --balance here. ICE weights are computed by the single
        # balance_mcool rule below (with explicit convergence checks), so that
        # no two rules ever write into the same mcool concurrently.
        cooler zoomify --nproc {params.threads} --resolutions {params.resolutions} {input.cool} -o {output.mcool} \
          > {log} 2>&1

        # cooler zoomify does NOT propagate the genome-assembly attribute from
        # the .cool root to the .mcool root or its resolution groups; stamp it
        # explicitly at both levels (P0-3 / P1-7).
        python - <<'PY' >> {log} 2>&1
import h5py
assembly = {ASSEMBLY!r}
if assembly:
    with h5py.File("{output.mcool}", "a") as f:
        f.attrs["genome-assembly"] = assembly
        for res in list(f.get("resolutions", {{}})):
            g = f["resolutions/" + res]
            g.attrs["genome-assembly"] = assembly
    print("stamped genome-assembly:", assembly)
PY
        """


rule balance_mcool:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    output:
        # a separate, balanced copy of the mcool (the input .mcool is never
        # modified: no hidden side effects, clean provenance (P1-6))
        balanced=f"{OUTDIR}/matrices/{{sample}}.balanced.mcool",
        marker=f"{OUTDIR}/matrices/{{sample}}.mcool.balanced.done",
        summary=f"{OUTDIR}/report/balance/{{sample}}.balance.tsv"
    params:
        base=lambda wc: sample_base_resolution(wc.sample),
        resolutions=lambda wc: sorted({sample_base_resolution(wc.sample), FEATURE_RESOLUTION} | set(BALANCE_RESOLUTIONS) | set(CONCORDANCE_RESOLUTIONS)),
        max_iters=BALANCE_MAX_ITERS,
        cis_only_fallback=BALANCE_CIS_ONLY_FALLBACK,
        assembly=ASSEMBLY
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/cooler/balance/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.balanced}) $(dirname {output.marker}) $(dirname {output.summary}) $(dirname {log})

        # Single serial point where ICE weights are computed. The input mcool
        # is copied first and only the copy is modified, so no two rules ever
        # write the same HDF5 file and provenance stays clean.
        # Convergence is enforced: --convergence-policy discard drops weights
        # if not converged, then an explicit HDF5 attribute check fails the
        # rule instead of letting downstream silently use non-converged
        # weights. Optionally falls back to cis-only balancing; the mode,
        # variance, masked-bin count and weight checksum are recorded per
        # resolution (P1-2).
        python - <<'PY' > {log} 2>&1
import subprocess, sys, shutil, hashlib
import numpy as np
import cooler, h5py

mcool = "{input.mcool}"
balanced = "{output.balanced}"
resolutions = {params.resolutions!r}
max_iters = {params.max_iters}
cis_only_fallback = {params.cis_only_fallback!r}
assembly = {params.assembly!r}

shutil.copyfile(mcool, balanced)

def convergence_status(path, res):
    # read the bins/weight HDF5 dataset attrs written by cooler balance
    # (converged/var may be scalars or per-bin arrays for cis-only balancing)
    with h5py.File(path, "r") as f:
        key = f"resolutions/{{res}}/bins/weight"
        if key not in f:
            return False, False, None
        ds = f[key]
        cv = np.asarray(ds.attrs.get("converged", False))
        converged = bool(cv.all()) if cv.size > 0 else False
        var = np.asarray(ds.attrs.get("var", -1))
        if var.size == 1:
            variance = float(var.ravel()[0])
        else:
            finite = var[np.isfinite(var)]
            variance = float(np.nanmean(finite)) if finite.size else float("nan")
    return True, converged, variance

def masked_bins_and_checksum(path, res):
    # number of non-finite weights + md5 of the weight vector
    with h5py.File(path, "r") as f:
        key = f"resolutions/{{res}}/bins/weight"
        if key not in f:
            return None, None
        w = f[key][...]
    w = np.asarray(w, dtype=float)
    n_masked = int(np.isnan(w).sum())
    checksum = hashlib.md5(np.ascontiguousarray(w)).hexdigest()
    return n_masked, checksum

def drop_weight(path, res):
    with h5py.File(path, "a") as f:
        key = f"resolutions/{{res}}/bins/weight"
        if key in f:
            del f[key]

def stamp_assembly(path, assembly):
    if not assembly:
        return
    with h5py.File(path, "a") as f:
        f.attrs["genome-assembly"] = assembly
        for res in list(f.get("resolutions", {{}})):
            g = f["resolutions/" + res]
            g.attrs["genome-assembly"] = assembly

# remove pre-existing weights from the copy at every resolution (the copy
# starts unbalanced; nothing may silently reuse stale weights)
all_res = [int(r.rsplit("/", 1)[-1]) for r in cooler.fileops.list_coolers(balanced)]
target = set(resolutions)
for res in all_res:
    drop_weight(balanced, res)

ok = True
rows = []
for res in sorted(resolutions):
    uri = f"{{balanced}}::/resolutions/{{res}}"
    drop_weight(balanced, res)  # paranoia: cooler balance errors if weight exists
    mode = "genome-wide"
    cmd = ["cooler", "balance", "--max-iters", str(max_iters),
           "--convergence-policy", "discard", uri]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"[balance {{res}}bp] command failed: {{r.stderr}}", file=sys.stderr)
        ok = False
        break
    has_weight, converged, variance = convergence_status(balanced, res)
    if not (has_weight and converged):
        print(f"[balance {{res}}bp] genome-wide not converged (weight={{has_weight}}, converged={{converged}})")
        if cis_only_fallback:
            print(f"[balance {{res}}bp] retrying with --cis-only")
            drop_weight(balanced, res)
            cmd2 = ["cooler", "balance", "--cis-only", "--max-iters", str(max_iters),
                    "--convergence-policy", "discard", uri]
            r2 = subprocess.run(cmd2, capture_output=True, text=True)
            has_weight, converged, variance = convergence_status(balanced, res)
            if r2.returncode != 0 or not (has_weight and converged):
                print(f"[balance {{res}}bp] cis-only also failed (weight={{has_weight}}, converged={{converged}})", file=sys.stderr)
                ok = False
                break
            mode = "cis-only"
            print(f"[balance {{res}}bp] cis-only converged (var={{variance:.3g}})")
        else:
            print(f"[balance {{res}}bp] not converged and cis-only fallback disabled", file=sys.stderr)
            ok = False
            break
    n_masked, checksum = masked_bins_and_checksum(balanced, res)
    rows.append((res, mode, variance, n_masked, checksum))
    print(f"[balance {{res}}bp] {{mode}} converged (var={{variance:.3g}}, masked={{n_masked}})")

if not ok:
    raise SystemExit(1)

stamp_assembly(balanced, assembly)
stamp_assembly(mcool, assembly)  # also stamp the (unbalanced) input for browsing

with open("{output.summary}", "w") as h:
    h.write("resolution\tmode\tconverged\tvariance\tn_masked_bins\tweight_md5\n")
    for res, mode, variance, n_masked, checksum in rows:
        h.write("\t".join(map(str, [res, mode, True, variance, n_masked, checksum])) + "\n")

with open("{output.marker}", "w") as h:
    h.write("sample\t{wildcards.sample}\n")
    h.write("balanced_mcool\t{output.balanced}\n")
    for res, mode, variance, n_masked, checksum in rows:
        h.write(f"{{res}}bp\t{{mode}}\tvar={{variance:.3g}}\tmasked={{n_masked}}\n")
print("balance_mcool done for {wildcards.sample}")
PY
        """


rule optional_hic_export:
    input:
        # export from the balanced mcool; depends on the balance marker so it
        # can never race with the balance rule (P1-6)
        mcool=f"{OUTDIR}/matrices/{{sample}}.balanced.mcool",
        balanced=f"{OUTDIR}/matrices/{{sample}}.mcool.balanced.done"
    output:
        hic=f"{OUTDIR}/matrices/{{sample}}.hic"
    params:
        enabled=EXPORT_HIC,
        threads=HIC_EXPORT_THREADS,
        tmpdir=HIC_EXPORT_TMPDIR
    threads: HIC_EXPORT_THREADS
    conda:
        "envs/hic_tools.yaml"
    log:
        f"logs/hic_export/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.hic}) $(dirname {log})
        if [ "{params.enabled}" = "True" ]; then
          tmp_opt=""
          if [ -n "{params.tmpdir}" ]; then
            mkdir -p "{params.tmpdir}"
            tmp_opt="--tmpdir {params.tmpdir}"
          fi
          hictk convert {input.mcool} {output.hic} --threads {params.threads} $tmp_opt > {log} 2>&1
        else
          echo "hic export disabled" > {log}
          touch {output.hic}
        fi
        """


rule multiqc:
    input:
        expand(f"{OUTDIR}/qc/fastp/{{sample}}/fastp.html", sample=SAMPLES),
        expand(f"{OUTDIR}/stats/pairtools/{{sample}}.filtered.stats.txt", sample=SAMPLES),
        expand(f"{OUTDIR}/stats/pairtools/{{sample}}.dedup.stats.txt", sample=SAMPLES),
        f"{OUTDIR}/qc/multiqc_custom/microc_library_qc_mqc.yaml",
        f"{OUTDIR}/qc/multiqc_custom/microc_matrix_concordance_mqc.yaml",
        f"{OUTDIR}/qc/multiqc_custom/microc_boundary_overlap_mqc.yaml",
        f"{OUTDIR}/qc/multiqc_custom/microc_compartment_correlation_mqc.yaml"
    output:
        html=f"{OUTDIR}/qc/multiqc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    conda:
        "envs/multiqc.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.html}) $(dirname {log})
        multiqc --force -o {OUTDIR}/qc/multiqc {OUTDIR} > {log} 2>&1
        """
