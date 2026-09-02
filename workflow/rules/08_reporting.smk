rule qc_gates:
    input:
        filtered_stats=expand(f"{OUTDIR}/stats/pairtools/{{sample}}.filtered.stats.txt", sample=SAMPLES),
        dedup_stats=expand(f"{OUTDIR}/stats/pairtools/{{sample}}.dedup.stats.txt", sample=SAMPLES),
        balance_markers=expand(f"{OUTDIR}/matrices/{{sample}}.mcool.balanced.done", sample=SAMPLES)
    output:
        tsv=f"{OUTDIR}/report/qc_gates.tsv"
    params:
        samples=SAMPLES,
        outdir=OUTDIR,
        gates=QC_GATES
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv})
        python - <<'PY' > {output.tsv}
import os

OUTDIR = {params.outdir!r}
SAMPLES = {params.samples!r}
gates = {params.gates!r}

def parse_stats(path):
    d = dict()
    if not os.path.isfile(path):
        return d
    for line in open(path):
        if "\t" in line:
            k, v = line.rstrip("\n").split("\t", 1)
            d[k] = v
    return d

print("sample\tmin_valid_pairs\tvalid_pairs\tpass_depth\tmin_frac_cis\tfrac_cis\tpass_cis\tmax_frac_dups\tfrac_dups\tpass_dups\tbalanced\tstatus")
for s in SAMPLES:
    fs = parse_stats(os.path.join(OUTDIR, "stats/pairtools", s + ".filtered.stats.txt"))
    ds = parse_stats(os.path.join(OUTDIR, "stats/pairtools", s + ".dedup.stats.txt"))
    total = int(fs.get("total", 0))
    cis = float(fs.get("summary/frac_cis", 0))
    dups = float(ds.get("summary/frac_dups", 1))
    marker = os.path.join(OUTDIR, "matrices", s + ".mcool.balanced.done")
    balanced = os.path.isfile(marker) and os.path.getsize(marker) > 0

    mv = int(gates.get("min_valid_pairs", 0))
    mc = float(gates.get("min_frac_cis", 0.0))
    md = float(gates.get("max_frac_dups", 1.0))
    p_depth = total >= mv
    p_cis = cis >= mc
    p_dups = dups <= md
    status = "pass" if (p_depth and p_cis and p_dups and balanced) else "failed/underpowered"
    print("\t".join(map(str, [s, mv, total, p_depth, mc, round(cis, 4), p_cis,
                              md, round(dups, 4), p_dups, balanced, status])))
PY
        """


rule filter_summary:
    input:
        filtered_stats=expand(f"{OUTDIR}/stats/pairtools/{{sample}}.filtered.stats.txt", sample=SAMPLES),
        dedup_stats=expand(f"{OUTDIR}/stats/pairtools/{{sample}}.dedup.stats.txt", sample=SAMPLES)
    output:
        tsv=f"{OUTDIR}/report/filter_summary.tsv"
    params:
        samples=SAMPLES,
        outdir=OUTDIR,
        info=lambda wc: {
            s: dict(
                protocol=sample_protocol(s),
                enzyme=str(_sample_cfg(s).get("enzyme", "")),
                base=sample_base_resolution(s),
                mapq=sample_min_mapq(s),
                maxcis=sample_max_cis_artifact_dist(s),
                samefrag=sample_same_fragment_max_dist(s),
                requnique=bool(_sample_cfg(s).get("require_unique", REQUIRE_UNIQUE)),
            )
            for s in SAMPLES
        }
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv})
        python - <<'PY' > {output.tsv}
import os

OUTDIR = {params.outdir!r}
SAMPLES = {params.samples!r}
info = {params.info!r}

def parse_stats(path):
    d = dict()
    if os.path.isfile(path):
        for line in open(path):
            if "\t" in line:
                k, v = line.rstrip("\n").split("\t", 1)
                d[k] = v
    return d

print("sample\tprotocol\tenzyme\tbase_resolution\tmin_mapq\tmax_cis_artifact_dist\tsame_fragment_max_dist\trequire_unique\tmapped\tdups\tnodups\tfrac_dups\tvalid_pairs\tfrac_cis\tfrac_cis_1kb+")
for s in SAMPLES:
    i = info[s]
    ds = parse_stats(os.path.join(OUTDIR, "stats/pairtools", s + ".dedup.stats.txt"))
    fs = parse_stats(os.path.join(OUTDIR, "stats/pairtools", s + ".filtered.stats.txt"))
    row = [
        s, i["protocol"], i["enzyme"], i["base"], i["mapq"], i["maxcis"],
        i["samefrag"], i["requnique"],
        int(ds.get("total_mapped", 0)),
        int(ds.get("total_dups", 0)),
        int(ds.get("total_nodups", 0)),
        round(float(ds.get("summary/frac_dups", 0)), 4),
        int(fs.get("total", 0)),
        round(float(fs.get("summary/frac_cis", 0)), 4),
        round(float(fs.get("summary/frac_cis_1kb+", 0)), 4),
    ]
    print("\t".join(map(str, row)))
PY
        """


rule run_manifest:
    input:
        validated=f"{OUTDIR}/metadata/sample_sheet.validated.tsv",
        filtered_stats=expand(f"{OUTDIR}/stats/pairtools/{{sample}}.filtered.stats.txt", sample=SAMPLES)
    output:
        tsv=f"{OUTDIR}/report/run_manifest.tsv",
        config_copy=f"{OUTDIR}/report/config.snapshot.yaml"
    params:
        samples=SAMPLES,
        outdir=OUTDIR,
        config_path="config.yaml",
        chrom_sizes=REF.get("chrom_sizes", ""),
        fasta=REF.get("bwa_indexed_fasta", ""),
        groups=GROUPS,
        info=lambda wc: {
            s: dict(
                protocol=sample_protocol(s),
                base=sample_base_resolution(s),
                maxcis=sample_max_cis_artifact_dist(s),
                samefrag=sample_same_fragment_max_dist(s),
            )
            for s in SAMPLES
        }
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {output.config_copy})
        cp {params.config_path} {output.config_copy}
        python - <<'PY' > {output.tsv}
import datetime, hashlib, os, subprocess

OUTDIR = {params.outdir!r}
SAMPLES = {params.samples!r}
config_path = {params.config_path!r}
chrom_sizes = {params.chrom_sizes!r}
fasta = {params.fasta!r}
groups = {params.groups!r}
info = {params.info!r}

def md5(path):
    # full-file streaming md5 (streams the whole file in 1MiB chunks, no
    # 64MiB cap; P0-2)
    try:
        h = hashlib.md5()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return "NA"

def git_rev():
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, timeout=10
        )
        if out.returncode == 0:
            return out.stdout.strip()
    except Exception:
        pass
    return "NA"

def git_dirty():
    try:
        out = subprocess.run(
            ["git", "status", "--porcelain"], capture_output=True, text=True, timeout=10
        )
        if out.returncode == 0:
            return "true" if out.stdout.strip() else "false"
    except Exception:
        pass
    return "NA"

def git_diff_hash():
    try:
        out = subprocess.run(["git", "diff"], capture_output=True, timeout=10)
        if out.returncode == 0:
            return hashlib.md5(out.stdout).hexdigest()
    except Exception:
        pass
    return "NA"

def size_or_na(path):
    try:
        return str(os.path.getsize(path))
    except OSError:
        return "NA"

print("key\tvalue")
print("timestamp\t" + datetime.datetime.now().isoformat())
print("git_commit\t" + git_rev())
print("git_dirty\t" + git_dirty())
print("git_diff_hash\t" + git_diff_hash())
print("config_md5\t" + md5(config_path))
print("chrom_sizes\t" + chrom_sizes)
print("chrom_sizes_md5\t" + md5(chrom_sizes))
print("reference_fasta\t" + fasta)
print("reference_fasta_md5\t" + md5(fasta))
print("n_samples\t" + str(len(SAMPLES)))
print("samples\t" + ",".join(SAMPLES))
print("groups\t" + str(groups))
for s in SAMPLES:
    i = info[s]
    print("sample/" + s + "/protocol\t" + i["protocol"])
    print("sample/" + s + "/base_resolution\t" + str(i["base"]))
    print("sample/" + s + "/max_cis_artifact_dist\t" + str(i["maxcis"]))
    print("sample/" + s + "/same_fragment_max_dist\t" + str(i["samefrag"]))

# ---- workflow hashes (P0-2): full-file md5 of every file under workflow/ ----
workflow_dir = "workflow"
if os.path.isdir(workflow_dir):
    for dirpath, dirnames, filenames in os.walk(workflow_dir):
        dirnames.sort()
        for fn in sorted(filenames):
            p = os.path.join(dirpath, fn)
            print(p + "\t" + md5(p))

# ---- balance summary (P1-2/P0-2): per resolution from balance.tsv ----
for s in SAMPLES:
    bpath = os.path.join(OUTDIR, "report/balance", s + ".balance.tsv")
    if not os.path.isfile(bpath):
        continue
    with open(bpath) as h:
        header = h.readline().rstrip("\n").split("\t")
        for line in h:
            row = line.rstrip("\n").split("\t")
            if len(row) < len(header):
                continue
            d = dict(zip(header, row))
            print("balance/" + s + "/" + str(d.get("resolution", "NA")) + "\t" + str(d.get("mode", "NA")) + "\t" + str(d.get("variance", "NA")) + "\t" + str(d.get("n_masked_bins", "NA")))

# ---- outputs manifest (P0-2): size in bytes of each key output ----
per_sample_outputs = [
    ("filtered.pairs.gz", "pairs/filtered", ".filtered.pairs.gz"),
    ("filtered.stats.txt", "stats/pairtools", ".filtered.stats.txt"),
    ("balanced.mcool", "matrices", ".balanced.mcool"),
    ("insulation.tsv", "features", "/insulation.tsv"),
    ("boundaries.bed", "features", "/boundaries.bed"),
    ("loops.bedpe", "features", "/loops.bedpe"),
    ("apa.tsv", "features", "/apa.tsv"),
]
for s in SAMPLES:
    for name, subdir, suffix in per_sample_outputs:
        p = os.path.join(OUTDIR, subdir, s + suffix)
        print("output/" + s + "/" + name + "\t" + size_or_na(p))
cross_sample_outputs = [
    ("qc_gates.tsv", "report/qc_gates.tsv"),
    ("filter_summary.tsv", "report/filter_summary.tsv"),
    ("run_manifest.tsv", "report/run_manifest.tsv"),
    ("integrative/matrix_compare.tsv", "integrative/matrix_compare.tsv"),
    ("integrative/compartment_flip.tsv", "integrative/compartment_flip.tsv"),
    ("integrative/differential_summary.tsv", "integrative/differential_summary.tsv"),
]
for name, rel in cross_sample_outputs:
    print("output/" + name + "\t" + size_or_na(os.path.join(OUTDIR, rel)))
PY
        """
rule multiqc_microc_metrics:
    input:
        fastp=expand(f"{OUTDIR}/qc/fastp/{{sample}}/fastp.json", sample=SAMPLES),
        dedup=expand(f"{OUTDIR}/stats/pairtools/{{sample}}.dedup.stats.txt", sample=SAMPLES),
        filtered=expand(f"{OUTDIR}/stats/pairtools/{{sample}}.filtered.stats.txt", sample=SAMPLES),
        boundaries=expand(f"{OUTDIR}/features/{{sample}}/boundaries.bed", sample=SAMPLES),
        loops=expand(f"{OUTDIR}/features/{{sample}}/loops.bedpe", sample=SAMPLES),
        apa=expand(f"{OUTDIR}/features/{{sample}}/apa.tsv", sample=SAMPLES),
        matrix=f"{OUTDIR}/qc/plots/replicate_concordance.tsv",
        boundary=f"{OUTDIR}/qc/plots/boundary_overlap.tsv",
        compartment=f"{OUTDIR}/qc/plots/compartment_correlation.tsv"
    output:
        library=f"{OUTDIR}/qc/multiqc_custom/microc_library_qc_mqc.yaml",
        matrix=f"{OUTDIR}/qc/multiqc_custom/microc_matrix_concordance_mqc.yaml",
        boundary=f"{OUTDIR}/qc/multiqc_custom/microc_boundary_overlap_mqc.yaml",
        compartment=f"{OUTDIR}/qc/multiqc_custom/microc_compartment_correlation_mqc.yaml"
    params:
        outdir=OUTDIR,
        samples=SAMPLES
    conda:
        "envs/multiqc.yaml"
    log:
        "logs/multiqc_microc_metrics.log"
    script:
        "../scripts/multiqc_microc_metrics.py"
