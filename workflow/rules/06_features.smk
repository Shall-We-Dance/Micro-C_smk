rule call_compartments:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.balanced.mcool",
        balanced=f"{OUTDIR}/matrices/{{sample}}.mcool.balanced.done"
    output:
        bedgraph=f"{OUTDIR}/features/{{sample}}/compartments_{{res}}bp.bedgraph"
    params:
        binsize=lambda wc: int(wc.res),
        fasta=lambda wc: config.get("reference", {}).get("bwa_indexed_fasta", ""),
        min_chrom_size=MIN_CHROM_SIZE
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/features/compartments/{{sample}}.{{res}}bp.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bedgraph}) $(dirname {log})
        python - <<'PY' > {log} 2>&1
import cooler
import cooltools
import bioframe
import pandas as pd

if not "{params.fasta}":
    raise ValueError("reference.bwa_indexed_fasta is required for compartment GC phasing")

clr = cooler.Cooler("{input.mcool}::/resolutions/{params.binsize}")

# Restrict to contigs >= min_chrom_size: tiny chrUn contigs produce noisy
# eigendecompositions and mostly-NaN tracks (dm6 has ~1860 of them).
chroms = [
    c for c, s in zip(clr.chromnames, clr.chromsizes)
    if s >= {params.min_chrom_size}
]
bins = clr.bins()[:][["chrom", "start", "end"]]
bins = bins[bins["chrom"].isin(chroms)].reset_index(drop=True)

genome = bioframe.load_fasta("{params.fasta}")
gc_cov = bioframe.frac_gc(bins, genome)

view_df = pd.DataFrame(
    {{
        "chrom": chroms,
        "start": 0,
        "end": clr.chromsizes[chroms].values,
        "name": chroms,
    }}
)

cis_eigs = cooltools.eigs_cis(
    clr,
    gc_cov,
    view_df=view_df,
    n_eigs=3,
)
eigenvector_track = cis_eigs[1][["chrom", "start", "end", "E1"]]
eigenvector_track = eigenvector_track.sort_values(["chrom", "start"])
eigenvector_track.to_csv("{output.bedgraph}", sep="\t", header=False, index=False)
print(f"Saved to: {output.bedgraph}")
PY
        """


rule compute_expected_cis:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.balanced.mcool",
        balanced=f"{OUTDIR}/matrices/{{sample}}.mcool.balanced.done"
    output:
        expected=f"{OUTDIR}/features/{{sample}}/expected_cis_{{res}}bp.tsv"
    params:
        binsize=lambda wc: int(wc.res),
        min_chrom_size=MIN_CHROM_SIZE
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/features/expected_cis/{{sample}}.{{res}}bp.log"
    threads: 8
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.expected}) $(dirname {log})
        # Only use contigs >= min_chrom_size: tiny chrUn contigs produce empty
        # bins that crash cooltools expected-cis/insulation (empty-index error).
        cooler dump -t chroms {input.mcool}::/resolutions/{params.binsize} \
          | awk 'BEGIN{{OFS="\t"}} $2>={params.min_chrom_size} {{print $1, 0, $2}}' > {output.expected}.view.bed
        cooltools expected-cis -p {threads} --view {output.expected}.view.bed \
          -o {output.expected} {input.mcool}::/resolutions/{params.binsize} > {log} 2>&1
        """


rule call_insulation_boundaries:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.balanced.mcool",
        balanced=f"{OUTDIR}/matrices/{{sample}}.mcool.balanced.done"
    output:
        bed=f"{OUTDIR}/features/{{sample}}/boundaries.bed",
        insulation=f"{OUTDIR}/features/{{sample}}/insulation.tsv"
    params:
        binsize=lambda wc: int(config.get("features", {}).get("feature_resolution", config.get("matrix", {}).get("base_resolution", 1000))),
        window=lambda wc: int(config.get("features", {}).get("insulation_window", 100000)),
        min_chrom_size=MIN_CHROM_SIZE
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/features/boundaries/{{sample}}.log"
    threads: 8
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bed}) $(dirname {output.insulation}) $(dirname {log})
        # Restrict to contigs >= min_chrom_size: tiny chrUn contigs crash
        # cooltools insulation with an empty-bins IndexError.
        cooler dump -t chroms {input.mcool}::/resolutions/{params.binsize} \
          | awk 'BEGIN{{OFS="\t"}} $2>={params.min_chrom_size} {{print $1, 0, $2}}' > {output.insulation}.view.bed
        # cooltools 0.7.x CLI: WINDOW is a positional arg; --windows does not exist.
        cooltools insulation -p {threads} --view {output.insulation}.view.bed -o {output.insulation} \
          {input.mcool}::/resolutions/{params.binsize} {params.window} >> {log} 2>&1
        python - <<'PY' >> {log} 2>&1
import os
import sys
import pandas as pd

ins_path = "{output.insulation}"
bed_path = "{output.bed}"

# Missing/empty insulation table or a missing is_boundary_* column means the
# insulation computation itself failed (schema absent): fail the rule instead
# of silently writing an empty boundaries file.
if not os.path.isfile(ins_path) or os.path.getsize(ins_path) == 0:
    print("[ERROR] insulation output missing/column absent - likely computation failure", file=sys.stderr)
    sys.exit(1)

df = pd.read_csv(ins_path, sep="\t")
bc = [c for c in df.columns if c.startswith("is_boundary_")]
if not bc:
    print("[ERROR] insulation output missing/column absent - likely computation failure", file=sys.stderr)
    sys.exit(1)

# Column exists but 0 rows are True: legitimately empty result (valid_empty).
bound = df.loc[df[bc[0]].astype(bool), ["chrom", "start", "end"]].sort_values(["chrom", "start"])
if len(bound) == 0:
    print("[WARN] boundary column present but 0 boundary rows; writing empty boundaries (valid_empty).")
    open(bed_path, "w").close()
    sys.exit(0)

bound.to_csv(bed_path, sep="\t", header=False, index=False)
print("Wrote " + str(len(bound)) + " boundaries")
PY
        if [[ ! -s {output.bed} ]]; then
          : > {output.bed}
        fi
        """


rule call_loops:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.balanced.mcool",
        expected=f"{OUTDIR}/features/{{sample}}/expected_cis_{FEATURE_RESOLUTION}bp.tsv",
        balanced=f"{OUTDIR}/matrices/{{sample}}.mcool.balanced.done"
    output:
        bedpe=f"{OUTDIR}/features/{{sample}}/loops.bedpe"
    params:
        binsize=FEATURE_RESOLUTION,
        min_chrom_size=MIN_CHROM_SIZE
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/features/loops/{{sample}}.log"
    threads: 8
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bedpe}) $(dirname {log})
        # cooltools 0.7.x CLI: dots requires EXPECTED_PATH positional and -o output.
        # --view must match the view used to compute expected (expected only covers
        # contigs >= min_chrom_size; without --view dots builds a view of ALL
        # contigs and the expected/view region-name check fails).
        cooler dump -t chroms {input.mcool}::/resolutions/{params.binsize} \
          | awk 'BEGIN{{OFS="\t"}} $2>={params.min_chrom_size} {{print $1, 0, $2}}' > {output.bedpe}.view.bed
        # dots must succeed; a non-zero exit fails the rule (no `|| true`).
        cooltools dots -p {threads} --view {output.bedpe}.view.bed -o {output.bedpe} \
          {input.mcool}::/resolutions/{params.binsize} {input.expected} >> {log} 2>&1
        # dots exited 0 but produced nothing -> legitimately empty result.
        if [[ ! -s {output.bedpe} ]]; then
          : > {output.bedpe}
          echo "[WARN] cooltools dots produced empty output; wrote empty loops file." >> {log}
        fi
        """


rule pileup_apa:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.balanced.mcool",
        loops=f"{OUTDIR}/features/{{sample}}/loops.bedpe",
        balanced=f"{OUTDIR}/matrices/{{sample}}.mcool.balanced.done"
    output:
        png=f"{OUTDIR}/features/{{sample}}/apa.png",
        tsv=f"{OUTDIR}/features/{{sample}}/apa.tsv"
    params:
        binsize=FEATURE_RESOLUTION,
        flank=APA_FLANK,
        nshifts=APA_NSHIFTS
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/features/apa/{{sample}}.log"
    threads: 8
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.png}) $(dirname {output.tsv}) $(dirname {log})
        python - <<'PY' > {log} 2>&1
import os
import sys
import math
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

uri = "{input.mcool}::/resolutions/{params.binsize}"
flank = {params.flank}
nshifts = {params.nshifts}


def fmt(x):
    return "NA" if not np.isfinite(x) else "{{:.6g}}".format(float(x))


def write_placeholder():
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.text(0.5, 0.5, "no loops called", ha="center", va="center")
    ax.axis("off")
    fig.savefig("{output.png}", dpi=200)
    with open("{output.tsv}", "w") as h:
        h.write("n_loops\tcenter_mean\tcorner_mean\tratio\tnshifts\n")
        h.write("0\tNA\tNA\tNA\t" + str(nshifts) + "\n")
    print("[WARN] loops.bedpe is empty; wrote placeholder APA")
    sys.exit(0)


# Legitimately empty result: no loops called -> placeholder APA + exit 0.
if os.path.getsize("{input.loops}") == 0:
    write_placeholder()

loops = pd.read_csv("{input.loops}", sep="\t")
if loops.shape[0] == 0:
    write_placeholder()

import cooler
from coolpuppy import coolpup

clr = cooler.Cooler(uri)
for c in ["start1", "end1", "start2", "end2"]:
    loops[c] = loops[c].astype(int)

# coolpuppy 1.1.0 has no pileup_with_controls(); the module-level pileup()
# enables random shifted-anchor controls via nshifts > 0 (it then uses the
# internal pileupsWithControl), so each pileup is normalized by its shifted
# control (O/E-style background subtraction). seed=42 makes the shifts
# reproducible. nproc=1 avoids a numpy2/coolpuppy multiprocessing pickling
# issue.
pups = coolpup.pileup(clr, loops, features_format="bedpe", flank=flank,
                      mindist=0, nproc=1, nshifts=nshifts, seed=42)
mat = pups["data"].iloc[0]
pup = np.asarray(mat, dtype=float)

# Center-vs-background APA enrichment:
#  - center_mean: mean of the central 25% of the pileup area
#    (square of side n//2 centered on the loop anchors).
#  - corner_mean: mean of the four corner squares, each ~12.5% of the area
#    (side sqrt(0.125)*n), clamped to n//4 so they stay outside the central
#    block; corner pixels are far from the anchors and represent background.
#  - ratio: center_mean / corner_mean.
n = pup.shape[0]
center_side = n // 2
lo = (n - center_side) // 2
hi = lo + center_side
center = pup[lo:hi, lo:hi]

corner_side = max(1, min(int(round(n * math.sqrt(0.125))), n // 4))
c = corner_side
corners = [
    pup[:c, :c],
    pup[:c, n - c:],
    pup[n - c:, :c],
    pup[n - c:, n - c:],
]
center_mean = float(np.nanmean(center))
corner_mean = float(np.nanmean(np.concatenate([cor.ravel() for cor in corners])))
if np.isfinite(center_mean) and np.isfinite(corner_mean) and corner_mean != 0:
    ratio = center_mean / corner_mean
else:
    ratio = float("nan")

fig, ax = plt.subplots(figsize=(4, 4))
im = ax.imshow(pup, cmap="Reds", interpolation="nearest")
ax.set_title("APA (n={{0}})".format(len(loops)))
plt.colorbar(im, ax=ax, fraction=0.046)
fig.savefig("{output.png}", dpi=200)

with open("{output.tsv}", "w") as h:
    h.write("n_loops\tcenter_mean\tcorner_mean\tratio\tnshifts\n")
    h.write(str(len(loops)) + "\t" + fmt(center_mean) + "\t" + fmt(corner_mean) + "\t" + fmt(ratio) + "\t" + str(nshifts) + "\n")
print("APA done: n_loops=" + str(len(loops)) + ", center_mean=" + fmt(center_mean) + ", corner_mean=" + fmt(corner_mean) + ", ratio=" + fmt(ratio))
PY
        """
