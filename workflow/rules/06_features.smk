rule call_compartments:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
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
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
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
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    output:
        bed=f"{OUTDIR}/features/{{sample}}/boundaries.bed",
        insulation=f"{OUTDIR}/features/{{sample}}/insulation.tsv"
    params:
        binsize=lambda wc: int(config.get("features", {}).get("feature_resolution", config.get("matrix", {}).get("base_resolution", 1000))),
        window=lambda wc: int(config.get("features", {}).get("insulation_window", 100000)),
        balance_max_iters=lambda wc: int(config.get("features", {}).get("balance_max_iters", 500)),
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
        cooler balance --max-iters {params.balance_max_iters} --convergence-policy store_final \
          {input.mcool}::/resolutions/{params.binsize} >> {log} 2>&1 || true
        # Restrict to contigs >= min_chrom_size: tiny chrUn contigs crash
        # cooltools insulation with an empty-bins IndexError.
        cooler dump -t chroms {input.mcool}::/resolutions/{params.binsize} \
          | awk 'BEGIN{{OFS="\t"}} $2>={params.min_chrom_size} {{print $1, 0, $2}}' > {output.insulation}.view.bed
        # cooltools 0.7.x CLI: WINDOW is a positional arg; --windows does not exist.
        cooltools insulation -p {threads} --view {output.insulation}.view.bed -o {output.insulation} \
          {input.mcool}::/resolutions/{params.binsize} {params.window} >> {log} 2>&1 || true
        python - <<'PY' >> {log} 2>&1
import pandas as pd
df = pd.read_csv("{output.insulation}", sep="\t")
bc = [c for c in df.columns if c.startswith("is_boundary_")]
if not bc:
    raise SystemExit("[WARN] no is_boundary_* column found in insulation table; writing empty boundaries.")
bound = df.loc[df[bc[0]].astype(bool), ["chrom", "start", "end"]].sort_values(["chrom", "start"])
bound.to_csv("{output.bed}", sep="\t", header=False, index=False)
print(f"Wrote {{len(bound)}} boundaries")
PY
        if [[ ! -s {output.bed} ]]; then
          : > {output.bed}
        fi
        """


rule call_loops:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool",
        expected=f"{OUTDIR}/features/{{sample}}/expected_cis_{FEATURE_RESOLUTION}bp.tsv"
    output:
        bedpe=f"{OUTDIR}/features/{{sample}}/loops.bedpe"
    params:
        binsize=FEATURE_RESOLUTION,
        balance_max_iters=lambda wc: int(config.get("features", {}).get("balance_max_iters", 500)),
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
        cooler balance --max-iters {params.balance_max_iters} --convergence-policy store_final \
          {input.mcool}::/resolutions/{params.binsize} >> {log} 2>&1 || true
        # cooltools 0.7.x CLI: dots requires EXPECTED_PATH positional and -o output.
        # --view must match the view used to compute expected (expected only covers
        # contigs >= min_chrom_size; without --view dots builds a view of ALL
        # contigs and the expected/view region-name check fails).
        cooler dump -t chroms {input.mcool}::/resolutions/{params.binsize} \
          | awk 'BEGIN{{OFS="\t"}} $2>={params.min_chrom_size} {{print $1, 0, $2}}' > {output.bedpe}.view.bed
        cooltools dots -p {threads} --view {output.bedpe}.view.bed -o {output.bedpe} \
          {input.mcool}::/resolutions/{params.binsize} {input.expected} >> {log} 2>&1 || true
        if [[ ! -s {output.bedpe} ]]; then
          : > {output.bedpe}
          echo "[WARN] cooltools dots failed or produced empty output; wrote empty loops file." >> {log}
        fi
        """


rule pileup_apa:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool",
        loops=f"{OUTDIR}/features/{{sample}}/loops.bedpe"
    output:
        png=f"{OUTDIR}/features/{{sample}}/apa.png",
        tsv=f"{OUTDIR}/features/{{sample}}/apa.tsv"
    params:
        binsize=FEATURE_RESOLUTION,
        flank=APA_FLANK
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
import numpy as np
import pandas as pd

uri = "{input.mcool}::/resolutions/{params.binsize}"
flank = {params.flank}

if os.path.getsize("{input.loops}") == 0:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.text(0.5, 0.5, "no loops called", ha="center", va="center")
    ax.axis("off")
    fig.savefig("{output.png}", dpi=200)
    with open("{output.tsv}", "w") as h:
        h.write("n_loops\tmean\tmax\n0\tNA\tNA\n")
    raise SystemExit("[WARN] loops.bedpe is empty; wrote placeholder APA")
import cooler
from coolpuppy import coolpup

clr = cooler.Cooler(uri)
loops = pd.read_csv("{input.loops}", sep="\t")
for c in ["start1", "end1", "start2", "end2"]:
    loops[c] = loops[c].astype(int)

# coolpuppy pileup requires balancing weights (zoomify --balance provides
# them). nproc=1 avoids a numpy2/coolpuppy multiprocessing pickling issue.
pups = coolpup.pileup(clr, loops, features_format="bedpe", flank=flank, nproc=1, nshifts=0)
mat = pups["data"].iloc[0]
pup = np.asarray(mat, dtype=float)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(4, 4))
im = ax.imshow(pup, cmap="Reds", interpolation="nearest")
ax.set_title("APA (n={{0}})".format(len(loops)))
plt.colorbar(im, ax=ax, fraction=0.046)
fig.savefig("{output.png}", dpi=200)

with open("{output.tsv}", "w") as h:
    h.write("n_loops\tmean\tmax\n")
    h.write(str(len(loops)) + "\t" + repr(float(np.nanmean(pup))) + "\t" + repr(float(np.nanmax(pup))) + "\n")
print("APA done: n_loops=" + str(len(loops)) + ", mean=" + str(np.nanmean(pup)))
PY
        """
