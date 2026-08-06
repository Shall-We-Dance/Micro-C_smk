rule qc_cis_trans:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    output:
        tsv=f"{OUTDIR}/qc/plots/{{sample}}/cis_trans.tsv"
    params:
        min_chrom_size=MIN_CHROM_SIZE
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/qc/cis_trans/{{sample}}.log"
    threads: 8
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        RES={config[matrix][base_resolution]}
        # A --view BED needs chrom,start,end (cooler dump -t chroms gives 2 cols).
        # Also drop contigs < min_chrom_size (tiny chrUn contigs crash cooltools).
        cooler dump -t chroms {input.mcool}::/resolutions/$RES | awk 'BEGIN{{OFS="\t"}} $2>={params.min_chrom_size} {{print $1, 0, $2}}' > {output.tsv}.view.bed
        cooltools expected-cis -p {threads} --view {output.tsv}.view.bed \
          -o {output.tsv} {input.mcool}::/resolutions/$RES > {log} 2>&1
        """


rule qc_distance_decay:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    output:
        png=f"{OUTDIR}/qc/plots/{{sample}}/distance_decay.png",
        expected=f"{OUTDIR}/qc/plots/{{sample}}/expected_cis.tsv"
    params:
        min_chrom_size=MIN_CHROM_SIZE
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/qc/distance_decay/{{sample}}.log"
    threads: 8
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.png}) $(dirname {output.expected}) $(dirname {log})
        RES={config[matrix][base_resolution]}
        cooler dump -t chroms {input.mcool}::/resolutions/$RES | awk 'BEGIN{{OFS="\t"}} $2>={params.min_chrom_size} {{print $1, 0, $2}}' > {output.expected}.view.bed
        cooltools expected-cis -p {threads} --view {output.expected}.view.bed \
          -o {output.expected} {input.mcool}::/resolutions/$RES > {log} 2>&1
        python - <<'PY'
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("{output.expected}", sep="\t")
plt.figure(figsize=(5,4))
plt.loglog(df['dist_bp'], df['balanced.avg'], lw=1)
plt.xlabel('distance (bp)')
plt.ylabel('contact frequency (balanced)')
plt.tight_layout()
plt.savefig('{output.png}', dpi=200)
print('distance_decay plot saved')
PY
        """


rule qc_replicate_concordance:
    input:
        mcools=expand(f"{OUTDIR}/matrices/{{sample}}.mcool", sample=SAMPLES)
    output:
        tsv=f"{OUTDIR}/qc/plots/replicate_concordance.tsv"
    params:
        resolutions=CONCORDANCE_RESOLUTIONS,
        min_chrom_size=MIN_CHROM_SIZE,
        samples=SAMPLES,
        outdir=OUTDIR
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/qc/concordance.log"
    threads: 8
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python - <<'PY' > {log} 2>&1
import os
import itertools
import cooler
import numpy as np
from scipy import stats

OUTDIR = {params.outdir!r}
SAMPLES = {params.samples!r}
params_resolutions = {params.resolutions!r}
min_chrom_size = {params.min_chrom_size}

rows = []
avail_per_sample = []
for s in SAMPLES:
    mc = OUTDIR + "/matrices/" + s + ".mcool"
    avail_per_sample.append({{int(r.rsplit("/", 1)[-1]) for r in cooler.fileops.list_coolers(mc)}})
for res in params_resolutions:
    if not all(res in av for av in avail_per_sample):
        for a, b in itertools.combinations(SAMPLES, 2):
            rows.append([a, b, res, "NA", "NA", 0])
        continue
    for a, b in itertools.combinations(SAMPLES, 2):
        clr_a = cooler.Cooler(OUTDIR + "/matrices/" + a + ".mcool::/resolutions/" + str(res))
        clr_b = cooler.Cooler(OUTDIR + "/matrices/" + b + ".mcool::/resolutions/" + str(res))
        chroms = [
            c
            for c, s in zip(clr_a.chromnames, clr_a.chromsizes)
            if s >= min_chrom_size
        ]
        vals_a, vals_b = [], []
        for chrom in chroms:
            ma = clr_a.matrix(balance=True).fetch(chrom)
            mb = clr_b.matrix(balance=True).fetch(chrom)
            tri_a = ma[np.triu_indices_from(ma, k=1)]
            tri_b = mb[np.triu_indices_from(mb, k=1)]
            finite = np.isfinite(tri_a) & np.isfinite(tri_b)
            vals_a.append(tri_a[finite])
            vals_b.append(tri_b[finite])
        if not vals_a:
            rows.append([a, b, res, "NA", "NA", 0])
            continue
        va = np.concatenate(vals_a)
        vb = np.concatenate(vals_b)
        pearson = stats.pearsonr(va, vb).statistic
        spearman = stats.spearmanr(va, vb).statistic
        rows.append([a, b, res, round(float(pearson), 4), round(float(spearman), 4), int(len(va))])

with open("{output.tsv}", "w") as handle:
    handle.write("sample1\tsample2\tresolution\tpearson\tspearman\tn_bins\n")
    for row in sorted(rows, key=lambda r: (r[1], r[0], r[2])):
        handle.write("\t".join(map(str, row)) + "\n")
print("computed concordance for " + str(len(rows)) + " sample-pair/resolution combinations")
PY
        """


rule qc_matrix_snapshot:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    output:
        png=f"{OUTDIR}/qc/plots/{{sample}}/matrix_snapshot.png"
    params:
        min_chrom_size=MIN_CHROM_SIZE
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/qc/matrix_snapshot/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.png}) $(dirname {log})
        python - <<'PY'
import cooler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

uri = '{input.mcool}::/resolutions/{config[matrix][base_resolution]}'
clr = cooler.Cooler(uri)
chrom = next(
    c for c, s in zip(clr.chromnames, clr.chromsizes)
    if s >= {params.min_chrom_size}
)
mat = clr.matrix(balance=True).fetch(chrom)
# Center the window on the chromosome: telomeric/centromeric edges are often
# unbalanced (NaN), producing all-NaN slices with a fixed vmax.
n = mat.shape[0]
mid = n // 2
win = slice(mid - 150, mid + 150)
snap = mat[win, win]
vmax = float(np.nanpercentile(snap, 98)) if np.isfinite(snap).any() else 1.0
plt.figure(figsize=(5,5))
plt.imshow(snap, cmap='Reds', vmax=vmax)
plt.title('{wildcards.sample} ' + chrom)
plt.tight_layout()
plt.savefig('{output.png}', dpi=200)
PY
        """
