rule qc_cis_trans:
    input:
        stats=f"{OUTDIR}/stats/pairtools/{{sample}}.filtered.stats.txt"
    output:
        tsv=f"{OUTDIR}/qc/plots/{{sample}}/cis_trans.tsv"
    log:
        f"logs/qc/cis_trans/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        # Real cis/trans summary (the expected-cis table lives in
        # qc_distance_decay's expected_cis.tsv).
        python - <<'PY' > {log} 2>&1
d = {{}}
for line in open("{input.stats}"):
    if "\\t" in line:
        k, v = line.rstrip("\\n").split("\\t", 1)
        d[k] = v
tot = int(d.get("total", 0))
cis = int(d.get("cis", 0))
trans = int(d.get("trans", 0))
c1k = int(d.get("cis_1kb+", 0))
c10k = int(d.get("cis_10kb+", 0))
rows = [
    ["metric", "value"],
    ["total_pairs", tot],
    ["cis", cis],
    ["trans", trans],
    ["frac_cis", round(cis / tot, 6) if tot else "NA"],
    ["frac_cis_1kb+", round(c1k / tot, 6) if tot else "NA"],
    ["frac_cis_10kb+", round(c10k / tot, 6) if tot else "NA"],
]
# per-chromosome cis/trans
for chrom in sorted(set(k.split("/")[1] for k in d if k.startswith("chrom_freq/") and len(k.split("/")) == 3)):
    c = int(d.get(f"chrom_freq/{{chrom}}/{{chrom}}", 0))
    all_c = sum(int(d.get(f"chrom_freq/{{chrom}}/{{x}}", 0)) for x in set(k.split("/")[2] for k in d if k.startswith(f"chrom_freq/{{chrom}}/")))
    rows.append([f"chr_cis/{{chrom}}", c])
    rows.append([f"chr_cis_frac/{{chrom}}", round(c / all_c, 6) if all_c else "NA"])
with open("{output.tsv}", "w") as h:
    for row in rows:
        h.write("\\t".join(map(str, row)) + "\\n")
print("cis/trans summary written")
PY
        """


rule qc_distance_decay:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    output:
        png=f"{OUTDIR}/qc/plots/{{sample}}/distance_decay.png",
        expected=f"{OUTDIR}/qc/plots/{{sample}}/expected_cis.tsv"
    params:
        min_chrom_size=MIN_CHROM_SIZE,
        base=lambda wc: sample_base_resolution(wc.sample)
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/qc/distance_decay/{{sample}}.log"
    threads: 8
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.png}) $(dirname {output.expected}) $(dirname {log})
        RES={params.base}
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
        min_chrom_size=MIN_CHROM_SIZE,
        base=lambda wc: sample_base_resolution(wc.sample)
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

uri = '{input.mcool}::/resolutions/{params.base}'
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


rule qc_boundary_overlap:
    input:
        boundaries=expand(f"{OUTDIR}/features/{{sample}}/boundaries.bed", sample=SAMPLES)
    output:
        tsv=f"{OUTDIR}/qc/plots/boundary_overlap.tsv"
    params:
        samples=SAMPLES,
        outdir=OUTDIR,
        # count as overlapping if boundary midpoints are within this distance
        window=100000
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/qc/boundary_overlap.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python - <<'PY' > {log} 2>&1
import itertools, os
import numpy as np

OUTDIR = {params.outdir!r}
SAMPLES = {params.samples!r}
window = {params.window}

def load_bounds(s):
    path = OUTDIR + "/features/" + s + "/boundaries.bed"
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return {{}}
    out = {{}}
    for line in open(path):
        p = line.rstrip("\\n").split("\\t")
        chrom = p[0]
        mid = (int(p[1]) + int(p[2])) // 2
        out.setdefault(chrom, []).append(mid)
    for c in out:
        out[c] = sorted(out[c])
    return out

def jaccard(a, b):
    # fraction of a's boundaries (per chrom) within +/-window of a b boundary
    hits = 0
    for chrom, mids in a.items():
        if chrom not in b:
            continue
        bm = b[chrom]
        j = 0
        for m in mids:
            lo, hi = m - window, m + window
            i = np.searchsorted(bm, lo)
            if i < len(bm) and bm[i] <= hi:
                hits += 1
    return hits / max(1, sum(len(v) for v in a.values()))

data = {{s: load_bounds(s) for s in SAMPLES}}
with open("{output.tsv}", "w") as h:
    h.write("sample1\\tsample2\\tn1\\tn2\\tjaccard_s1_over_s2\\tjaccard_s2_over_s1\\n_overlap\\n_union\\n_jaccard\\n")
    for a, b in itertools.combinations(SAMPLES, 2):
        n1 = sum(len(v) for v in data[a].values())
        n2 = sum(len(v) for v in data[b].values())
        j12 = jaccard(data[a], data[b])
        j21 = jaccard(data[b], data[a])
        over = round(j12 * n1)
        uni = n1 + n2 - over
        jac = round(over / max(1, uni), 4)
        h.write("\\t".join(map(str, [a, b, n1, n2, round(j12,4), round(j21,4), over, uni, jac])) + "\\n")
print("boundary overlap done for", len(SAMPLES), "samples")
PY
        """


rule qc_compartment_correlation:
    input:
        compartments=expand(
            f"{OUTDIR}/features/{{sample}}/compartments_{{res}}bp.bedgraph",
            sample=SAMPLES,
            res=COMPARTMENT_RESOLUTIONS,
        )
    output:
        tsv=f"{OUTDIR}/qc/plots/compartment_correlation.tsv"
    params:
        samples=SAMPLES,
        outdir=OUTDIR,
        resolutions=COMPARTMENT_RESOLUTIONS
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/qc/compartment_correlation.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python - <<'PY' > {log} 2>&1
import itertools, os
import numpy as np
import pandas as pd
from scipy import stats

OUTDIR = {params.outdir!r}
SAMPLES = {params.samples!r}
resolutions = {params.resolutions!r}

def load_e1(s, res):
    path = OUTDIR + "/features/" + s + "/compartments_" + str(res) + "bp.bedgraph"
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return None
    df = pd.read_csv(path, sep="\\t", header=None, names=["chrom", "start", "end", "E1"])
    df = df.dropna(subset=["E1"])
    return df.set_index(["chrom", "start"])["E1"]

with open("{output.tsv}", "w") as h:
    h.write("sample1\\tsample2\\tresolution\\tchrom\\tpearson\\tspearman\\tn_bins\\n")
    for res in resolutions:
        tracks = {{s: load_e1(s, res) for s in SAMPLES}}
        for a, b in itertools.combinations(SAMPLES, 2):
            ta, tb = tracks[a], tracks[b]
            if ta is None or tb is None:
                continue
            both = pd.concat([ta, tb], axis=1, keys=["a", "b"]).dropna()
            for chrom, grp in both.groupby(level="chrom"):
                if len(grp) < 10:
                    continue
                va = grp["a"].values
                vb = grp["b"].values
                pr = stats.pearsonr(va, vb).statistic
                sr = stats.spearmanr(va, vb).statistic
                h.write("\\t".join(map(str, [a, b, res, chrom, round(float(pr),4), round(float(sr),4), len(grp)])) + "\\n")
print("compartment correlation done")
PY
        """
