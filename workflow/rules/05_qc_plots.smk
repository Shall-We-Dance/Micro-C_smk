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
    if "\t" in line:
        k, v = line.rstrip("\n").split("\t", 1)
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
# per-chromosome cis/trans: chrom_freq keys are a->b pair counts and both
# orders may appear, so a chromosome can occur as chrom1 or chrom2. Total
# counts all unique keys where the chrom is either side (no double-counting).
chrom_pairs = set(
    tuple(k.split("/")[1:3])
    for k in d
    if k.startswith("chrom_freq/") and len(k.split("/")) == 3
)
for chrom in sorted({{c for a, b in chrom_pairs for c in (a, b)}}):
    cis_c = int(d.get(f"chrom_freq/{{chrom}}/{{chrom}}", 0))
    total_c = sum(
        int(d[f"chrom_freq/{{a}}/{{b}}"])
        for a, b in chrom_pairs
        if a == chrom or b == chrom
    )
    rows.append([f"chr_cis/{{chrom}}", cis_c])
    rows.append([f"chr_cis_frac/{{chrom}}", round(cis_c / total_c, 6) if total_c else "NA"])
with open("{output.tsv}", "w") as h:
    for row in rows:
        h.write("\t".join(map(str, row)) + "\n")
print("cis/trans summary written")
PY
        """


rule qc_distance_decay:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.balanced.mcool",
        balanced=f"{OUTDIR}/matrices/{{sample}}.mcool.balanced.done"
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
        mcools=expand(f"{OUTDIR}/matrices/{{sample}}.balanced.mcool", sample=SAMPLES),
        balanced=expand(f"{OUTDIR}/matrices/{{sample}}.mcool.balanced.done", sample=SAMPLES)
    output:
        tsv=f"{OUTDIR}/qc/plots/replicate_concordance.tsv"
    params:
        resolutions=CONCORDANCE_RESOLUTIONS,
        strata=CONCORDANCE_STRATA,
        min_chrom_size=MIN_CHROM_SIZE,
        samples=SAMPLES,
        outdir=OUTDIR,
        sample_meta={s: _sample_cfg(s).get("biological_sample", "NA") for s in SAMPLES}
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
import itertools
import os
import sys
import warnings

import cooler
import numpy as np
from scipy import stats

warnings.filterwarnings("ignore")

OUTDIR = {params.outdir!r}
SAMPLES = {params.samples!r}
params_resolutions = {params.resolutions!r}
strata = {params.strata!r}
min_chrom_size = {params.min_chrom_size}
sample_meta = {params.sample_meta!r}

rng = np.random.default_rng(1234)
rows = []

avail_per_sample = []
for s in SAMPLES:
    mc = OUTDIR + "/matrices/" + s + ".balanced.mcool"
    avail_per_sample.append(
        {{int(r.rsplit("/", 1)[-1]) for r in cooler.fileops.list_coolers(mc)}}
    )

def rep_of(s):
    v = sample_meta.get(s) or "NA"
    return str(v)

def schema_fail(msg):
    sys.stderr.write("schema check failed for replicate_concordance.tsv: " + msg + "\n")
    raise SystemExit(1)

for res in params_resolutions:
    # A stratum [lo, hi) is only reachable at resolution R if it spans bins
    # of that resolution: open-ended strata (hi == 0) are always reachable,
    # closed strata only when hi > R. Unreachable strata are skipped entirely.
    reachable = [(lo, hi) for lo, hi in strata if hi == 0 or hi > res]
    available = all(res in av for av in avail_per_sample)
    for a, b in itertools.combinations(SAMPLES, 2):
        is_rep = (
            "true"
            if (a != b and rep_of(a) != "NA" and rep_of(b) != "NA"
                and rep_of(a) == rep_of(b))
            else "false"
        )
        if not available:
            for lo, hi in reachable:
                rows.append([a, b, res, lo, hi, "oe", "NA", "NA", 0, is_rep])
            continue
        clr_a = cooler.Cooler(OUTDIR + "/matrices/" + a + ".balanced.mcool::/resolutions/" + str(res))
        clr_b = cooler.Cooler(OUTDIR + "/matrices/" + b + ".balanced.mcool::/resolutions/" + str(res))
        chroms_a = {{c for c, s in zip(clr_a.chromnames, clr_a.chromsizes) if s >= min_chrom_size}}
        chroms_b = {{c for c, s in zip(clr_b.chromnames, clr_b.chromsizes) if s >= min_chrom_size}}
        va_all, vb_all, dist_all = [], [], []
        fa_all, fb_all = [], []
        for chrom in sorted(chroms_a & chroms_b):
            ma = clr_a.matrix(balance=True).fetch(chrom)
            mb = clr_b.matrix(balance=True).fetch(chrom)
            iu, ju = np.triu_indices_from(ma, k=1)
            ta = ma[iu, ju]
            tb = mb[iu, ju]
            fa_all.append(np.isfinite(ta))
            fb_all.append(np.isfinite(tb))
            va_all.append(ta)
            vb_all.append(tb)
            dist_all.append((ju - iu) * res)
        if not va_all:
            for lo, hi in reachable:
                rows.append([a, b, res, lo, hi, "oe", "NA", "NA", 0, is_rep])
            continue
        va_all = np.concatenate(va_all)
        vb_all = np.concatenate(vb_all)
        fa_all = np.concatenate(fa_all)
        fb_all = np.concatenate(fb_all)
        dist_all = np.concatenate(dist_all)
        # per-stratum shared finite pixels plus per-sample finite totals;
        # the finite masks may differ per sample even on shared bins, so the
        # totals are counted per sample.
        stratum_shared = []
        total_a = 0
        total_b = 0
        for lo, hi in reachable:
            # stratum = distance in [lo, hi); hi == 0 means >= lo (open-ended)
            mask = dist_all >= lo
            if hi > 0:
                mask &= dist_all < hi
            ta_s = va_all[mask]
            tb_s = vb_all[mask]
            fa_s = fa_all[mask]
            fb_s = fb_all[mask]
            finite_a_s = int(fa_s.sum())
            finite_b_s = int(fb_s.sum())
            total_a += finite_a_s
            total_b += finite_b_s
            stratum_shared.append((lo, hi, ta_s, tb_s, fa_s, fb_s))
        # matched depth: K = min(total_a, total_b); each stratum contributes
        # n_s = round(K * count_s / max(total_a, total_b)) pixels drawn from
        # its shared finite pixels. When the totals are equal, K == max and
        # n_s == count_s, i.e. no subsampling happens.
        K = min(total_a, total_b)
        M = max(total_a, total_b)
        for lo, hi, ta_s, tb_s, fa_s, fb_s in stratum_shared:
            shared = fa_s & fb_s
            count_s = int(shared.sum())
            n_s = int(round(K * count_s / M)) if M > 0 else 0
            if count_s == 0 or n_s < 50:
                rows.append([a, b, res, lo, hi, "oe", "NA", "NA", count_s, is_rep])
                continue
            # O/E within the stratum: each sample is scaled by its own stratum
            # mean (over its own finite pixels) before correlation, removing
            # the within-stratum P(s) gradient.
            oe_a = ta_s[shared] / ta_s[fa_s].mean()
            oe_b = tb_s[shared] / tb_s[fb_s].mean()
            if M > K:
                idx = rng.choice(count_s, size=n_s, replace=False)
                oe_a = oe_a[idx]
                oe_b = oe_b[idx]
            pr = stats.pearsonr(oe_a, oe_b).statistic
            sr = stats.spearmanr(oe_a, oe_b).statistic
            rows.append([
                a, b, res, lo, hi, "oe",
                round(float(pr), 4) if np.isfinite(pr) else "NA",
                round(float(sr), 4) if np.isfinite(sr) else "NA",
                n_s,
                is_rep,
            ])

with open("{output.tsv}", "w") as handle:
    handle.write("sample1\tsample2\tresolution\tstratum_lo\tstratum_hi\tnorm\tpearson\tspearman\tn_pixels\tis_replicate\n")
    for row in sorted(rows, key=lambda r: (r[1], r[0], r[2], r[3])):
        handle.write("\t".join(map(str, row)) + "\n")

expected_header = ["sample1", "sample2", "resolution", "stratum_lo", "stratum_hi", "norm", "pearson", "spearman", "n_pixels", "is_replicate"]
with open("{output.tsv}") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
if header != expected_header:
    schema_fail("header mismatch: " + repr(header))
if len(lines) <= 1:
    schema_fail("expected > 1 data rows, got " + str(len(lines)))
for ln in lines:
    p = ln.split("\t")
    if len(p) != len(header):
        schema_fail("expected " + str(len(header)) + " columns, got " + str(len(p)) + ": " + ln)
    try:
        int(p[2])
        int(p[3])
        int(p[4])
        int(p[8])
        for k in (6, 7):
            if p[k] != "NA":
                float(p[k])
    except ValueError:
        schema_fail("non-numeric value in row: " + ln)
    if p[5] != "oe":
        schema_fail("norm must be 'oe', got " + p[5])
    if p[9] not in ("true", "false"):
        schema_fail("is_replicate must be true/false, got " + p[9])
print("computed distance-stratified concordance for " + str(len(rows)) + " sample-pair/resolution/stratum combinations")
PY
        """


rule qc_matrix_snapshot:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.balanced.mcool",
        balanced=f"{OUTDIR}/matrices/{{sample}}.mcool.balanced.done"
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
        tolerance_bp=BOUNDARY_TOLERANCE_BP,
        scripts_dir=os.path.abspath(os.path.join(workflow.basedir, "scripts"))
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/qc/boundary_overlap.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python - <<'PY' > {log} 2>&1
import itertools
import os
import sys

OUTDIR = {params.outdir!r}
SAMPLES = {params.samples!r}
tolerance_bp = {params.tolerance_bp}
sys.path.insert(0, {params.scripts_dir!r})
from boundary_matching import boundary_metrics

def load_bounds(s):
    path = OUTDIR + "/features/" + s + "/boundaries.bed"
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return {{}}
    out = {{}}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        chrom = p[0]
        mid = (int(p[1]) + int(p[2])) // 2
        out.setdefault(chrom, []).append(mid)
    for c in out:
        out[c] = sorted(out[c])
    return out

def schema_fail(msg):
    sys.stderr.write("schema check failed for boundary_overlap.tsv: " + msg + "\n")
    raise SystemExit(1)

data = {{s: load_bounds(s) for s in SAMPLES}}
with open("{output.tsv}", "w") as h:
    h.write("sample1\tsample2\tn1\tn2\tn_matched\tprecision\trecall\tf1\tjaccard\ttolerance_bp\n")
    for a, b in itertools.combinations(SAMPLES, 2):
        n1 = sum(len(v) for v in data[a].values())
        n2 = sum(len(v) for v in data[b].values())
        # one-to-one greedy nearest-neighbor matching per chromosome via the
        # shared workflow/scripts/boundary_matching.py implementation (the
        # same code the unit tests exercise)
        m = 0
        for chrom, mids_a in data[a].items():
            if chrom not in data[b]:
                continue
            m += boundary_metrics(mids_a, data[b][chrom], tolerance_bp)[0]
        precision = m / n1 if n1 else 0.0
        recall = m / n2 if n2 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
        denom = n1 + n2 - m
        jaccard = m / denom if denom else 0.0
        if not (0.0 <= jaccard <= 1.0):
            raise ValueError(
                "impossible jaccard "
                + str(jaccard)
                + " for "
                + str(a)
                + " vs "
                + str(b)
                + " (m=" + str(m) + ", n1=" + str(n1) + ", n2=" + str(n2) + ")"
            )
        h.write("\t".join(map(str, [
            a, b, n1, n2, m,
            round(precision, 6), round(recall, 6), round(f1, 6),
            round(jaccard, 6), tolerance_bp,
        ])) + "\n")

expected_header = ["sample1", "sample2", "n1", "n2", "n_matched", "precision", "recall", "f1", "jaccard", "tolerance_bp"]
with open("{output.tsv}") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
if header != expected_header:
    schema_fail("header mismatch: " + repr(header))
if len(lines) <= 1:
    schema_fail("expected > 1 data rows, got " + str(len(lines)))
for ln in lines:
    p = ln.split("\t")
    if len(p) != len(header):
        schema_fail("expected " + str(len(header)) + " columns, got " + str(len(p)) + ": " + ln)
    try:
        int(p[2])
        int(p[3])
        int(p[4])
        float(p[5])
        float(p[6])
        float(p[7])
        float(p[8])
        int(p[9])
    except ValueError:
        schema_fail("non-numeric value in row: " + ln)
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
        resolutions=COMPARTMENT_RESOLUTIONS,
        min_bins=COMPARTMENT_MIN_BINS_PER_CHROM
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/qc/compartment_correlation.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python - <<'PY' > {log} 2>&1
import itertools
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

OUTDIR = {params.outdir!r}
SAMPLES = {params.samples!r}
resolutions = {params.resolutions!r}
min_bins = {params.min_bins!r}

def chrom_type_of(chrom):
    # sex chromosomes = chrX/chrY (case-insensitive prefix)
    cl = str(chrom).lower()
    if cl.startswith("chrx") or cl.startswith("chry"):
        return "sex"
    return "autosome"

def load_e1(s, res):
    path = OUTDIR + "/features/" + s + "/compartments_" + str(res) + "bp.bedgraph"
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return None
    df = pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end", "E1"])
    df = df.dropna(subset=["E1"])
    return df.set_index(["chrom", "start"])["E1"]

def schema_fail(msg):
    sys.stderr.write("schema check failed for compartment_correlation.tsv: " + msg + "\n")
    raise SystemExit(1)

with open("{output.tsv}", "w") as h:
    h.write("sample1\tsample2\tresolution\tchrom\tpearson\tspearman\tn_bins\tchrom_type\n")
    for res in resolutions:
        tracks = {{s: load_e1(s, res) for s in SAMPLES}}
        for a, b in itertools.combinations(SAMPLES, 2):
            ta, tb = tracks[a], tracks[b]
            if ta is None or tb is None:
                continue
            both = pd.concat([ta, tb], axis=1, keys=["a", "b"]).dropna()
            for chrom, grp in both.groupby(level="chrom"):
                # drop unstable chromosomes (e.g. chrY) with too few finite bins
                if len(grp) < min_bins:
                    continue
                va = grp["a"].values
                vb = grp["b"].values
                pr = stats.pearsonr(va, vb).statistic
                sr = stats.spearmanr(va, vb).statistic
                h.write("\t".join(map(str, [
                    a, b, res, chrom,
                    round(float(pr), 4), round(float(sr), 4),
                    len(grp), chrom_type_of(chrom),
                ])) + "\n")

expected_header = ["sample1", "sample2", "resolution", "chrom", "pearson", "spearman", "n_bins", "chrom_type"]
with open("{output.tsv}") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
if header != expected_header:
    schema_fail("header mismatch: " + repr(header))
if len(lines) <= 1:
    schema_fail("expected > 1 data rows, got " + str(len(lines)))
for ln in lines:
    p = ln.split("\t")
    if len(p) != len(header):
        schema_fail("expected " + str(len(header)) + " columns, got " + str(len(p)) + ": " + ln)
    try:
        int(p[2])
        float(p[4])
        float(p[5])
        int(p[6])
    except ValueError:
        schema_fail("non-numeric value in row: " + ln)
    if p[7] not in ("autosome", "sex"):
        schema_fail("chrom_type must be autosome or sex, got " + p[7])
print("compartment correlation done")
PY
        """
