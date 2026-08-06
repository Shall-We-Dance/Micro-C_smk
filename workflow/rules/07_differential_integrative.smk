if len(GROUPS) >= 2:

    rule merge_group_coolers:
        input:
            cool=lambda wc: expand(f"{OUTDIR}/matrices/{{sample}}.cool", sample=GROUPS[wc.group])
        output:
            mcool=f"{OUTDIR}/integrative/{{group}}.mcool"
        params:
            resolutions=PAIR_RES_CSV,
            threads=ZOOMIFY_THREADS
        threads: ZOOMIFY_THREADS
        conda:
            "envs/cooltools.yaml"
        log:
            f"logs/integrative/merge_{{group}}.log"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.mcool}) $(dirname {log})
            cooler merge {input.cool} -o {output.mcool}.base.cool >> {log} 2>&1
            cooler zoomify --nproc {params.threads} --resolutions {params.resolutions} --balance \
              {output.mcool}.base.cool -o {output.mcool} >> {log} 2>&1
            """

    rule differential_compare:
        input:
            group_mcools=expand(f"{OUTDIR}/integrative/{{group}}.mcool", group=list(GROUPS.keys())),
            compartments=expand(
                f"{OUTDIR}/features/{{sample}}/compartments_{{res}}bp.bedgraph",
                sample=SAMPLES,
                res=COMPARTMENT_RESOLUTIONS,
            )
        output:
            matrix=f"{OUTDIR}/integrative/matrix_compare.tsv",
            flip=f"{OUTDIR}/integrative/compartment_flip.tsv"
        params:
            resolutions=CONCORDANCE_RESOLUTIONS,
            min_chrom_size=MIN_CHROM_SIZE,
            compartment_resolutions=COMPARTMENT_RESOLUTIONS,
            samples=SAMPLES,
            groups=GROUPS,
            outdir=OUTDIR
        conda:
            "envs/cooltools.yaml"
        log:
            f"logs/integrative/differential_compare.log"
        threads: 8
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.matrix}) $(dirname {output.flip}) $(dirname {log})
            python - <<'PY' > {log} 2>&1
import os
import itertools
import cooler
import numpy as np
import pandas as pd
from scipy import stats

OUTDIR = {params.outdir!r}
SAMPLES = {params.samples!r}
GROUPS = {params.groups!r}
params_resolutions = {params.resolutions!r}
min_chrom_size = {params.min_chrom_size}
params_compartment_resolutions = {params.compartment_resolutions!r}

group_names = sorted(GROUPS.keys())

# ---- 1) genome-wide matrix comparison between group aggregates ----
with open("{output.matrix}", "w") as h:
    h.write("group1\tgroup2\tresolution\tpearson\tspearman\tn_bins\n")
for g1, g2 in itertools.combinations(group_names, 2):
    avail1 = {{int(r.rsplit("/", 1)[-1]) for r in cooler.fileops.list_coolers(OUTDIR + "/integrative/" + g1 + ".mcool")}}
    avail2 = {{int(r.rsplit("/", 1)[-1]) for r in cooler.fileops.list_coolers(OUTDIR + "/integrative/" + g2 + ".mcool")}}
    for res in params_resolutions:
        if res not in avail1 or res not in avail2:
            with open("{output.matrix}", "a") as h:
                h.write(g1 + "\t" + g2 + "\t" + str(res) + "\tNA\tNA\t0\n")
            continue
        clr1 = cooler.Cooler(OUTDIR + "/integrative/" + g1 + ".mcool::/resolutions/" + str(res))
        clr2 = cooler.Cooler(OUTDIR + "/integrative/" + g2 + ".mcool::/resolutions/" + str(res))
        chroms = [
            c
            for c, s in zip(clr1.chromnames, clr1.chromsizes)
            if s >= min_chrom_size
        ]
        v1, v2 = [], []
        for chrom in chroms:
            m1 = clr1.matrix(balance=True).fetch(chrom)
            m2 = clr2.matrix(balance=True).fetch(chrom)
            t1 = m1[np.triu_indices_from(m1, k=1)]
            t2 = m2[np.triu_indices_from(m2, k=1)]
            finite = np.isfinite(t1) & np.isfinite(t2)
            v1.append(t1[finite])
            v2.append(t2[finite])
        if not v1:
            with open("{output.matrix}", "a") as h:
                h.write(g1 + "\t" + g2 + "\t" + str(res) + "\tNA\tNA\t0\n")
            continue
        a = np.concatenate(v1)
        b = np.concatenate(v2)
        pearson = stats.pearsonr(a, b).statistic
        spearman = stats.spearmanr(a, b).statistic
        with open("{output.matrix}", "a") as h:
            h.write(g1 + "\t" + g2 + "\t" + str(res) + "\t" + format(pearson, ".4f") + "\t" + format(spearman, ".4f") + "\t" + str(len(a)) + "\n")

# ---- 2) per-chromosome compartment (E1) flip between groups ----
with open("{output.flip}", "w") as h:
    h.write("group1\tgroup2\tresolution\tchrom\tpearson\tflip_fraction\tn_bins\n")
sample_groups = {{s: g for g, ss in GROUPS.items() for s in ss}}
for g1, g2 in itertools.combinations(group_names, 2):
    for res in params_compartment_resolutions:
        tracks = {{}}
        for sample in SAMPLES:
            path = OUTDIR + "/features/" + sample + "/compartments_" + str(res) + "bp.bedgraph"
            if not os.path.isfile(path) or os.path.getsize(path) == 0:
                continue
            df = pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end", "E1"])
            df = df.dropna(subset=["E1"])
            tracks[sample] = df.set_index(["chrom", "start"])["E1"]
        if not tracks:
            continue
        samples_a = [s for s in SAMPLES if sample_groups[s] == g1]
        samples_b = [s for s in SAMPLES if sample_groups[s] == g2]
        if not samples_a or not samples_b:
            continue
        mean_a = pd.concat([tracks[s] for s in samples_a], axis=1).mean(axis=1)
        mean_b = pd.concat([tracks[s] for s in samples_b], axis=1).mean(axis=1)
        table = pd.DataFrame({{"A": mean_a, "B": mean_b}}).dropna()
        for chrom, grp in table.groupby(level="chrom"):
            if len(grp) < 10:
                continue
            va = grp["A"].values
            vb = grp["B"].values
            pearson = stats.pearsonr(va, vb).statistic
            flip = float(np.mean(np.sign(va) != np.sign(vb)))
            with open("{output.flip}", "a") as h:
                h.write(g1 + "\t" + g2 + "\t" + str(res) + "\t" + str(chrom) + "\t" + format(pearson, ".4f") + "\t" + format(flip, ".4f") + "\t" + str(len(grp)) + "\n")

print("differential compare done for groups: " + ", ".join(group_names))
PY
            """

    rule differential_summary:
        input:
            matrix=f"{OUTDIR}/integrative/matrix_compare.tsv",
            flip=f"{OUTDIR}/integrative/compartment_flip.tsv"
        output:
            summary=f"{OUTDIR}/integrative/differential_summary.tsv"
        run:
            import os

            os.makedirs(os.path.dirname(output.summary), exist_ok=True)
            n_matrix = sum(1 for _ in open(input.matrix)) - 1 if os.path.getsize(input.matrix) > 0 else 0
            n_flip = sum(1 for _ in open(input.flip)) - 1 if os.path.getsize(input.flip) > 0 else 0
            with open(output.summary, "w") as handle:
                handle.write("analysis_module\tstatus\tnote\n")
                handle.write(f"matrix_compare\tready\t{n_matrix} group-pair/resolution combinations\n")
                handle.write(f"compartment_flip\tready\t{n_flip} group-pair/chromosome combinations\n")
                handle.write("boundaries\tready\tPer-sample insulation boundaries available\n")
                handle.write("loops\tready\tPer-sample dots/loop calls available\n")
                handle.write("integration\tplanned\tJoin with ATAC/RNA/ChIP features in downstream notebooks\n")

else:

    rule differential_summary:
        output:
            summary=f"{OUTDIR}/integrative/differential_summary.tsv"
        run:
            import os

            os.makedirs(os.path.dirname(output.summary), exist_ok=True)
            with open(output.summary, "w") as handle:
                handle.write("analysis_module\tstatus\tnote\n")
                handle.write("compartments\tready\tPer-sample E1 tracks available\n")
                handle.write("boundaries\tready\tPer-sample insulation boundaries available\n")
                handle.write("loops\tready\tPer-sample dots/loop calls available\n")
                handle.write("differential\tskipped\tDefine >=2 groups in config['groups'] to enable group comparisons\n")
                handle.write("integration\tplanned\tJoin with ATAC/RNA/ChIP features in downstream notebooks\n")
