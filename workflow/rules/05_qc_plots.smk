rule qc_cis_trans:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    output:
        tsv=f"{OUTDIR}/qc/plots/{{sample}}/cis_trans.tsv"
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/qc/cis_trans/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        cooltools expected-cis {input.mcool}::/resolutions/{config[matrix][base_resolution]} --view <(cooler dump -t chroms {input.mcool}::/resolutions/{config[matrix][base_resolution]}) > {output.tsv} 2> {log}
        """


rule qc_distance_decay:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    output:
        png=f"{OUTDIR}/qc/plots/{{sample}}/distance_decay.png"
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/qc/distance_decay/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.png}) $(dirname {log})
        python - <<'PY'
import matplotlib.pyplot as plt
import pandas as pd
import subprocess

cmd = "cooltools expected-cis {input.mcool}::/resolutions/{config[matrix][base_resolution]}"
df = pd.read_csv(subprocess.check_output(cmd, shell=True), sep='\t')
plt.figure(figsize=(5,4))
plt.loglog(df['dist'], df['balanced.avg'], lw=1)
plt.xlabel('distance')
plt.ylabel('contact frequency')
plt.tight_layout()
plt.savefig('{output.png}', dpi=200)
PY
        """


rule qc_replicate_concordance:
    input:
        mcools=expand(f"{OUTDIR}/matrices/{{sample}}.mcool", sample=SAMPLES)
    output:
        tsv=f"{OUTDIR}/qc/plots/replicate_concordance.tsv"
    run:
        import os
        import itertools

        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)
        with open(output.tsv, "w") as handle:
            handle.write("sample1\tsample2\tstatus\n")
            for a, b in itertools.combinations(SAMPLES, 2):
                handle.write(f"{a}\t{b}\tTODO_compute_with_genomewide_correlation\n")


rule qc_matrix_snapshot:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    output:
        png=f"{OUTDIR}/qc/plots/{{sample}}/matrix_snapshot.png"
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
import matplotlib.pyplot as plt

uri = '{input.mcool}::/resolutions/{config[matrix][base_resolution]}'
clr = cooler.Cooler(uri)
chrom = clr.chromnames[0]
mat = clr.matrix(balance=True).fetch(chrom)
plt.figure(figsize=(5,5))
plt.imshow(mat[:300,:300], cmap='fall', vmax=0.02)
plt.title('{wildcards.sample} ' + chrom)
plt.tight_layout()
plt.savefig('{output.png}', dpi=200)
PY
        """
