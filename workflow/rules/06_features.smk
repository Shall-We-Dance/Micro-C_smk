rule call_compartments:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    output:
        bedgraph=f"{OUTDIR}/features/{{sample}}/compartments.bedgraph"
    params:
        binsize=lambda wc: int(config.get("matrix", {}).get("base_resolution", 1000))
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/features/compartments/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bedgraph}) $(dirname {log})
        cooltools eigs-cis --n-eigs 3 --out-prefix {OUTDIR}/features/{wildcards.sample}/eigs \
          {input.mcool}::/resolutions/{params.binsize} > {log} 2>&1
        awk 'BEGIN{{OFS="\t"}} NR>1 {{print $1,$2,$3,$5}}' {OUTDIR}/features/{wildcards.sample}/eigs.cis.vecs.tsv > {output.bedgraph}
        """


rule call_insulation_boundaries:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    output:
        bed=f"{OUTDIR}/features/{{sample}}/boundaries.bed"
    params:
        binsize=lambda wc: int(config.get("matrix", {}).get("base_resolution", 1000)),
        window=lambda wc: int(config.get("features", {}).get("insulation_window", 100000))
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/features/boundaries/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bed}) $(dirname {log})
        cooltools insulation --windows {params.window} --output {OUTDIR}/features/{wildcards.sample}/insulation.tsv \
          {input.mcool}::/resolutions/{params.binsize} > {log} 2>&1
        awk 'BEGIN{{OFS="\t"}} NR>1 && $8==1 {{print $1,$2,$3}}' {OUTDIR}/features/{wildcards.sample}/insulation.tsv > {output.bed}
        """


rule call_loops:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    output:
        bedpe=f"{OUTDIR}/features/{{sample}}/loops.bedpe"
    params:
        binsize=lambda wc: int(config.get("matrix", {}).get("base_resolution", 1000))
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/features/loops/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bedpe}) $(dirname {log})
        cooltools dots --outname {output.bedpe} {input.mcool}::/resolutions/{params.binsize} > {log} 2>&1
        """


rule pileup_apa:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool",
        loops=f"{OUTDIR}/features/{{sample}}/loops.bedpe"
    output:
        png=f"{OUTDIR}/features/{{sample}}/apa.png"
    params:
        binsize=lambda wc: int(config.get("matrix", {}).get("base_resolution", 1000))
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/features/apa/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.png}) $(dirname {log})
        python - <<'PY'
import matplotlib.pyplot as plt
plt.figure(figsize=(4,4))
plt.text(0.5, 0.5, 'APA placeholder\nreplace with coolpuppy/cooltools pileup', ha='center', va='center')
plt.axis('off')
plt.savefig('{output.png}', dpi=200)
PY
        echo 'APA placeholder generated' > {log}
        """
