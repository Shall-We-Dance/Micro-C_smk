rule call_compartments:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    output:
        bedgraph=f"{OUTDIR}/features/{{sample}}/compartments_{{res}}bp.bedgraph"
    params:
        binsize=lambda wc: int(wc.res),
        fasta=lambda wc: config.get("reference", {}).get("bwa_indexed_fasta", "")
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
bins = clr.bins()[:][["chrom", "start", "end"]]
genome = bioframe.load_fasta("{params.fasta}")
gc_cov = bioframe.frac_gc(bins, genome)

view_df = pd.DataFrame(
    {{
        "chrom": clr.chromnames,
        "start": 0,
        "end": clr.chromsizes.values,
        "name": clr.chromnames,
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


rule call_insulation_boundaries:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    output:
        bed=f"{OUTDIR}/features/{{sample}}/boundaries.bed"
    params:
        binsize=lambda wc: int(config.get("features", {}).get("feature_resolution", config.get("matrix", {}).get("base_resolution", 1000))),
        window=lambda wc: int(config.get("features", {}).get("insulation_window", 100000)),
        balance_max_iters=lambda wc: int(config.get("features", {}).get("balance_max_iters", 500))
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/features/boundaries/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bed}) $(dirname {log})
        cooler balance --max-iters {params.balance_max_iters} --convergence-policy store_final \
          {input.mcool}::/resolutions/{params.binsize} >> {log} 2>&1 || true
        cooltools insulation --windows {params.window} --output {OUTDIR}/features/{wildcards.sample}/insulation.tsv \
          {input.mcool}::/resolutions/{params.binsize} >> {log} 2>&1 || true
        if [[ -s {OUTDIR}/features/{wildcards.sample}/insulation.tsv ]]; then
          awk 'BEGIN{{OFS="\t"}} NR>1 && $8==1 {{print $1,$2,$3}}' {OUTDIR}/features/{wildcards.sample}/insulation.tsv > {output.bed}
        else
          : > {output.bed}
          echo "[WARN] cooltools insulation failed or produced empty output; wrote empty boundaries file." >> {log}
        fi
        """


rule call_loops:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    output:
        bedpe=f"{OUTDIR}/features/{{sample}}/loops.bedpe"
    params:
        binsize=lambda wc: int(config.get("features", {}).get("feature_resolution", config.get("matrix", {}).get("base_resolution", 1000))),
        balance_max_iters=lambda wc: int(config.get("features", {}).get("balance_max_iters", 500))
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/features/loops/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bedpe}) $(dirname {log})
        cooler balance --max-iters {params.balance_max_iters} --convergence-policy store_final \
          {input.mcool}::/resolutions/{params.binsize} >> {log} 2>&1 || true
        cooltools dots --outname {output.bedpe} {input.mcool}::/resolutions/{params.binsize} >> {log} 2>&1 || true
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
