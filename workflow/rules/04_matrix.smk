rule pairs_stats:
    input:
        pairs=f"{OUTDIR}/pairs/filtered/{{sample}}.filtered.pairs.gz"
    output:
        stats=f"{OUTDIR}/stats/pairtools/{{sample}}.filtered.stats.txt"
    conda:
        "envs/pairtools.yaml"
    log:
        f"logs/pairtools/stats/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.stats}) $(dirname {log})
        pairtools stats {input.pairs} > {output.stats} 2> {log}
        """


rule pairs_to_cool:
    input:
        pairs=f"{OUTDIR}/pairs/filtered/{{sample}}.filtered.pairs.gz"
    output:
        cool=f"{OUTDIR}/matrices/{{sample}}.cool"
    params:
        binsize=lambda wc: int(config.get("matrix", {}).get("base_resolution", 1000)),
        chromsizes=lambda wc: REF["chrom_sizes"]
    threads: int(THREADS.get("cooler", 8))
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/cooler/load/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.cool}) $(dirname {log})
        cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 {params.chromsizes}:{params.binsize} {input.pairs} {output.cool} \
          > {log} 2>&1
        """


rule zoomify_mcool:
    input:
        cool=f"{OUTDIR}/matrices/{{sample}}.cool"
    output:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    params:
        resolutions=PAIR_RES_CSV
    conda:
        "envs/cooltools.yaml"
    log:
        f"logs/cooler/zoomify/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.mcool}) $(dirname {log})
        cooler zoomify --resolutions {params.resolutions} --balance {input.cool} -o {output.mcool} > {log} 2>&1
        """


rule optional_hic_export:
    input:
        mcool=f"{OUTDIR}/matrices/{{sample}}.mcool"
    output:
        hic=f"{OUTDIR}/matrices/{{sample}}.hic"
    params:
        enabled=EXPORT_HIC,
        threads=HIC_EXPORT_THREADS,
        tmpdir=HIC_EXPORT_TMPDIR
    threads: HIC_EXPORT_THREADS
    conda:
        "envs/hic_tools.yaml"
    log:
        f"logs/hic_export/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.hic}) $(dirname {log})
        if [ "{params.enabled}" = "True" ]; then
          tmp_opt=""
          if [ -n "{params.tmpdir}" ]; then
            mkdir -p "{params.tmpdir}"
            tmp_opt="--tmpdir {params.tmpdir}"
          fi
          hictk convert {input.mcool} {output.hic} --threads {params.threads} $tmp_opt > {log} 2>&1
        else
          echo "hic export disabled" > {log}
          touch {output.hic}
        fi
        """


rule multiqc:
    input:
        expand(f"{OUTDIR}/qc/fastp/{{sample}}/fastp.html", sample=SAMPLES),
        expand(f"{OUTDIR}/stats/pairtools/{{sample}}.filtered.stats.txt", sample=SAMPLES),
        expand(f"{OUTDIR}/stats/pairtools/{{sample}}.dedup.stats.txt", sample=SAMPLES)
    output:
        html=f"{OUTDIR}/qc/multiqc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    conda:
        "envs/multiqc.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.html}) $(dirname {log})
        multiqc --force -o {OUTDIR}/qc/multiqc {OUTDIR} > {log} 2>&1
        """
