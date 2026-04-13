def units(sample):
    return list(range(len(config["samples"][sample]["R1"])))


rule merge_raw_fastq_per_sample:
    input:
        r1=lambda wc: [config["samples"][wc.sample]["R1"][i] for i in units(wc.sample)],
        r2=lambda wc: [config["samples"][wc.sample]["R2"][i] for i in units(wc.sample)]
    output:
        merged_r1=temp(f"{OUTDIR}/tmp/merged_raw/{{sample}}_R1.fastq.gz"),
        merged_r2=temp(f"{OUTDIR}/tmp/merged_raw/{{sample}}_R2.fastq.gz")
    threads: 2
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.merged_r1})
        cat {input.r1} > {output.merged_r1}
        cat {input.r2} > {output.merged_r2}
        """


rule fastp_sample_level:
    input:
        r1=f"{OUTDIR}/tmp/merged_raw/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/tmp/merged_raw/{{sample}}_R2.fastq.gz"
    output:
        clean_r1=temp(f"{OUTDIR}/tmp/fastp/{{sample}}_R1.fastq.gz"),
        clean_r2=temp(f"{OUTDIR}/tmp/fastp/{{sample}}_R2.fastq.gz"),
        html=f"{OUTDIR}/qc/fastp/{{sample}}/fastp.html",
        json=f"{OUTDIR}/qc/fastp/{{sample}}/fastp.json"
    log:
        f"logs/fastp/{{sample}}.log"
    threads: int(THREADS.get("fastp", 8))
    conda:
        "envs/qc.yaml"
    params:
        dedup_arg=("--dedup" if bool(config.get("fastp", {}).get("dedup_adapter", {}).get("dedup", True)) else ""),
        detect_adapter=("--detect_adapter_for_pe" if bool(config.get("fastp", {}).get("dedup_adapter", {}).get("detect_adapter_for_pe", True)) else "")
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.clean_r1}) $(dirname {output.html}) $(dirname {log})

        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.clean_r1} -O {output.clean_r2} \
          --thread {threads} \
          {params.dedup_arg} \
          {params.detect_adapter} \
          --html {output.html} --json {output.json} \
          > {log} 2>&1
        """


rule seqkit_stats_sample_level:
    input:
        raw_r1=f"{OUTDIR}/tmp/merged_raw/{{sample}}_R1.fastq.gz",
        raw_r2=f"{OUTDIR}/tmp/merged_raw/{{sample}}_R2.fastq.gz",
        clean_r1=f"{OUTDIR}/tmp/fastp/{{sample}}_R1.fastq.gz",
        clean_r2=f"{OUTDIR}/tmp/fastp/{{sample}}_R2.fastq.gz"
    output:
        raw_stats=f"{OUTDIR}/qc/seqkit/{{sample}}/raw.seqkit.stats.tsv",
        clean_stats=f"{OUTDIR}/qc/seqkit/{{sample}}/clean.seqkit.stats.tsv"
    log:
        f"logs/seqkit/{{sample}}.log"
    threads: int(THREADS.get("seqkit", 2))
    conda:
        "envs/qc.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.raw_stats}) $(dirname {log})

        seqkit stats -a -T -j {threads} {input.raw_r1} {input.raw_r2} > {output.raw_stats} 2> {log}
        seqkit stats -a -T -j {threads} {input.clean_r1} {input.clean_r2} > {output.clean_stats} 2>> {log}
        """
