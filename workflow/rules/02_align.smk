rule align_bwa_mem:
    input:
        r1=f"{OUTDIR}/tmp/fastp/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/tmp/fastp/{{sample}}_R2.fastq.gz"
    output:
        bam=temp(f"{OUTDIR}/bam/{{sample}}.name_sorted.bam")
    params:
        ref=lambda wc: REF["bwa_indexed_fasta"],
        aligner=lambda wc: config.get("alignment", {}).get("aligner", "bwa-mem2")
    threads: int(THREADS.get("bwa", 8))
    log:
        f"logs/bwa/{{sample}}.log"
    conda:
        "envs/alignment.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bam}) $(dirname {log})

        if [ "{params.aligner}" = "bwa-mem2" ]; then
          bwa-mem2 mem -SP -t {threads} {params.ref} {input.r1} {input.r2} 2> {log} \
            | samtools view -@ {threads} -bS - \
            | samtools sort -@ {threads} -n -o {output.bam} -
        else
          bwa mem -SP -t {threads} {params.ref} {input.r1} {input.r2} 2> {log} \
            | samtools view -@ {threads} -bS - \
            | samtools sort -@ {threads} -n -o {output.bam} -
        fi
        """
