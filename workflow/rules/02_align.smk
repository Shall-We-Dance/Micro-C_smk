rule align_bwa_mem:
    input:
        r1=f"{OUTDIR}/tmp/fastp/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/tmp/fastp/{{sample}}_R2.fastq.gz"
    output:
        bam=temp(f"{OUTDIR}/bam/{{sample}}.name_sorted.bam")
    params:
        ref=lambda wc: REF["bwa_indexed_fasta"],
        aligner=lambda wc: config.get("alignment", {}).get("aligner", "bwa-mem2")
    # Thread budget = bwa + samtools view + samtools sort (all run in parallel
    # within the same rule). Previously {threads} was passed to all three,
    # tripling the real CPU usage vs. what snakemake accounted for (e.g. 4
    # concurrent jobs used ~192 cores on a 72-core node). bwa gets the bulk;
    # view/sort get the remainder.
    threads: 18
    log:
        f"logs/bwa/{{sample}}.log"
    conda:
        "envs/alignment.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bam}) $(dirname {log})

        # samtools sort -o writes chunk files named {output}.tmp.*; leftovers
        # from a killed run make reruns fail with "File exists". Clear them.
        rm -f {output.bam}.tmp.*

        if [ "{params.aligner}" = "bwa-mem2" ]; then
          bwa-mem2 mem -SP -t 10 {params.ref} {input.r1} {input.r2} 2> {log} \
            | samtools view -@ 4 -bS - \
            | samtools sort -@ 4 -n -o {output.bam} -
        else
          bwa mem -SP -t 10 {params.ref} {input.r1} {input.r2} 2> {log} \
            | samtools view -@ 4 -bS - \
            | samtools sort -@ 4 -n -o {output.bam} -
        fi
        """
