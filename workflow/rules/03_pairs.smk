rule parse_bam_to_pairs:
    input:
        bam=f"{OUTDIR}/bam/{{sample}}.name_sorted.bam"
    output:
        pairsam=temp(f"{OUTDIR}/pairs/raw/{{sample}}.pairsam.gz")
    params:
        chromsizes=lambda wc: REF["chrom_sizes"]
    threads: int(THREADS.get("pairtools", 8))
    log:
        f"logs/pairtools/parse/{{sample}}.log"
    conda:
        "envs/pairtools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.pairsam}) $(dirname {log})
        pairtools parse --chroms-path {params.chromsizes} --drop-sam --nproc-in {threads} --nproc-out {threads} \
          {input.bam} -o {output.pairsam} > {log} 2>&1
        """


rule sort_pairs:
    input:
        f"{OUTDIR}/pairs/raw/{{sample}}.pairsam.gz"
    output:
        pairsam=temp(f"{OUTDIR}/pairs/sorted/{{sample}}.sorted.pairsam.gz")
    threads: int(THREADS.get("pairtools", 8))
    log:
        f"logs/pairtools/sort/{{sample}}.log"
    conda:
        "envs/pairtools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.pairsam}) $(dirname {log})
        pairtools sort --nproc {threads} {input} -o {output.pairsam} > {log} 2>&1
        """


rule dedup_pairs:
    input:
        f"{OUTDIR}/pairs/sorted/{{sample}}.sorted.pairsam.gz"
    output:
        dedup=temp(f"{OUTDIR}/pairs/dedup/{{sample}}.dedup.pairsam.gz"),
        stats=f"{OUTDIR}/stats/pairtools/{{sample}}.dedup.stats.txt"
    threads: int(THREADS.get("pairtools", 8))
    log:
        f"logs/pairtools/dedup/{{sample}}.log"
    conda:
        "envs/pairtools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.dedup}) $(dirname {output.stats}) $(dirname {log})
        pairtools dedup --mark-dups --output-stats {output.stats} --nproc-in {threads} --nproc-out {threads} \
          {input} -o {output.dedup} > {log} 2>&1
        """


rule filter_pairs:
    input:
        dedup=f"{OUTDIR}/pairs/dedup/{{sample}}.dedup.pairsam.gz"
    output:
        pairs=f"{OUTDIR}/pairs/filtered/{{sample}}.filtered.pairs.gz"
    log:
        f"logs/pairtools/filter/{{sample}}.log"
    conda:
        "envs/pairtools.yaml"
    params:
        mapq=MIN_MAPQ,
        min_dist=(MAX_CIS_DISTANCE_ARTIFACT if MAX_CIS_DISTANCE_ARTIFACT is not None else 0),
        blacklist=(BLACKLIST_BED if BLACKLIST_ENABLED else "")
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.pairs}) $(dirname {log})

        EXPR='(pair_type=="UU") and (mapq1>={params.mapq}) and (mapq2>={params.mapq})'

        if [ "{params.min_dist}" -gt 0 ]; then
          EXPR="$EXPR and ((chrom1!=chrom2) or (abs(pos1-pos2)>={params.min_dist}))"
        fi

        pairtools select "$EXPR" {input.dedup} \
          | pairtools split --output-pairs {output.pairs} > {log} 2>&1

        if [ -n "{params.blacklist}" ]; then
          pairtools restrict -f {params.blacklist} {output.pairs} -o {output.pairs}.tmp >> {log} 2>&1
          mv {output.pairs}.tmp {output.pairs}
        fi
        """
