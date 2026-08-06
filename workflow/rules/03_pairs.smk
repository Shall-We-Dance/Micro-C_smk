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
          --add-columns mapq --walks-policy mask \
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
    params:
        max_mismatch=DEDUP_MAX_MISMATCH
    threads: int(THREADS.get("pairtools", 8))
    log:
        f"logs/pairtools/dedup/{{sample}}.log"
    conda:
        "envs/pairtools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.dedup}) $(dirname {output.stats}) $(dirname {log})
        pairtools dedup --mark-dups --max-mismatch {params.max_mismatch} --output-stats {output.stats} \
          --nproc-in {threads} --nproc-out {threads} \
          {input} -o {output.dedup} > {log} 2>&1
        """


rule filter_pairs:
    input:
        dedup=f"{OUTDIR}/pairs/dedup/{{sample}}.dedup.pairsam.gz"
    output:
        pairs=f"{OUTDIR}/pairs/filtered/{{sample}}.filtered.pairs.gz",
        stats=f"{OUTDIR}/stats/pairtools/{{sample}}.filtered.stats.txt"
    log:
        f"logs/pairtools/filter/{{sample}}.log"
    conda:
        "envs/pairtools.yaml"
    params:
        mapq=MIN_MAPQ,
        min_dist=(MAX_CIS_DISTANCE_ARTIFACT if MAX_CIS_DISTANCE_ARTIFACT else 0),
        require_unique="UU" if REQUIRE_UNIQUE else "",
        blacklist=(BLACKLIST_BED if BLACKLIST_ENABLED else "")
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.pairs}) $(dirname {output.stats}) $(dirname {log})

        # distiller-parity default: keep pair types UU, UR, RU (rescued
        # single-ligation products); drop M (multi-mapped/low-mapq) and
        # W (masked walkers) pairs. If require_unique is true, keep UU only.
        if [ -n "{params.require_unique}" ]; then
          EXPR='(pair_type=="UU") and (mapq1>={params.mapq}) and (mapq2>={params.mapq})'
        else
          EXPR='((pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")) and (mapq1>={params.mapq}) and (mapq2>={params.mapq})'
        fi

        if [ "{params.min_dist}" -gt 0 ]; then
          # Absolute-distance check without Python builtins (abs() may be
          # unavailable in pairtools select's restricted eval context).
          EXPR="$EXPR and ((chrom1!=chrom2) or ((pos1-pos2>={params.min_dist}) or (pos2-pos1>={params.min_dist})))"
        fi

        # parse_bam_to_pairs uses --drop-sam, so the stream has pairs columns
        # only. Write filtered pairs directly instead of `pairtools split`.
        pairtools select "$EXPR" {input.dedup} -o {output.pairs} > {log} 2>&1

        if [ -n "{params.blacklist}" ]; then
          pairtools restrict -f {params.blacklist} {output.pairs} -o {output.pairs}.tmp >> {log} 2>&1
          mv {output.pairs}.tmp {output.pairs}
        fi

        pairtools stats {output.pairs} > {output.stats} 2>> {log}
        """


rule index_pairs:
    input:
        f"{OUTDIR}/pairs/filtered/{{sample}}.filtered.pairs.gz"
    output:
        px2=f"{OUTDIR}/pairs/filtered/{{sample}}.filtered.pairs.gz.px2"
    conda:
        "envs/pairtools.yaml"
    log:
        f"logs/pairtools/index/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.px2}) $(dirname {log})
        pairix {input} > {log} 2>&1
        """
