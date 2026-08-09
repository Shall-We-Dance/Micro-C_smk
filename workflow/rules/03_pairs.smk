rule split_bam_to_buckets:
    input:
        bam=f"{OUTDIR}/bam/{{sample}}.name_sorted.bam"
    output:
        buckets=[f"{OUTDIR}/pairs/buckets/{{sample}}/bucket.{i}.bam" for i in range(PARSE_BUCKETS)]
    params:
        n_buckets=PARSE_BUCKETS,
        read_threads=int(THREADS.get("pairtools", 12)),
        outdir=f"{OUTDIR}/pairs/buckets/{{sample}}"
    threads: int(THREADS.get("pairtools", 12))
    log:
        f"logs/pairtools/split/{{sample}}.log"
    conda:
        "envs/pairtools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})
        # Single pass over the name-sorted bam: group alignments by read name
        # (mates are adjacent), hash the name into PARSE_BUCKETS buckets, write
        # per-bucket BGZF bams. pysam reads with `threads` (parallel BGZF
        # decompression); the parse step below then runs one pairtools parse
        # process per bucket, giving ~Nx parse throughput.
        python - <<'PY' > {log} 2>&1
import zlib, os, pysam
src = pysam.AlignmentFile("{input.bam}", "rb", threads={params.read_threads})
outdir = "{params.outdir}"
outs = [
    pysam.AlignmentFile(os.path.join(outdir, f"bucket.{{i}}.bam"), "wb", template=src, threads=2)
    for i in range({params.n_buckets})
]
cur = None
b = 0
n = 0
for aln in src:
    if aln.query_name != cur:
        cur = aln.query_name
        b = zlib.crc32(cur.encode()) % {params.n_buckets}
        n += 1
    outs[b].write(aln)
for o in outs:
    o.close()
print("split", n, "read names into", {params.n_buckets}, "buckets")
PY
        """


rule parse_bam_bucket:
    input:
        bam=f"{OUTDIR}/pairs/buckets/{{sample}}/bucket.{{i}}.bam"
    output:
        pairsam=f"{OUTDIR}/pairs/buckets/{{sample}}/bucket.{{i}}.pairsam.gz"
    params:
        chromsizes=lambda wc: REF["chrom_sizes"]
    threads: 3   # 1 decompressor + 2 compressors
    log:
        f"logs/pairtools/parse/{{sample}}/bucket.{{i}}.log"
    conda:
        "envs/pairtools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.pairsam}) $(dirname {log})
        pairtools parse --chroms-path {params.chromsizes} --drop-sam --nproc-in 1 --nproc-out 2 \
          --add-columns mapq --walks-policy mask \
          {input.bam} -o {output.pairsam} > {log} 2>&1
        """


rule parse_bam_to_pairs:
    input:
        buckets=lambda wc: [f"{OUTDIR}/pairs/buckets/{wc.sample}/bucket.{i}.pairsam.gz" for i in range(PARSE_BUCKETS)]
    output:
        pairsam=temp(f"{OUTDIR}/pairs/raw/{{sample}}.pairsam.gz")
    log:
        f"logs/pairtools/merge_parse/{{sample}}.log"
    conda:
        "envs/pairtools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.pairsam}) $(dirname {log})
        # Concatenate bucket pairsams: keep the header of the first bucket,
        # strip the rest. The file is NOT position-sorted (buckets are in
        # read-name order) -- pairtools sort below restores the sort order.
        python - <<'PY' > {log} 2>&1
import gzip
buckets = {input.buckets!r}
with gzip.open("{output.pairsam}", "wt") as out:
    for bi, path in enumerate(buckets):
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    if bi == 0:
                        out.write(line)
                else:
                    out.write(line)
print("merged", len(buckets), "buckets into {output.pairsam}")
PY
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
        # NOT temp(): kept on disk so re-filtering with new filter parameters
        # (mapq, max_cis_artifact_dist, same_fragment_max_dist) only re-runs
        # filter_pairs instead of the whole parse->sort->dedup chain.
        dedup=f"{OUTDIR}/pairs/dedup/{{sample}}.dedup.pairsam.gz",
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
        # Same split as parse: in/out pools each get half the declared threads.
        pairtools dedup --mark-dups --max-mismatch {params.max_mismatch} --output-stats {output.stats} \
          --nproc-in 6 --nproc-out 6 \
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
        mapq=lambda wc: sample_min_mapq(wc.sample),
        min_dist=lambda wc: sample_max_cis_artifact_dist(wc.sample),
        require_unique=lambda wc: "UU" if (_sample_cfg(wc.sample).get("require_unique", REQUIRE_UNIQUE)) else "",
        blacklist=(BLACKLIST_BED if BLACKLIST_ENABLED else ""),
        same_frag=lambda wc: sample_same_fragment_max_dist(wc.sample)
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

        if [ "{params.same_frag}" -gt 0 ]; then
          # Drop unligated/self-ligated same-fragment pairs: cis, convergent
          # (+-), closer than the cutoff. Under upper-triangle normalization
          # the stored orientation of a same-fragment pair is always "+-".
          # NOTE: use AND between the two distance inequalities: with
          # upper-triangle normalization pos1<=pos2, so (pos1-pos2<D) is
          # always true; OR would match every pair. AND correctly tests
          # |pos1-pos2| < D in both normalized and raw files.
          EXPR="$EXPR and not ((chrom1==chrom2) and (strand1=='+') and (strand2=='-') and (pos1-pos2<{params.same_frag}) and (pos2-pos1<{params.same_frag}))"
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
        # -f: overwrite a stale index when re-filtering produces a new pairs file
        pairix -f {input} > {log} 2>&1
        """
