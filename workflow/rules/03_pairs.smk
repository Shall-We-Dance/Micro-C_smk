import os

rule split_bam_to_buckets:
    input:
        bam=f"{OUTDIR}/bam/{{sample}}.name_sorted.bam"
    output:
        buckets=[f"{OUTDIR}/pairs/buckets/{{sample}}/bucket.{i}.bam" for i in range(PARSE_BUCKETS)]
    params:
        n_buckets=PARSE_BUCKETS,
        read_threads=int(THREADS.get("pairtools", 12)),
        outdir=f"{OUTDIR}/pairs/buckets/{{sample}}"
    # Declared threads = pysam reader (read_threads) + each of the
    # PARSE_BUCKETS writers using 2 BGZF compression threads each, matching
    # the real CPU usage inside the rule.
    threads: int(THREADS.get("pairtools", 12)) + 2 * PARSE_BUCKETS
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
        chromsizes=lambda wc: REF["chrom_sizes"],
        # threads-1 compressors (1 decompressor), derived from the same config
        # key as the threads directive (the `threads` variable is unavailable
        # in params lambdas on current snakemake).
        nproc_out=lambda wc: max(2, int(THREADS.get("pairtools", 12)) // 4) - 1
    # 1 decompressor + (threads-1) compressors, so declared threads equal the
    # real pairtools nproc usage.
    threads: max(2, int(THREADS.get("pairtools", 12)) // 4)
    log:
        f"logs/pairtools/parse/{{sample}}/bucket.{{i}}.log"
    conda:
        "envs/pairtools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.pairsam}) $(dirname {log})
        pairtools parse --chroms-path {params.chromsizes} --drop-sam --nproc-in 1 --nproc-out {params.nproc_out} \
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
        # IMPORTANT: compress with pbgzip (BGZF) -- pairtools' readers require
        # BGZF for .gz files; plain gzip (python gzip.open) is rejected by
        # `pairtools sort` ("no header or is empty").
        python - <<'PY' 2>> {log} | pbgzip -c -n 4 > {output.pairsam}
import gzip, sys
buckets = {input.buckets!r}
for bi, path in enumerate(buckets):
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                if bi == 0:
                    sys.stdout.write(line)
            else:
                sys.stdout.write(line)
sys.stderr.write("merged " + str(len(buckets)) + " buckets into {output.pairsam}\n")
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
        max_mismatch=DEDUP_MAX_MISMATCH,
        # Split the declared threads between the input and output pools;
        # nproc_in + nproc_out == threads.
        nproc_in=max(2, int(THREADS.get("pairtools", 8)) // 2),
        nproc_out=int(THREADS.get("pairtools", 8)) - max(2, int(THREADS.get("pairtools", 8)) // 2)
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
          --nproc-in {params.nproc_in} --nproc-out {params.nproc_out} \
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
        same_frag=lambda wc: sample_same_fragment_max_dist(wc.sample),
        scripts_dir=lambda wc: os.path.abspath("workflow/scripts")
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

        # NOTE: blacklist filtering below is INTERVAL filtering -- a pair is
        # dropped when chrom1/pos1 or chrom2/pos2 falls inside a blacklisted
        # BED interval. Intervals are merged per chromosome BEFORE filtering
        # so overlapping/nested intervals cannot shadow each other (review
        # P1-4). `pairtools restrict` only ANNOTATES restriction fragments
        # and does NOT remove pairs, so it must not be used here.
        if [ -n "{params.blacklist}" ]; then
          python - <<'PY' >> {log} 2>&1
import gzip, sys

sys.path.insert(0, {params.scripts_dir!r})
import blacklist as bl

# Load BED intervals and MERGE them per chromosome BEFORE filtering: the
# old per-interval lookup only checked the last interval with start <= pos,
# so positions covered by nested/overlapping intervals were missed (e.g.
# [0,100) and [50,60): pos 80 was reported as not blacklisted) -- review
# P1-4. bl.merge_intervals collapses overlapping/nested intervals into
# disjoint merged ones; bl.is_blacklisted then does an O(log n) lookup.
raw = {{}}
for bed in [{params.blacklist!r}]:
    with open(bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                continue
            raw.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))

merged = {{chrom: bl.merge_intervals(iv) for chrom, iv in raw.items()}}

def in_blacklist(chrom, pos):
    iv = merged.get(chrom)
    return iv is not None and bl.is_blacklisted(iv, pos)

kept = removed = 0
with gzip.open("{output.pairs}", "rt") as src, \
     gzip.open("{output.pairs}.tmp", "wt", compresslevel=1) as dst:
    for line in src:
        if line.startswith("#"):
            dst.write(line)
            continue
        fields = line.split("\t")
        if in_blacklist(fields[1], int(fields[2])) or in_blacklist(fields[3], int(fields[4])):
            removed += 1
            continue
        dst.write(line)
        kept += 1
print("blacklist interval filtering: kept %d pairs, removed %d pairs" % (kept, removed))
PY
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
