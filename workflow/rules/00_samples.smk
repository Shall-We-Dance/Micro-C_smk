rule validate_sample_sheet:
    output:
        f"{OUTDIR}/metadata/sample_sheet.validated.tsv"
    run:
        import os
        import csv

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)

        # Reference prerequisites: bwa index + .fai + chrom sizes are required
        # for every sample, so validate them once here.
        ref_fasta = config["reference"]["bwa_indexed_fasta"]
        chrom_sizes = config["reference"]["chrom_sizes"]
        for path in (ref_fasta, ref_fasta + ".fai", chrom_sizes):
            if not os.path.isfile(path):
                raise ValueError(
                    f"validate_sample_sheet: missing reference file {path}; "
                    f"check config reference.bwa_indexed_fasta / reference.chrom_sizes"
                )

        def _is_gzip(path):
            with open(path, "rb") as fh:
                return fh.read(2) == b"\x1f\x8b"

        def _first_read_name(path):
            # read the first read name from a (possibly multi-member) gz FASTQ
            import gzip
            with gzip.open(path, "rt") as fh:
                for line in fh:
                    if line.startswith("@"):
                        return line.strip()
            return ""

        def _sync_ok(r1_path, r2_path):
            n1, n2 = _first_read_name(r1_path), _first_read_name(r2_path)
            if not n1 or not n2:
                # an empty (gzip) FASTQ is invalid: validation must fail
                return False

            def _cluster_id(name):
                # Full cluster identifier of a mate read name. The ONLY
                # legitimate R1/R2 difference is the read number, which can
                # appear as the last colon field
                # (@INSTRUMENT:FLOWCELL:LANE:TILE:X:Y:1/2) or inside the
                # trailing CASAVA suffix (" 1:N:0:index" / " 2:N:0:index")
                # and the SRA " length=N" annotation. Strip those
                # read-number/length fields and compare everything else.
                tokens = name.split()
                if len(tokens) > 1 and tokens[-1].startswith(("1:", "2:")):
                    tokens.pop()  # CASAVA " <read>:Y/N:0:<index>" suffix
                if len(tokens) > 1 and tokens[-1].startswith("length="):
                    tokens.pop()  # SRA " length=<bp>" annotation
                fields = ":".join(tokens).split(":")
                if fields[-1] in ("1", "2"):
                    fields = fields[:-1]
                return ":".join(fields)

            return _cluster_id(n1) == _cluster_id(n2)

        def _check_fastq_structure(sample, lane_idx, path, max_reads, label):
            # Spot-check FASTQ record structure: each read must be exactly 4
            # lines (@header / sequence / + / quality) and the sequence and
            # quality lengths must match. Truncated or malformed records
            # (including stray blank lines) fail validation.
            import gzip
            n_reads = 0
            with gzip.open(path, "rt") as fh:
                while n_reads < max_reads:
                    header = fh.readline()
                    if not header:
                        break
                    seq, plus, qual = fh.readline(), fh.readline(), fh.readline()
                    if not seq or not plus or not qual:
                        raise ValueError(
                            f"Sample {sample} lane {lane_idx}: {label} FASTQ truncated "
                            f"at read {n_reads + 1}: {path}"
                        )
                    if not header.startswith("@") or not plus.startswith("+"):
                        raise ValueError(
                            f"Sample {sample} lane {lane_idx}: {label} FASTQ read "
                            f"{n_reads + 1} is not a 4-line record "
                            f"(header={header[:30]!r}, plus={plus[:30]!r}): {path}"
                        )
                    seq = seq.rstrip("\r\n")
                    qual = qual.rstrip("\r\n")
                    if len(seq) != len(qual):
                        raise ValueError(
                            f"Sample {sample} lane {lane_idx}: {label} FASTQ read "
                            f"{n_reads + 1} sequence length {len(seq)} != quality "
                            f"length {len(qual)}: {path}"
                        )
                    n_reads += 1
            return n_reads

        with open(output[0], "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow([
                "sample", "lane", "r1", "r2", "protocol", "enzyme",
                "base_resolution", "max_cis_artifact_dist",
                "same_fragment_max_dist", "min_mapq", "require_unique",
                "biological_sample", "condition", "library", "batch",
            ])

            for sample in SAMPLES:
                r1 = config["samples"][sample].get("R1", [])
                r2 = config["samples"][sample].get("R2", [])

                if not r1 or not r2:
                    raise ValueError(f"Sample {sample} has empty R1 or R2 entries in config.yaml")
                if len(r1) != len(r2):
                    raise ValueError(f"Sample {sample} has unequal R1 ({len(r1)}) and R2 ({len(r2)}) files")

                cfg = _sample_cfg(sample)
                enzyme = str(cfg.get("enzyme", "NA"))
                require_unique = "true" if cfg.get("require_unique", REQUIRE_UNIQUE) else "false"
                meta = {k: str(cfg.get(k, "NA")) for k in
                        ("biological_sample", "condition", "library", "batch")}

                for lane_idx, (r1_path, r2_path) in enumerate(zip(r1, r2), start=1):
                    for path in (r1_path, r2_path):
                        if not os.path.isfile(path):
                            raise ValueError(
                                f"Sample {sample} lane {lane_idx}: R1/R2 file does not exist: {path}"
                            )
                        if not _is_gzip(path):
                            raise ValueError(
                                f"Sample {sample} lane {lane_idx}: not a gzip file (bad magic bytes): {path}"
                            )
                    if not _sync_ok(r1_path, r2_path):
                        raise ValueError(
                            f"Sample {sample} lane {lane_idx}: R1 and R2 read names do not match "
                            f"(mate synchronization check failed): {r1_path} vs {r2_path}"
                        )
                    # FASTQ structure spot check: first 1000 R1 / 500 R2 reads
                    _check_fastq_structure(sample, lane_idx, r1_path, 1000, "R1")
                    _check_fastq_structure(sample, lane_idx, r2_path, 500, "R2")
                    writer.writerow([
                        sample, lane_idx, r1_path, r2_path,
                        sample_protocol(sample), enzyme,
                        sample_base_resolution(sample),
                        sample_max_cis_artifact_dist(sample),
                        sample_same_fragment_max_dist(sample),
                        sample_min_mapq(sample), require_unique,
                        meta["biological_sample"], meta["condition"],
                        meta["library"], meta["batch"],
                    ])
