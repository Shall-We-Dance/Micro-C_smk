rule validate_sample_sheet:
    output:
        f"{OUTDIR}/metadata/sample_sheet.validated.tsv"
    run:
        import os
        import csv

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)

        with open(output[0], "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "lane", "r1", "r2"])

            for sample in SAMPLES:
                r1 = config["samples"][sample].get("R1", [])
                r2 = config["samples"][sample].get("R2", [])

                if not r1 or not r2:
                    raise ValueError(f"Sample {sample} has empty R1 or R2 entries in config.yaml")
                if len(r1) != len(r2):
                    raise ValueError(f"Sample {sample} has unequal R1 ({len(r1)}) and R2 ({len(r2)}) files")

                for lane_idx, (r1_path, r2_path) in enumerate(zip(r1, r2), start=1):
                    writer.writerow([sample, lane_idx, r1_path, r2_path])
