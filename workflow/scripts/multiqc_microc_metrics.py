"""Build MultiQC custom-content tables from Micro-C workflow outputs."""

import csv
import json
import math
from pathlib import Path
from typing import Any

import yaml


OUTDIR = Path(str(snakemake.params.outdir))
SAMPLES = [str(s) for s in snakemake.params.samples]


def parse_pairtools(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open() as handle:
        for line in handle:
            if "\t" in line:
                key, value = line.rstrip("\n").split("\t", 1)
                values[key] = value
    return values


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(value: Any, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def as_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def pct(numerator: float, denominator: float) -> float | str:
    if denominator <= 0 or not math.isfinite(numerator):
        return "NA"
    return round(100.0 * numerator / denominator, 3)


def finite_or_na(value: Any, digits: int = 4) -> float | str:
    number = as_float(value)
    return round(number, digits) if math.isfinite(number) else "NA"


def count_nonempty(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open() as handle:
        return sum(1 for line in handle if line.strip())


def dump_custom(path: str, payload: dict[str, Any]) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False, allow_unicode=True)


def metric_header(title: str, description: str, suffix: str = "", number_format: str = "") -> dict[str, Any]:
    header: dict[str, Any] = {
        "title": title,
        "description": description,
        "scale": False,
    }
    if suffix:
        header["suffix"] = suffix
    if number_format:
        header["format"] = number_format
    return header


library_data: dict[str, dict[str, Any]] = {}
for sample in SAMPLES:
    with (OUTDIR / "qc/fastp" / sample / "fastp.json").open() as handle:
        fastp = json.load(handle)
    before = fastp["summary"]["before_filtering"]
    after = fastp["summary"]["after_filtering"]
    dedup = parse_pairtools(OUTDIR / "stats/pairtools" / f"{sample}.dedup.stats.txt")
    filtered = parse_pairtools(OUTDIR / "stats/pairtools" / f"{sample}.filtered.stats.txt")
    apa_rows = read_tsv(OUTDIR / "features" / sample / "apa.tsv")
    apa = apa_rows[0] if apa_rows else {}

    raw_pairs = as_int(before.get("total_reads")) // 2
    retained_pairs = as_int(after.get("total_reads")) // 2
    total_parsed = as_int(dedup.get("total"))
    mapped = as_int(dedup.get("total_mapped"))
    unmapped = as_int(dedup.get("total_unmapped"))
    one_sided = as_int(dedup.get("total_single_sided_mapped"))
    duplicate_pairs = as_int(dedup.get("total_dups"))
    unique_pairs = as_int(dedup.get("total_nodups"))
    valid_pairs = as_int(filtered.get("total"))
    cis = as_int(filtered.get("cis"))
    cis_1kb = as_int(filtered.get("cis_1kb+"))
    cis_10kb = as_int(filtered.get("cis_10kb+"))
    q30_rate = as_float(after.get("q30_rate"))

    library_data[sample] = {
        "raw_pairs": raw_pairs,
        "q30_pct": round(100.0 * q30_rate, 3),
        "retained_pct": pct(retained_pairs, raw_pairs),
        "paired_mapped_pct": pct(mapped, total_parsed),
        "unmapped_pct": pct(unmapped, total_parsed),
        "single_sided_pct": pct(one_sided, total_parsed),
        "duplicate_pct": pct(duplicate_pairs, mapped),
        "unique_pairs": unique_pairs,
        "valid_pairs": valid_pairs,
        "valid_yield_pct": pct(valid_pairs, raw_pairs),
        "cis_pct": pct(cis, valid_pairs),
        "cis_1kb_pct": pct(cis_1kb, valid_pairs),
        "cis_10kb_pct": pct(cis_10kb, valid_pairs),
        "boundaries": count_nonempty(OUTDIR / "features" / sample / "boundaries.bed"),
        "loops": as_int(apa.get("n_loops"), max(0, count_nonempty(OUTDIR / "features" / sample / "loops.bedpe") - 1)),
        "apa_ratio": finite_or_na(apa.get("ratio")),
        "apa_mean": finite_or_na(apa.get("mean", apa.get("center_mean")), 6),
        "apa_max": finite_or_na(apa.get("max"), 6),
    }

library_payload = {
    "id": "microc_library_qc",
    "section_name": "Micro-C library QC",
    "description": (
        "Per-library sequencing, mapping, complexity, valid-contact and structural-feature metrics. "
        "Pairtools duplicate percentage uses mapped pairs as the denominator; valid yield uses raw read pairs."
    ),
    "plot_type": "table",
    "pconfig": {
        "id": "microc_library_qc_table",
        "title": "Micro-C library quality summary",
        "col1_header": "Sample",
        "sort_rows": False,
    },
    "headers": {
        "raw_pairs": metric_header("Raw pairs", "Paired-end fragments before fastp"),
        "q30_pct": metric_header("Q30", "Bases at or above Q30 after fastp", "%"),
        "retained_pct": metric_header("Reads retained", "Read pairs retained by fastp", "%"),
        "paired_mapped_pct": metric_header("Paired mapped", "Pairs with both ends mapped", "%"),
        "unmapped_pct": metric_header("Unmapped", "Fully unmapped parsed pairs", "%"),
        "single_sided_pct": metric_header("Single-sided", "Pairs with only one mapped end", "%"),
        "duplicate_pct": metric_header("Duplicates", "Pairtools duplicates divided by mapped pairs", "%"),
        "unique_pairs": metric_header("Unique pairs", "Mapped non-duplicate pairs before final selection"),
        "valid_pairs": metric_header("Valid pairs", "Pairs retained after final Micro-C filters"),
        "valid_yield_pct": metric_header("Valid yield", "Final valid pairs divided by raw read pairs", "%"),
        "cis_pct": metric_header("Cis", "Cis fraction among final valid pairs", "%"),
        "cis_1kb_pct": metric_header("Cis >=1 kb", "Final valid cis pairs separated by at least 1 kb", "%"),
        "cis_10kb_pct": metric_header("Cis >=10 kb", "Final valid cis pairs separated by at least 10 kb", "%"),
        "boundaries": metric_header("Boundaries", "Number of called insulation boundaries"),
        "loops": metric_header("Loops", "Number of loop calls used for APA"),
        "apa_ratio": metric_header("APA ratio", "APA center-to-corner enrichment when available", number_format="{:.6g}"),
        "apa_mean": metric_header("APA mean", "Legacy APA mean or center mean when available", number_format="{:.6g}"),
        "apa_max": metric_header("APA max", "Legacy maximum APA signal when available", number_format="{:.6g}"),
    },
    "data": library_data,
}
dump_custom(snakemake.output.library, library_payload)


matrix_data: dict[str, dict[str, Any]] = {}
for i, row in enumerate(read_tsv(OUTDIR / "qc/plots/replicate_concordance.tsv"), start=1):
    resolution = as_int(row.get("resolution"))
    label = f"{row.get('sample1', 'sample1')} vs {row.get('sample2', 'sample2')} @ {resolution:,} bp"
    if "stratum_lo" in row:
        label += f" [{row.get('stratum_lo')}, {row.get('stratum_hi')})"
    if label in matrix_data:
        label += f" #{i}"
    matrix_data[label] = {
        "resolution": resolution,
        "stratum_lo": as_int(row.get("stratum_lo")) if "stratum_lo" in row else "NA",
        "stratum_hi": as_int(row.get("stratum_hi")) if "stratum_hi" in row else "NA",
        "pearson": finite_or_na(row.get("pearson")),
        "spearman": finite_or_na(row.get("spearman")),
        "n_observations": as_int(row.get("n_pixels", row.get("n_bins"))),
        "is_replicate": row.get("is_replicate", "NA"),
    }

matrix_payload = {
    "id": "microc_matrix_concordance",
    "section_name": "Micro-C matrix concordance",
    "description": "Pairwise contact-matrix Pearson and Spearman correlations.",
    "plot_type": "table",
    "pconfig": {"id": "microc_matrix_concordance_table", "col1_header": "Comparison"},
    "headers": {
        "resolution": metric_header("Resolution", "Matrix resolution", " bp"),
        "stratum_lo": metric_header("Distance from", "Lower distance bound when stratified", " bp"),
        "stratum_hi": metric_header("Distance to", "Upper distance bound; 0 means open-ended", " bp"),
        "pearson": metric_header("Pearson", "Pearson correlation", number_format="{:.4f}"),
        "spearman": metric_header("Spearman", "Spearman rank correlation", number_format="{:.4f}"),
        "n_observations": metric_header("N", "Compared bins or pixels"),
        "is_replicate": metric_header("Replicate", "Whether sample metadata marks the pair as biological replicates"),
    },
    "data": matrix_data,
}
dump_custom(snakemake.output.matrix, matrix_payload)


boundary_data: dict[str, dict[str, Any]] = {}
for i, row in enumerate(read_tsv(OUTDIR / "qc/plots/boundary_overlap.tsv"), start=1):
    label = f"{row.get('sample1', 'sample1')} vs {row.get('sample2', 'sample2')}"
    if label in boundary_data:
        label += f" #{i}"
    n1 = as_int(row.get("n1"))
    n2 = as_int(row.get("n2"))
    matched = as_int(row.get("n_matched", row.get("n_overlap")))
    precision = as_float(row.get("precision", row.get("jaccard_s1_over_s2")))
    recall = as_float(row.get("recall", row.get("jaccard_s2_over_s1")))
    f1 = as_float(row.get("f1"))
    if not math.isfinite(f1) and math.isfinite(precision) and math.isfinite(recall) and precision + recall:
        f1 = 2.0 * precision * recall / (precision + recall)
    boundary_data[label] = {
        "n1": n1,
        "n2": n2,
        "matched": matched,
        "precision": finite_or_na(precision),
        "recall": finite_or_na(recall),
        "f1": finite_or_na(f1),
        "jaccard": finite_or_na(row.get("jaccard", row.get("n_jaccard"))),
        "tolerance_bp": as_int(row.get("tolerance_bp")) if row.get("tolerance_bp") else "NA",
    }

boundary_payload = {
    "id": "microc_boundary_overlap",
    "section_name": "Micro-C boundary overlap",
    "description": "Pairwise insulation-boundary overlap statistics.",
    "plot_type": "table",
    "pconfig": {"id": "microc_boundary_overlap_table", "col1_header": "Comparison"},
    "headers": {
        "n1": metric_header("Boundaries 1", "Boundary count in sample 1"),
        "n2": metric_header("Boundaries 2", "Boundary count in sample 2"),
        "matched": metric_header("Matched", "Matched or overlapping boundaries"),
        "precision": metric_header("Precision", "Matched divided by boundaries in sample 1", number_format="{:.4f}"),
        "recall": metric_header("Recall", "Matched divided by boundaries in sample 2", number_format="{:.4f}"),
        "f1": metric_header("F1", "Harmonic mean of precision and recall", number_format="{:.4f}"),
        "jaccard": metric_header("Jaccard", "Boundary-set Jaccard index", number_format="{:.4f}"),
        "tolerance_bp": metric_header("Tolerance", "Matching tolerance", " bp"),
    },
    "data": boundary_data,
}
dump_custom(snakemake.output.boundary, boundary_payload)


compartment_data: dict[str, dict[str, Any]] = {}
for i, row in enumerate(read_tsv(OUTDIR / "qc/plots/compartment_correlation.tsv"), start=1):
    label = (
        f"{row.get('sample1', 'sample1')} vs {row.get('sample2', 'sample2')} "
        f"| {row.get('chrom', 'genome')}"
    )
    if label in compartment_data:
        label += f" #{i}"
    compartment_data[label] = {
        "resolution": as_int(row.get("resolution")),
        "chromosome": row.get("chrom", "NA"),
        "pearson": finite_or_na(row.get("pearson")),
        "spearman": finite_or_na(row.get("spearman")),
        "n_bins": as_int(row.get("n_bins")),
    }

compartment_payload = {
    "id": "microc_compartment_correlation",
    "section_name": "Micro-C compartment correlation",
    "description": "Chromosome-level correlation of compartment eigenvector tracks.",
    "plot_type": "table",
    "pconfig": {"id": "microc_compartment_correlation_table", "col1_header": "Comparison"},
    "headers": {
        "resolution": metric_header("Resolution", "Compartment resolution", " bp"),
        "chromosome": metric_header("Chromosome", "Chromosome or contig"),
        "pearson": metric_header("Pearson", "Pearson correlation of E1 values", number_format="{:.4f}"),
        "spearman": metric_header("Spearman", "Spearman correlation of E1 values", number_format="{:.4f}"),
        "n_bins": metric_header("Bins", "Finite bins included in the comparison"),
    },
    "data": compartment_data,
}
dump_custom(snakemake.output.compartment, compartment_payload)

print(f"Wrote MultiQC custom content for {len(SAMPLES)} Micro-C libraries")
