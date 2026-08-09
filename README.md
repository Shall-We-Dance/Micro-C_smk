# Micro-C Snakemake Pipeline

A modular Snakemake workflow for end-to-end Micro-C (and Hi-C) data processing, from raw FASTQ files to normalized contact matrices, quality metrics, and downstream feature outputs.

## Overview

- Protocol-aware: `protocol: micro-c | hi-c` per sample (see `config.yaml`).
- Parallel `pairtools parse` via read-name bucketing (~5.7x speedup measured on dm6).
- Optional artifact filters: Hi-C short-cis/self-ligation removal, Micro-C same-fragment (unligated/self-ligated mono-nucleosome) removal.
- QC: fastp, pairtools stats, distance decay, matrix snapshots, replicate concordance, cross-sample boundary overlap and A/B compartment correlation, MultiQC.

## Inputs

- `config.yaml` with sample metadata and pipeline parameters.
- Reference resources (`reference.*`), including indexed FASTA and chromosome sizes.
- Paired-end FASTQ files declared under `samples.*`.

## Pipeline Steps

1. Validate sample sheet and metadata.
2. Trim adapters and low-quality bases (`fastp`).
3. Align reads to the reference genome (`bwa mem` or `bwa-mem2 mem`).
4. Split the name-sorted BAM into read-name buckets and parse in parallel (`pairtools parse`), then concatenate.
5. Sort, deduplicate, and filter pairs (`pairtools sort/dedup/select`).
6. Build contact matrices and multiresolution outputs (`cooler`).
7. Feature-level outputs: compartments, insulation/boundaries, loops (dots), APA pileups.
8. QC summaries, cross-sample analyses, and a MultiQC report.

## Repository Structure

```text
workflow/
  Snakefile
  rules/
    00_samples.smk
    01_qc_trim.smk
    02_align.smk
    03_pairs.smk
    04_matrix.smk
    05_qc_plots.smk
    06_features.smk
    07_differential_integrative.smk
    envs/
```

## Configuration

Edit `config.yaml` before running the workflow:

- `reference.*`: genome references and chromosome size files.
- `samples.*`: per-sample settings and R1/R2 FASTQ lane lists.
  - `protocol: micro-c | hi-c` — Hi-C samples default to 5kb base resolution
    (`matrix.hic_base_resolution`) and a 1kb short-cis artifact cutoff
    (`matrix.hic_cis_artifact_dist`, self-ligation products of 4-cutter REs,
    cf. HiC-Pro/distiller); Micro-C uses the global 1kb base and no distance
    filter.
  - `enzyme`, `cut_sites`, `min_frag_size`, `max_frag_size`: restriction-enzyme
    metadata for Hi-C samples.
  - Per-sample overrides: `base_resolution`, `resolutions`,
    `max_cis_artifact_dist`, `same_fragment_max_dist`, `min_mapq`,
    `require_unique`.
  - `same_fragment_max_dist`: drops cis "+-" (convergent) pairs closer than
    the cutoff — unligated/self-ligated mono-nucleosome fragments in Micro-C
    (~147-200 bp); 0 = off.
- `alignment.aligner`: `bwa-mem2` (default) or `bwa-mem`.
- `pairs.parse_buckets`: number of read-name buckets for parallel parse.
- `matrix.*`: matrix resolutions, protocol defaults, optional `.hic` export
  (`hictk convert`).
- `features.*`: compartment/loop/boundary calling parameters.
- `qc.concordance_resolutions`: replicate concordance resolutions.

## Outputs

- `results/qc/`: fastp reports, pairtools stats, per-sample QC tables/plots,
  replicate concordance, cross-sample boundary overlap and compartment
  correlation tables, MultiQC report.
- `results/pairs/`: filtered pairs (+ pairix index) and pairtools stats.
- `results/matrices/`: `.cool` and `.mcool` contact matrices.
- `results/features/`: A/B compartments (per-resolution bedgraph), expected
  cis, insulation/boundaries, loops (BEDPE), APA pileups.
- `results/integrative/`: differential summaries when `groups` are defined.

## Quick Start

```bash
snakemake -s workflow/Snakefile --use-conda --cores 32
```

Dry run:

```bash
snakemake -s workflow/Snakefile -n
```

## Notes

- Intermediate files (merged FASTQ, alignments, raw/sorted pairsam) are
  cleaned up automatically; the deduplicated pairsam is kept so re-filtering
  with new `pairs.filter.*` parameters only re-runs `filter_pairs` and
  downstream steps.
- UR/RU (rescued) pairs are kept by default; they require read lengths longer
  than the fragment size, so very short-read libraries will have ~0 rescued
  pairs.
