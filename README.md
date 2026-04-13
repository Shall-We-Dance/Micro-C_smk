# Micro-C Snakemake Analysis Pipeline

A modular **Snakemake** workflow for processing Micro-C paired-end FASTQ files from raw reads to contact matrices and downstream chromatin-structure analyses.

## Overview

This project is being organized to follow a distiller-nf-like logical flow while keeping each stage explicit, reproducible, and easy to extend.

### Inputs

- Paired-end FASTQ files (`R1`, `R2`) for each sample/lane.
- Reference genome and aligner index resources.
- Optional auxiliary tracks (e.g., blacklist BED).

### Outputs

- Sample-level QC reports (`fastp`, `seqkit`, `MultiQC`).
- Aligned BAM/SAM (optional retention).
- Pairwise contacts and filtered/deduplicated pairs.
- Binned contact matrices (`.cool`, multi-resolution `.mcool`, optional `.hic`).
- Normalized matrices and QC/diagnostic plots.
- Downstream feature calls and integrative analysis outputs.

### Intermediate file handling

- Raw lane FASTQs are merged per sample.
- Temporary intermediates are written under `output/tmp/` and marked as Snakemake `temp()` where appropriate.
- BAM retention is controlled by config (`output.keep_bam`).

## Pipeline blueprint (distiller-nf-style)

The target pipeline is organized into the following stages:

1. **Sample sheet / metadata validation**
   - Validate sample manifest schema (sample ID, read pairs, replicate/group labels).
   - Ensure R1/R2 lane counts are matched per sample.

2. **FASTQ QC**
   - Generate pre-processing QC summaries.
   - Capture base quality, length distribution, duplication hints, and overrepresented sequences.

3. **Adapter / quality trimming**
   - Perform adapter removal and quality trimming (currently `fastp`).
   - Optional read-level dedup where appropriate for library type.

4. **Alignment to genome (`bwa mem` / `bwa-mem2`)**
   - Align cleaned reads to the reference genome.
   - Record alignment metrics and logs for each sample.

5. **Parse BAM/SAM to contact pairs**
   - Convert alignment records into pair-format contacts.
   - Preserve required columns for downstream filtering and matrix generation.

6. **Sort pairs**
   - Coordinate- or pair-key-aware sort for deterministic downstream operations.

7. **Deduplicate pairs**
   - Remove PCR/optical duplicates at pair level.
   - Emit duplication metrics.

8. **Filter pairs**
   - Keep unique/high-quality contacts.
   - Optional blacklist filtering.
   - Optional short-distance artifact filtering.

9. **Generate stats and MultiQC**
   - Aggregate step-wise metrics into run-level summaries.
   - Produce a unified MultiQC report.

10. **Bin to contact matrices**
    - Build single-resolution `.cool` matrices.
    - Build multi-resolution `.mcool` pyramids.
    - Optionally export `.hic` for Juicebox ecosystem compatibility.

11. **Balance / normalize matrices**
    - Apply balancing/normalization (e.g., ICE/KR-style where supported).
    - Track convergence and excluded bins.

12. **QC plots**
    - cis/trans ratios.
    - distance-decay curves.
    - replicate concordance diagnostics.
    - matrix snapshots at representative loci.

13. **Downstream feature calling**
    - compartments.
    - insulation / boundaries.
    - loops / dots.
    - pileups / APA.

14. **Differential / integrative analysis**
    - Condition-level comparison.
    - Replicate-aware differential chromatin contact analysis.
    - Integration with orthogonal assays.

## Current implementation status

Implemented now:

- Sample-level lane merge for FASTQ inputs.
- `fastp` processing and report generation.
- `seqkit` stats on raw vs cleaned reads.
- `MultiQC` aggregation rule scaffold.

Planned next modules:

- Alignment (`bwa mem`/`bwa-mem2`) with configurable references.
- Pair extraction, sorting, deduplication, and filtering.
- Matrix generation (`cooler`) and normalization.
- Plotting and downstream 3D-genome feature analyses.

## Requirements

- Snakemake (recommended with conda/mamba integration).
- Conda/Mamba for rule-specific environments.
- Toolchain per module (current modules use `fastp`, `seqkit`, `multiqc`).

## Installation

Recommended: create environments on-the-fly via Snakemake.

```bash
snakemake --use-conda --cores 16
```

To speed up conda solves, consider using mamba:

```bash
snakemake --use-conda --conda-frontend mamba --cores 16
```

## Configuration

Edit `config.yaml`.

Key fields:

- `threads`: per-tool threading.
- `reference`: genome resources (to be expanded as alignment module is added).
- `samples`: sample -> list of R1/R2 FASTQs.
- `fastp`: trimming/dedup-related options.

Example:

```yaml
reference:

samples:
  sampleA:
    R1:
      - "raw/sampleA_L001_R1.fastq.gz"
      - "raw/sampleA_L002_R1.fastq.gz"
    R2:
      - "raw/sampleA_L001_R2.fastq.gz"
      - "raw/sampleA_L002_R2.fastq.gz"
```

Notes:

- The workflow assumes paired-end reads and requires both R1 and R2 lists to be the same length per sample.
- For compressed FASTQ, stream concatenation by lane is valid.
- Keep output structure stable (`output/qc`, `output/tmp`) to simplify later module integration.

## Running the workflow

```bash
snakemake -s workflow/Snakefile --use-conda --cores 16
```

Dry run:

```bash
snakemake -s workflow/Snakefile -n
```

## Contact

For questions, please open an issue or contact the pipeline maintainer.
