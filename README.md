# Micro-C Snakemake Analysis Pipeline

This repository provides 

## Overview

### Inputs


### Outputs

### Intermediate file handling


## Pipeline steps

1. **Merge FASTQs (per-sample)**  
   Raw R1 files are concatenated into a single sample-level R1 file, and raw R2 files are concatenated into a single sample-level R2 file. (For gzip-compressed FASTQ files, stream concatenation with `cat` is valid.)

2. **fastp QC (per-sample)**  
   The merged FASTQs are processed by fastp once per sample. This produces sample-level HTML/JSON reports that are consistent with the reads used for alignment. The cleaned FASTQs from this step are temporary.

3. 

## Requirements



## Installation

Recommended: create environments on-the-fly via Snakemake.

Example:
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

* The workflow assumes paired-end reads and requires both R1 and R2 lists to be the same length per sample.
* STAR genome index construction is **not** performed in this workflow; build the STAR index in advance and point `reference.star_index` to that directory.
* BWA index construction is **not** performed in this workflow; provide prebuilt index files for `reference.bwa_indexed_fasta`.
  Example:
  `bwa index /path/to/genome.fa`

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
