# Micro-C Snakemake Pipeline (distiller-nf style)

A modular Snakemake implementation for end-to-end Micro-C processing, from raw FASTQ to matrices, QC, and downstream feature calling.

## Implemented workflow modules

1. **Sample sheet / metadata validation**
   - Validates paired-end lane definitions from `config.yaml`.
   - Writes `results/metadata/sample_sheet.validated.tsv`.

2. **Adapter / quality trimming**
   - `fastp` at sample-level (after lane concatenation).

3. **Alignment to genome (BWA-MEM/BWA-MEM2)**
   - Configurable aligner: `bwa mem` or `bwa-mem2 mem`.
   - Name-sorted BAM output.

4. **Parse BAM/SAM to contact pairs**
   - `pairtools parse`.

5. **Sort pairs**
   - `pairtools sort`.

6. **Deduplicate pairs**
   - `pairtools dedup` with per-sample dedup stats.

7. **Filter pairs**
   - Keep unique/high-quality contacts (`pair_type==UU`, MAPQ threshold).
   - Optional blacklist restriction.
   - Optional short-cis artifact filtering.

8. **Generate stats and MultiQC**
   - `pairtools stats` + fastp aggregation via MultiQC.

9. **Bin to contact matrices**
    - `.cool` via `cooler cload pairs`.
    - `.mcool` multi-resolution via `cooler zoomify`.
    - Optional `.hic` export hook via Juicer Tools.

10. **Balance / normalize matrices**
    - balancing enabled during `cooler zoomify --balance`.

11. **QC plots**
    - cis/trans proxy table.
    - distance-decay plot.
    - replicate concordance table scaffold.
    - matrix snapshot plot.

12. **Downstream feature calling**
    - compartments (`cooltools eigs-cis`).
    - insulation/boundaries (`cooltools insulation`).
    - loops/dots (`cooltools dots`).
    - APA/pileup placeholder hook.

13. **Differential / integrative analysis scaffold**
    - Produces a summary table with hooks to extend condition-wise differential analyses.

---

## Repository layout

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

Edit `config.yaml`:

- `reference.bwa_indexed_fasta`: indexed FASTA path.
- `reference.chrom_sizes`: chromosome sizes for pairtools/cooler.
- `samples.<sample>.R1/R2`: lane-level FASTQ lists.
- `alignment.aligner`: `bwa-mem2` (default) or `bwa-mem`.
- `pairs.filter.*`: filtering logic (MAPQ, blacklist, short-distance artifact threshold).
- `matrix.*`: base resolution, multires resolutions, optional hic export.

## Run

```bash
snakemake -s workflow/Snakefile --use-conda --cores 16
```

Dry-run:

```bash
snakemake -s workflow/Snakefile -n
```

## Notes

- This is designed to be extensible in a distiller-nf-like style while staying idiomatic to Snakemake.
- Some advanced modules (replicate concordance metrics, APA, differential statistics) are provided as explicit scaffolds/placeholders and can be upgraded with project-specific methods.
