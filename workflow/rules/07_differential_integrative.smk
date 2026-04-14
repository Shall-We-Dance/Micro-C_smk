rule differential_integrative_summary:
    input:
        mcools=expand(f"{OUTDIR}/matrices/{{sample}}.mcool", sample=SAMPLES),
        loops=expand(f"{OUTDIR}/features/{{sample}}/loops.bedpe", sample=SAMPLES),
        boundaries=expand(f"{OUTDIR}/features/{{sample}}/boundaries.bed", sample=SAMPLES),
        compartments=expand(
            f"{OUTDIR}/features/{{sample}}/compartments_{{res}}bp.bedgraph",
            sample=SAMPLES,
            res=COMPARTMENT_RESOLUTIONS,
        ),
        concordance=f"{OUTDIR}/qc/plots/replicate_concordance.tsv"
    output:
        summary=f"{OUTDIR}/integrative/differential_summary.tsv"
    run:
        import os

        os.makedirs(os.path.dirname(output.summary), exist_ok=True)
        with open(output.summary, "w") as handle:
            handle.write("analysis_module\tstatus\tnote\n")
            handle.write("compartments\tready\tPer-sample E1 tracks available\n")
            handle.write("boundaries\tready\tPer-sample insulation boundaries available\n")
            handle.write("loops\tready\tPer-sample dots/loop calls available\n")
            handle.write("differential\tplanned\tHook in cooltools/cooler compare over condition groups\n")
            handle.write("integration\tplanned\tJoin with ATAC/RNA/ChIP features in downstream notebooks\n")
