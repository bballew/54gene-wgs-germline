rule run_summary:
    """
    Create a run summary report aggregating information
    about overall run that is poorly captured by multiqc.
    """
    input:
        r_resources="workflow/scripts/run_summary_resources.R",
    output:
        "results/run_summary/run_summary.html",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/run_summary.Rmd"
