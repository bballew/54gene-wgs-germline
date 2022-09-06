rule bcftools_extract_final_subject_list:
    """
    Use bcftools to extract subject set in final
    vcf file
    """
    input:
        vcf="results/post_qc_exclusions/samples_excluded.HC_variants.hardfiltered.vcf.gz",
    output:
        "results/run_summary/final_subject_list.tsv",
    conda:
        "../envs/bcftools_tabix.yaml"
    benchmark:
        "results/performance_benchmarks/bcftools_extract_final_subject_list/bcftools_extract_final_subject_list.tsv"
    shell:
        "bcftools query -l {input} > {output}"


if full:

    checkpoint run_summary:
        """
        Create a run summary report aggregating information
        about overall run that is poorly captured by multiqc for the full run mode.
        """
        input:
            start_time="results/tat_tracking/start_time.txt",
            output_subject_list="results/run_summary/final_subject_list.tsv",
            r_resources="workflow/scripts/run_summary_resources.R",
            exclude_list="results/post_qc_exclusions/exclude_list_with_annotation.tsv",
            relatedness="results/qc/relatedness/somalier.pairs.tsv",
            sex="results/qc/sex_check/ploidy.txt",
            ped_file="results/qc/relatedness/sex_linker.ped",
            fastqc="results/multiqc/multiqc_data/multiqc_fastqc_1.txt",
            bcftools_stats="results/qc/bcftools_stats/joint_called_stats.out",
        output:
            report="results/run_summary/run_summary.html",
        params:
            input_samples=SAMPLES,
            out_prefix="results/run_summary/run_summary",
            somalier=config["somalier"],
            run_mode="full",
        conda:
            "../envs/r.yaml"
        benchmark:
            "results/performance_benchmarks/run_summary/run_summary.tsv"
        script:
            "../scripts/run_summary.Rmd"


if jointgeno:

    checkpoint run_summary:
        """
        Create a run summary report aggregating information
        about overall run that is poorly captured by multiqc for the joint-genotyping run mode.
        """
        input:
            start_time="results/tat_tracking/start_time.txt",
            output_subject_list="results/run_summary/final_subject_list.tsv",
            r_resources="workflow/scripts/run_summary_resources.R",
            exclude_list="results/post_qc_exclusions/exclude_list_with_annotation.tsv",
            relatedness="results/qc/relatedness/somalier.pairs.tsv",
            sex="results/qc/sex_check/ploidy.txt",
            ped_file="results/qc/relatedness/sex_linker.ped",
            bcftools_stats="results/qc/bcftools_stats/joint_called_stats.out",
        output:
            report="results/run_summary/run_summary.html",
        params:
            input_samples=SAMPLES,
            out_prefix="results/run_summary/run_summary",
            somalier=config["somalier"],
            run_mode="jointgeno",
        conda:
            "../envs/r.yaml"
        benchmark:
            "results/performance_benchmarks/run_summary/run_summary.tsv"
        script:
            "../scripts/run_summary.Rmd"


if fastq_qc_only:

    checkpoint run_summary:
        """
        Create a run summary report aggregating information
        about overall run that is poorly captured by multiqc for the fastqc only mode.
        """
        input:
            start_time="results/tat_tracking/start_time.txt",
            r_resources="workflow/scripts/run_summary_resources.R",
            fastqc="results/multiqc/multiqc_data/multiqc_fastqc_1.txt",
        output:
            report="results/run_summary/run_summary.html",
        params:
            input_samples=SAMPLES,
            out_prefix="results/run_summary/run_summary",
            somalier=False,
            run_mode="fastq_qc_only",
        conda:
            "../envs/r.yaml"
        benchmark:
            "results/performance_benchmarks/run_summary/run_summary.tsv"
        script:
            "../scripts/run_initial_summary.Rmd"


## To resolve the situation where performance benchmark placeholder files are
## not available at DAG construction, and to prevent race conditions, defer
## evaluation of the input to combine_benchmarks until the final report is complete.


def aggregate_benchmark_files(wildcards):
    """
    defer evaluation of benchmark file targets until after the pipeline is run
    """
    checkpoint_output = checkpoints.run_summary.get().output[0]
    return [j for i in expand("results/performance_benchmarks/*/*.tsv") for j in glob(i)]


rule combine_benchmarks:
    """Create a concatenated file with all the benchmarking stats generated for
    a specified set of rules defined in a list within the config as 'benchmarks',
    and append the rule name and process (i.e. sample name) to the file as columns.
    """
    input:
        tsv=aggregate_benchmark_files,
    output:
        benchmarks="results/performance_benchmarks/combined_benchmarks.tsv",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/combine_benchmarks.R"


rule benchmarking_report:
    """Take the concatenated benchmark file and generate a standard HTML report for it
    with standard plots for each metric. This report is definitely a work in progress
    and there is plenty of room for improvement in the visualizations.

    The concatenated benchmark file is passed as an input using the snakemake object in R.
    """
    input:
        benchmarks="results/performance_benchmarks/combined_benchmarks.tsv",
        start_time="results/tat_tracking/start_time.txt",
    output:
        "results/performance_benchmarks/benchmarking_report.html",
    params:
        threshold=config["time_threshold"],
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/benchmarking_report.Rmd"
