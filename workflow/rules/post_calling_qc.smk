from glob import glob


rule variant_stats:
    input:
        r="resources/Homo_sapiens_assembly38.fasta",
        f="resources/Homo_sapiens_assembly38.fasta.fai",
        vcf="results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz",
        i="results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz.tbi",
    output:
        "results/qc/bcftools_stats/joint_called_stats.out",
    benchmark:
        "results/performance_benchmarks/variant_stats/variant_stats.tsv"
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "bcftools stats "
        "--af-bins 0.01,0.05,0.1,1 "
        "-F {input.r} "
        "-s- {input.vcf} > {output}"


rule plot_variant_stats:
    input:
        "results/qc/bcftools_stats/joint_called_stats.out",
    output:
        "results/qc/bcftools_stats/plots/counts_by_af.indels.dat",
        "results/qc/bcftools_stats/plots/indels.0.dat",
        "results/qc/bcftools_stats/plots/substitutions.0.png",
        "results/qc/bcftools_stats/plots/counts_by_af.snps.dat",
        "results/qc/bcftools_stats/plots/depth.0.dat",
        "results/qc/bcftools_stats/plots/indels.0.png",
        "results/qc/bcftools_stats/plots/tstv_by_af.0.dat",
        "results/qc/bcftools_stats/plots/depth.0.png",
        "results/qc/bcftools_stats/plots/indels_by_sample.0.png",
        "results/qc/bcftools_stats/plots/tstv_by_qual.0.dat",
        "results/qc/bcftools_stats/plots/plot.py",
        "results/qc/bcftools_stats/plots/dp_by_sample.0.png",
        "results/qc/bcftools_stats/plots/plot-vcfstats.log",
        "results/qc/bcftools_stats/plots/tstv_by_sample.0.dat",
        "results/qc/bcftools_stats/plots/hets_by_sample.0.png",
        "results/qc/bcftools_stats/plots/singletons_by_sample.0.png",
        "results/qc/bcftools_stats/plots/hwe.0.dat",
        "results/qc/bcftools_stats/plots/tstv_by_sample.0.png",
        "results/qc/bcftools_stats/plots/snps_by_sample.0.png",
    params:
        d="results/qc/bcftools_stats/plots/",
    benchmark:
        "results/performance_benchmarks/plot_variant_stats/plot_variant_stats.tsv"
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "plot-vcfstats -P -p {params.d} {input}"


# ID and exclude samples with het/hom above....x?  Make tunable for WGS vs WES?  Use some outlier threshold?


rule create_ped:
    input:
        linker=config["sexLinker"],
    output:
        ped="results/qc/relatedness/sex_linker.ped",
    params:
        prefix="results/qc/relatedness/sex_linker",
    benchmark:
        "results/performance_benchmarks/create_ped/create_ped.tsv"
    shell:
        "python workflow/scripts/generate_ped.py {input} {params.prefix}"


if config["somalier"]:

    rule check_relatedness:
        input:
            vcf="results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz",
            r="resources/Homo_sapiens_assembly38.fasta",
            ped="results/qc/relatedness/sex_linker.ped",
        output:
            o1="results/qc/relatedness/somalier.html",
            o2="results/qc/relatedness/somalier.pairs.tsv",
            o3="results/qc/relatedness/somalier.samples.tsv",
        params:
            d="results/qc/relatedness/extracted/",
            o="results/qc/relatedness/somalier",
        benchmark:
            "results/performance_benchmarks/check_relatedness/check_relatedness.tsv"
        conda:
            "../envs/somalier.yaml"
        shell:
            "somalier extract "
            "-d {params.d} "
            "--sites $CONDA_PREFIX/share/somalier/sites.hg38.vcf.gz "
            "-f {input.r} {input.vcf} && "
            "somalier relate --ped {input.ped} -o {params.o} {params.d}/*.somalier"


else:

    rule mock_somalier_outputs:
        """"""
        output:
            o1=temp("results/qc/relatedness/somalier.html"),
            o2=temp("results/qc/relatedness/somalier.pairs.tsv"),
            o3=temp("results/qc/relatedness/somalier.samples.tsv"),
        shell:
            "touch {output}"


# rule per_base_coverage:
#     input:
#     output:
#     params:
#     benchmark:
#     conda:
#         "../envs/bedtools.yaml"
#     resources:
#     shell:
#         "bedtools genomecov ...."


rule sex_check:
    input:
        vcf="results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz",
        i="results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz.tbi",
    output:
        txt="results/qc/sex_check/ploidy.txt",
        png="results/qc/sex_check/ploidy.png",
    params:
        p="results/qc/sex_check/ploidy",
    benchmark:
        "results/performance_benchmarks/sex_check/sex_check.tsv"
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "bcftools +guess-ploidy -g hg38 {input.vcf} -v > {output.txt} && "
        "guess-ploidy.py {output.txt} {params.p}"


# rule sex_discordance:
#     input:
#     output:
#     shell:

if full:

    rule create_exclude_list:
        input:
            v="results/qc/contamination_check/summary.txt",
            b="results/qc/bcftools_stats/joint_called_stats.out",
        output:
            l="results/post_qc_exclusions/exclude_list.tsv",
            a="results/post_qc_exclusions/exclude_list_with_annotation.tsv",
        params:
            out="results/post_qc_exclusions/exclude_list",
            r=config["max_het_ratio"],
            d=config["min_avg_depth"],
            c=config["max_contam"],
        benchmark:
            "results/performance_benchmarks/create_exclude_list/create_exclude_list.tsv"
        conda:
            "../envs/python.yaml"
        shell:
            "python workflow/scripts/create_exclude_list.py {input.b} {params.out} --verify {input.v} -r {params.r} -d {params.d} -c {params.c}"


else:

    rule create_exclude_list:
        input:
            b="results/qc/bcftools_stats/joint_called_stats.out",
        output:
            l="results/post_qc_exclusions/exclude_list.tsv",
            a="results/post_qc_exclusions/exclude_list_with_annotation.tsv",
        params:
            out="results/post_qc_exclusions/exclude_list",
            r=config["max_het_ratio"],
            d=config["min_avg_depth"],
        benchmark:
            "results/performance_benchmarks/create_exclude_list/create_exclude_list.tsv"
        conda:
            "../envs/python.yaml"
        shell:
            "python workflow/scripts/create_exclude_list.py {input.b} {params.out} -r {params.r} -d {params.d}"


rule exclude_samples:
    input:
        v="results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz",
        i="results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz.tbi",
        l="results/post_qc_exclusions/exclude_list.tsv",
    output:
        v="results/post_qc_exclusions/samples_excluded.HC_variants.hardfiltered.vcf.gz",
        i="results/post_qc_exclusions/samples_excluded.HC_variants.hardfiltered.vcf.gz.tbi",
    benchmark:
        "results/performance_benchmarks/exclude_samples/exclude_samples.tsv"
    threads: config["bcftools"]["threads"]
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "bcftools view -S ^{input.l} --threads {threads} -Ou {input.v} | "
        "bcftools annotate --threads {threads} --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o {output.v} && "
        "tabix -p vcf {output.v}"


if full:

    rule multiqc:
        """Generate one multiQC report for all input fastqs.
        Should add samtools stats output and possibly others eventually,
        dedup metrics, ...

        update 13jan2022: combine previously split out pre-trimming qc
        data into the same report, in a second fastqc processing pass,
        and inform multiqc of how to handle this using a config yaml.
        """
        input:
            expand("results/fastqc/{rg}_r1_fastqc.zip", rg=sampleDict.keys()),
            expand("results/fastqc/{rg}_r2_fastqc.zip", rg=sampleDict.keys()),
            expand("results/post_trimming_fastqc/{rg}_r1_fastqc.zip", rg=sampleDict.keys()),
            expand("results/post_trimming_fastqc/{rg}_r2_fastqc.zip", rg=sampleDict.keys()),
            "results/qc/sex_check/ploidy.txt",
            "results/qc/relatedness/somalier.pairs.tsv",
            "results/qc/bcftools_stats/joint_called_stats.out",
            expand("results/paired_trimmed_reads/{rg}_fastp.json", rg=sampleDict.keys()),
            expand("results/dedup/{sample}.metrics.txt", sample=SAMPLES),
            expand("results/bqsr/{sample}.recal_table", sample=SAMPLES),
            expand("results/alignment_stats/{sample}.txt", sample=SAMPLES),
            expand("results/qc/contamination_check/{sample}.selfSM", sample=SAMPLES),
            "results/HaplotypeCaller/filtered/HC.variant_calling_detail_metrics",
            "results/HaplotypeCaller/filtered/HC.variant_calling_summary_metrics",
            mqc_config="config/multiqc.yaml",
        output:
            "results/multiqc/multiqc.html",
            "results/multiqc/multiqc_data/multiqc_fastqc_1.txt",
        benchmark:
            "results/performance_benchmarks/multiqc/benchmarks.tsv"
        params:
            outDir="results/multiqc/",
            outName="multiqc.html",
            inDirs="results/fastqc results/post_trimming_fastqc results/qc/sex_check results/qc/bcftools_stats results/qc/contamination_check results/paired_trimmed_reads results/dedup results/bqsr results/alignment_stats results/HaplotypeCaller/filtered",
            relatedness="results/qc/relatedness" if config["somalier"] else "",
        conda:
            "../envs/fastqc_multiqc.yaml"
        shell:
            "multiqc --force -o {params.outDir} -n {params.outName} --config {input.mqc_config} {params.inDirs} {params.relatedness}"


if jointgeno:

    rule multiqc:
        """Generate one multiQC report for joint genotyping run mode.
        Should add samtools stats output and possibly others eventually,
        dedup metrics, ...

        note that the multiqc configuration file config/multiqc.yaml
        is not meant for use with this variant of the multiqc rule.
        depending on later use cases, there may need to be two separate
        config yamls for the two different instances of the rule.
        """
        input:
            "results/qc/sex_check/ploidy.txt",
            "results/qc/relatedness/somalier.pairs.tsv",
            "results/qc/bcftools_stats/joint_called_stats.out",
            "results/HaplotypeCaller/filtered/HC.variant_calling_detail_metrics",
            "results/HaplotypeCaller/filtered/HC.variant_calling_summary_metrics",
            mqc_config="config/multiqc.yaml",
        output:
            "results/multiqc/multiqc.html",
        benchmark:
            "results/performance_benchmarks/multiqc/benchmarks.tsv"
        params:
            outDir="results/multiqc/",
            outName="multiqc.html",
            inDirs="results/qc/sex_check results/qc/bcftools_stats results/HaplotypeCaller/filtered",
            relatedness="results/qc/relatedness" if config["somalier"] else "",
        conda:
            "../envs/fastqc_multiqc.yaml"
        shell:
            "multiqc --force -o {params.outDir} --config {input.mqc_config} -n {params.outName} {params.inDirs} {params.relatedness}"


if fastq_qc_only:

    rule multiqc:
        """Generate one multiQC report for all input fastqs."""
        input:
            expand("results/fastqc/{rg}_r1_fastqc.zip", rg=sampleDict.keys()),
            expand("results/fastqc/{rg}_r2_fastqc.zip", rg=sampleDict.keys()),
            expand("results/post_trimming_fastqc/{rg}_r1_fastqc.zip", rg=sampleDict.keys()),
            expand("results/post_trimming_fastqc/{rg}_r2_fastqc.zip", rg=sampleDict.keys()),
            mqc_config="config/multiqc.yaml",
        output:
            "results/multiqc/multiqc.html",
            "results/multiqc/multiqc_data/multiqc_fastqc_1.txt",
        benchmark:
            "results/performance_benchmarks/multiqc/benchmarks.tsv"
        params:
            outDir="results/multiqc/",
            outName="multiqc.html",
            inDirs="results/fastqc results/post_trimming_fastqc",
        conda:
            "../envs/fastqc_multiqc.yaml"
        shell:
            "multiqc --force -o {params.outDir} --config {input.mqc_config} -n {params.outName} {params.inDirs}"


## To resolve the situation where performance benchmark placeholder files are
## not available at DAG construction, and to prevent race conditions, defer
## evaluation of the input to combine_benchmarks until the final report is complete.


def aggregate_benchmark_files(wildcards):
    """
    defer evaluation of benchmark file targets until after the pipeline is run
    """
    checkpoint_output = checkpoints.run_summary.get(**wildcards).output[0]
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
    output:
        "results/performance_benchmarks/benchmarking_report.html",
    params:
        threshold=config["time_threshold"],
        clusterVersion=clusterconfig,
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/benchmarking_report.Rmd"
