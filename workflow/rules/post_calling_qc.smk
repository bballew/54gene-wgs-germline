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


rule check_relatedness:
    input:
        vcf="results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz",
        r="resources/Homo_sapiens_assembly38.fasta",
    output:
        sites="results/qc/relatedness/sites.hg38.vcf.gz",
        s="results/qc/relatedness/somalier",
        o1="results/qc/relatedness/somalier.html",
        o2="results/qc/relatedness/somalier.pairs.tsv",
        o3="results/qc/relatedness/somalier.samples.tsv",
    params:
        d="results/qc/relatedness/extracted/",
        o="results/qc/relatedness/somalier",
    benchmark:
        "results/performance_benchmarks/check_relatedness/check_relatedness.tsv"
    shell:
        "wget -O {output.s} https://github.com/brentp/somalier/releases/download/v0.2.12/somalier && "
        "wget -O {output.sites} https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz && "
        "chmod +x {output.s} && "
        "./results/qc/relatedness/somalier extract "
        "-d {params.d} "
        "--sites {output.sites} "
        "-f {input.r} {input.vcf} && "
        "./results/qc/relatedness/somalier relate -o {params.o} {params.d}/*.somalier"


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

rule combine_benchmarks:
    input:
        tsv=expand("results/performance_benchmarks/{rule}/*.tsv",rule=bench_rules),
    output: 
        "results/performance_benchmarks/combined_benchmarks.tsv"
    script:
        "scripts/combine_benchmarks.R"

rule benchmarking_report:
    input:
        benchmarks="results/performance_benchmarks/combined_benchmarks.tsv"
    output: 
        "results/performance_benchmarks/benchmarking_report.html"
    script:
        "scripts/create_benchmarking_report.Rmd"

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
        """
        input:
            expand("results/fastqc/{rg}_r1_fastqc.zip", rg=sampleDict.keys()),
            expand("results/fastqc/{rg}_r2_fastqc.zip", rg=sampleDict.keys()),
            "results/qc/sex_check/ploidy.txt",
            "results/qc/relatedness/somalier.pairs.tsv",
            "results/qc/bcftools_stats/joint_called_stats.out",
            expand("results/paired_trimmed_reads/{rg}_fastp.json", rg=sampleDict.keys()),
            expand("results/dedup/{sample}.metrics.txt", sample=SAMPLES),
            expand("results/bqsr/{sample}.recal_table",sample=SAMPLES),
            expand("results/alignment_stats/{sample}.txt", sample=SAMPLES),
            # expand("results/qc/contamination_check/{sample}.selfSM"), sample=SAMPLES,  # only if full mode
            "results/HaplotypeCaller/filtered/HC.variant_calling_detail_metrics",
            "results/HaplotypeCaller/filtered/HC.variant_calling_summary_metrics",
        output:
            "results/multiqc/multiqc.html",
        benchmark:
            "results/performance_benchmarks/multiqc/benchmarks.tsv"
        params:
            outDir="results/multiqc/",
            outName="multiqc.html",
            inDirs="results/fastqc results/qc results/paired_trimmed_reads results/dedup results/bqsr results/alignment_stats results/HaplotypeCaller/filtered",
        conda:
            "../envs/fastqc_multiqc.yaml"
        shell:
            "multiqc --force -o {params.outDir} -n {params.outName} {params.inDirs}"

if jointgeno:
    rule multiqc:
        """Generate one multiQC report for all input fastqs.
        Should add samtools stats output and possibly others eventually,
        dedup metrics, ...
        """
        input:
            "results/qc/sex_check/ploidy.txt",
            "results/qc/relatedness/somalier.pairs.tsv",
            "results/qc/bcftools_stats/joint_called_stats.out",
            "results/HaplotypeCaller/filtered/HC.variant_calling_detail_metrics",
            "results/HaplotypeCaller/filtered/HC.variant_calling_summary_metrics",
        output:
            "results/multiqc/multiqc.html",
        benchmark:
            "results/performance_benchmarks/multiqc/benchmarks.tsv"
        params:
            outDir="results/multiqc/",
            outName="multiqc.html",
            inDirs="results/qc results/HaplotypeCaller/filtered",
        conda:
            "../envs/fastqc_multiqc.yaml"
        shell:
            "multiqc --force -o {params.outDir} -n {params.outName} {params.inDirs}"
