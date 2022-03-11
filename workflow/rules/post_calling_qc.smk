rule variant_stats:
    input:
        r="resources/Homo_sapiens_assembly38.fasta",
        f="resources/Homo_sapiens_assembly38.fasta.fai",
        vcfList=expand(
            "results/HaplotypeCaller/filtered/{chrom}.hardfiltered.vcf.gz", chrom=chromList
        ),
        indexList=expand(
            "results/HaplotypeCaller/filtered/{chrom}.hardfiltered.vcf.gz.tbi", chrom=chromList
        ),
    output:
        "results/qc/bcftools_stats/{chrom}/joint_called_stats.out",
    benchmark:
        "results/performance_benchmarks/variant_stats/{chrom}/variant_stats.tsv"
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "bcftools stats "
        "--af-bins 0.01,0.05,0.1,1 "
        "-F {input.r} "
        "-s- {input.vcfList} > {output}"


rule plot_variant_stats:
    input:
        i=expand("results/qc/bcftools_stats/{chrom}/joint_called_stats.out", chrom=chromList),
    output:
        "results/qc/bcftools_stats/plots/{chrom}/counts_by_af.indels.dat",
        "results/qc/bcftools_stats/plots/{chrom}/indels.0.dat",
        "results/qc/bcftools_stats/plots/{chrom}/substitutions.0.png",
        "results/qc/bcftools_stats/plots/{chrom}/counts_by_af.snps.dat",
        "results/qc/bcftools_stats/plots/{chrom}/depth.0.dat",
        "results/qc/bcftools_stats/plots/{chrom}/indels.0.png",
        "results/qc/bcftools_stats/plots/{chrom}/tstv_by_af.0.dat",
        "results/qc/bcftools_stats/plots/{chrom}/depth.0.png",
        "results/qc/bcftools_stats/plots/{chrom}/indels_by_sample.0.png",
        "results/qc/bcftools_stats/plots/{chrom}/tstv_by_qual.0.dat",
        "results/qc/bcftools_stats/plots/{chrom}/plot.py",
        "results/qc/bcftools_stats/plots/{chrom}/dp_by_sample.0.png",
        "results/qc/bcftools_stats/plots/{chrom}/plot-vcfstats.log",
        "results/qc/bcftools_stats/plots/{chrom}/tstv_by_sample.0.dat",
        "results/qc/bcftools_stats/plots/{chrom}/hets_by_sample.0.png",
        "results/qc/bcftools_stats/plots/{chrom}/singletons_by_sample.0.png",
        "results/qc/bcftools_stats/plots/{chrom}/hwe.0.dat",
        "results/qc/bcftools_stats/plots/{chrom}/tstv_by_sample.0.png",
        "results/qc/bcftools_stats/plots/{chrom}/snps_by_sample.0.png",
    params:
        d="results/qc/bcftools_stats/plots/{chrom}",
    benchmark:
        "results/performance_benchmarks/plot_variant_stats/{chrom}/plot_variant_stats.tsv"
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "plot-vcfstats -P -p {params.d} {input.i}"


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
            vcf=(
                expand("results/bqsr/{sample}.bam", sample=SAMPLES)
                if full
                else expand(
                    "results/HaplotypeCaller/called/{sample}_all_regions.g.vcf.gz", sample=SAMPLES
                )
            ),
            r="resources/Homo_sapiens_assembly38.fasta",
            ped="results/qc/relatedness/sex_linker.ped",
        output:
            o1="results/qc/relatedness/somalier.html",
            o2="results/qc/relatedness/somalier.pairs.tsv",
            o3="results/qc/relatedness/somalier.samples.tsv",
        params:
            d="results/qc/relatedness/extracted",
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

if full:

    rule create_exclude_list:
        input:
            v=expand(
                "results/qc/contamination_check/{region}/summary.txt",
                region=config["verifyBamID"]["region"],
            ),
            b=expand("results/qc/bcftools_stats/{chrom}/joint_called_stats.out", chrom=chromList),
        output:
            l="results/post_qc_exclusions/{chrom}/exclude_list.tsv",
            a="results/post_qc_exclusions/{chrom}/exclude_list_with_annotation.tsv",
        params:
            out="results/post_qc_exclusions/{chrom}/exclude_list",
            r=config["max_het_ratio"],
            d=config["min_avg_depth"],
            c=config["max_contam"],
        benchmark:
            "results/performance_benchmarks/create_exclude_list/{chrom}/create_exclude_list.tsv"
        conda:
            "../envs/python.yaml"
        shell:
            "python workflow/scripts/create_exclude_list.py {input.b} {params.out} --verify {input.v} -r {params.r} -d {params.d} -c {params.c}"


else:

    rule create_exclude_list:
        input:
            b=expand("results/qc/bcftools_stats/{chrom}/joint_called_stats.out", chrom=chromList),
        output:
            l="results/post_qc_exclusions/{chrom}/exclude_list.tsv",
            a="results/post_qc_exclusions/{chrom}/exclude_list_with_annotation.tsv",
        params:
            out="results/post_qc_exclusions/{chrom}/exclude_list",
            r=config["max_het_ratio"],
            d=config["min_avg_depth"],
        benchmark:
            "results/performance_benchmarks/create_exclude_list/{chrom}/create_exclude_list.tsv"
        conda:
            "../envs/python.yaml"
        shell:
            "python workflow/scripts/create_exclude_list.py {input.b} {params.out} -r {params.r} -d {params.d}"


rule exclude_samples:
    input:
        vcfList=expand(
            "results/HaplotypeCaller/filtered/{chrom}.hardfiltered.vcf.gz", chrom=chromList
        ),
        indexList=expand(
            "results/HaplotypeCaller/filtered/{chrom}.hardfiltered.vcf.gz.tbi", chrom=chromList
        ),
        l=expand("results/post_qc_exclusions/{chrom}/exclude_list.tsv", chrom=chromList),
    output:
        v="results/post_qc_exclusions/samples_excluded.{chrom}.hardfiltered.vcf.gz",
        i="results/post_qc_exclusions/samples_excluded.{chrom}.hardfiltered.vcf.gz.tbi",
    benchmark:
        "results/performance_benchmarks/exclude_samples/{chrom}/exclude_samples.tsv"
    threads: config["bcftools"]["threads"]
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "bcftools view -S ^{input.l} --threads {threads} -Ou {input.vcfList} | "
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
            "results/qc/relatedness/somalier.pairs.tsv",
            expand("results/qc/bcftools_stats/{chrom}/joint_called_stats.out", chrom=chromList),
            expand("results/paired_trimmed_reads/{rg}_fastp.json", rg=sampleDict.keys()),
            expand("results/dedup/{sample}.metrics.txt", sample=SAMPLES),
            expand("results/bqsr/{sample}.recal_table", sample=SAMPLES),
            expand("results/alignment_stats/{sample}.txt", sample=SAMPLES),
            expand(
                "results/qc/contamination_check/{region}/{sample}.selfSM",
                sample=SAMPLES,
                region=config["verifyBamID"]["region"],
            ),
            expand(
                "results/HaplotypeCaller/filtered/{chrom}.variant_calling_detail_metrics",
                chrom=chromList,
            ),
            expand(
                "results/HaplotypeCaller/filtered/{chrom}.variant_calling_summary_metrics",
                chrom=chromList,
            ),
            mqc_config="config/multiqc.yaml",
        output:
            "results/multiqc/multiqc.html",
            "results/multiqc/multiqc_data/multiqc_fastqc_1.txt",
        benchmark:
            "results/performance_benchmarks/multiqc/benchmarks.tsv"
        params:
            outDir="results/multiqc/",
            outName="multiqc.html",
            inDirs="results/fastqc results/post_trimming_fastqc results/qc/bcftools_stats results/qc/contamination_check results/paired_trimmed_reads results/dedup results/bqsr results/alignment_stats results/HaplotypeCaller/filtered",
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
            "results/qc/relatedness/somalier.pairs.tsv",
            expand("results/qc/bcftools_stats/{chrom}/joint_called_stats.out", chrom=chromList),
            expand(
                "results/HaplotypeCaller/filtered/{chrom}.variant_calling_detail_metrics",
                chrom=chromList,
            ),
            expand(
                "results/HaplotypeCaller/filtered/{chrom}.variant_calling_summary_metrics",
                chrom=chromList,
            ),
            mqc_config="config/multiqc.yaml",
        output:
            "results/multiqc/multiqc.html",
        benchmark:
            "results/performance_benchmarks/multiqc/benchmarks.tsv"
        params:
            outDir="results/multiqc/",
            outName="multiqc.html",
            inDirs="results/qc/bcftools_stats results/HaplotypeCaller/filtered",
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
