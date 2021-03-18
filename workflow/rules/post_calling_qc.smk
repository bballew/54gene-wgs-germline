rule variant_stats:
    input:
        r="resources/Homo_sapiens_assembly38.fasta",
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
        # "results/qc/bcftools_stats/plots/tstv_by_qual.0.png",
        "results/qc/bcftools_stats/plots/tstv_by_sample.0.dat",
        "results/qc/bcftools_stats/plots/hets_by_sample.0.png",
        "results/qc/bcftools_stats/plots/singletons_by_sample.0.png",
        "results/qc/bcftools_stats/plots/hwe.0.dat",
        "results/qc/bcftools_stats/plots/tstv_by_sample.0.png",
        "results/qc/bcftools_stats/plots/snps_by_sample.0.png",
        # "results/qc/bcftools_stats/plots/hwe.0.png",
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

if full:

    rule create_exclude_list:
        input:
            v="results/qc/contamination_check/summary.txt",
            b="results/qc/bcftools_stats/joint_called_stats.out",
        output:
            l="results/post_qc_exclusions/exclude_list.tsv",
            a="results/post_qc_exclusions/exclude_list_with_annotation.tsv",
        params:
            "results/post_qc_exclusions/exclude_list",
        benchmark:
            "results/performance_benchmarks/create_exclude_list/create_exclude_list.tsv"
        conda:
            "../envs/python.yaml"
        shell:
            "python workflow/scripts/create_exclude_list.py {input.b} {params} --verify {input.v}"


else:

    rule create_exclude_list:
        input:
            b="results/qc/bcftools_stats/joint_called_stats.out",
        output:
            l="results/post_qc_exclusions/exclude_list.tsv",
            a="results/post_qc_exclusions/exclude_list_with_annotation.tsv",
        params:
            "results/post_qc_exclusions/exclude_list",
        benchmark:
            "results/performance_benchmarks/create_exclude_list/create_exclude_list.tsv"
        conda:
            "../envs/python.yaml"
        shell:
            "python workflow/scripts/create_exclude_list.py {input.b} {params}"


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
    threads: 4
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "bcftools view -S ^{input.l} --threads {threads} -Oz -o {output.v} {input.v} && "
        "tabix -p vcf {output.v}"
