rule variant_stats:
    input:
        r="resources/Homo_sapiens_assembly38.fasta",
        vcf="results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz",
        i="results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz.tbi",
    output:
        "results/bcftools_stats/joint_called_stats.out",
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
        "results/bcftools_stats/joint_called_stats.out",
    output:
        "results/bcftools_stats/plots/counts_by_af.indels.dat",
        "results/bcftools_stats/plots/indels.0.dat",
        "results/bcftools_stats/plots/substitutions.0.png",
        "results/bcftools_stats/plots/counts_by_af.snps.dat",
        "results/bcftools_stats/plots/indels.0.pdf",
        "results/bcftools_stats/plots/summary.log",
        "results/bcftools_stats/plots/depth.0.dat",
        "results/bcftools_stats/plots/indels.0.png",
        "results/bcftools_stats/plots/summary.tex",
        "results/bcftools_stats/plots/depth.0.pdf",
        "results/bcftools_stats/plots/indels_by_sample.0.pdf",
        "results/bcftools_stats/plots/tstv_by_af.0.dat",
        "results/bcftools_stats/plots/depth.0.png",
        "results/bcftools_stats/plots/indels_by_sample.0.png",
        "results/bcftools_stats/plots/tstv_by_qual.0.dat",
        "results/bcftools_stats/plots/dp_by_sample.0.pdf",
        "results/bcftools_stats/plots/plot.py",
        "results/bcftools_stats/plots/tstv_by_qual.0.pdf",
        "results/bcftools_stats/plots/dp_by_sample.0.png",
        "results/bcftools_stats/plots/plot-vcfstats.log",
        "results/bcftools_stats/plots/tstv_by_qual.0.png",
        "results/bcftools_stats/plots/hets_by_sample.0.pdf",
        "results/bcftools_stats/plots/singletons_by_sample.0.pdf",
        "results/bcftools_stats/plots/tstv_by_sample.0.dat",
        "results/bcftools_stats/plots/hets_by_sample.0.png",
        "results/bcftools_stats/plots/singletons_by_sample.0.png",
        "results/bcftools_stats/plots/tstv_by_sample.0.pdf",
        "results/bcftools_stats/plots/hwe.0.dat",
        "results/bcftools_stats/plots/snps_by_sample.0.pdf",
        "results/bcftools_stats/plots/tstv_by_sample.0.png",
        "results/bcftools_stats/plots/hwe.0.pdf",
        "results/bcftools_stats/plots/snps_by_sample.0.png",
        "results/bcftools_stats/plots/hwe.0.png",
        "results/bcftools_stats/plots/substitutions.0.pdf",
    params:
        d="results/bcftools_stats/plots/",
    benchmark:
        "results/performance_benchmarks/plot_variant_stats/plot_variant_stats.tsv"
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "plot-vcfstats -P -p {params.d} {input}"


# assemble and vcftools-plot
# ID and exclude samples with het/hom above....x?  Make tunable for WGS vs WES?  Use some outlier threshold?

# rule concordance_and_sample_swaps:
#     input:
#         vcf="",
#         i=""
#     output:
#     params:
#     benchmark:
#     conda:
#         "../envs/bcftools_tabix.yaml"
#     resources:
#     shell:
#         "bcftools gtcheck -g {input.vcf} {input.vcf}"


rule check_relatedness:
    input:
        vcf="results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz",
        r="resources/Homo_sapiens_assembly38.fasta",
    output:
        sites="results/relatedness/sites.hg38.vcf.gz",
        s="results/relatedness/somalier",
        o1="results/relatedness/somalier.html",
        o2="results/relatedness/somalier.pairs.tsv",
        o3="results/relatedness/somalier.samples.tsv",
        o4="results/relatedness/somalier.groups.tsv",
    params:
        d="results/relatedness/extracted/",
    benchmark:
        "results/performance_benchmarks/check_relatedness/check_relatedness.tsv"
    shell:
        "wget -O {output.s} https://github.com/brentp/somalier/releases/download/v0.2.12/somalier && "
        "wget -O {output.sites} https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz && "
        "./results/relatedness/somalier extract "
        "-d {params.d} "
        "--sites {output.sites} "
        "-f {input.r} {input.vcf} && "
        "./results/relatedness/somalier relate {params.d}/*.somalier"


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
        txt="results/sex_check/ploidy.txt",
        png="results/sex_check/ploidy.png",
    params:
        p="ploidy",
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


rule contamination_check:
    input:
        vcf="results/HaplotypeCaller/filtered/snps.hardfiltered.vcf.gz",
        i1="results/HaplotypeCaller/filtered/snps.hardfiltered.vcf.gz.tbi",
        bam="results/bqsr/{sample}.bam",
        i2="results/bqsr/{sample}.bam.bai",
    output:
        "results/contamination_check/{sample}.selfSM",
        "results/contamination_check/{sample}.selfRG",
        "results/contamination_check/{sample}.log",
        "results/contamination_check/{sample}.depthSM",
        "results/contamination_check/{sample}.depthRG",
    params:
        d="results/contamination_check/{sample}",
    benchmark:
        "results/performance_benchmarks/contamination_check/{sample}.tsv"
    conda:
        "../envs/verifybamid.yaml"
    shell:
        "verifyBamID "
        "--vcf {input.vcf} "
        "--bam {input.bam} "
        "--out {params.d}"
