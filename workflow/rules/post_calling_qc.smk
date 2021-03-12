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
        "results/qc/bcftools_stats/plots/tstv_by_qual.0.png",
        "results/qc/bcftools_stats/plots/tstv_by_sample.0.dat",
        "results/qc/bcftools_stats/plots/hets_by_sample.0.png",
        "results/qc/bcftools_stats/plots/singletons_by_sample.0.png",
        "results/qc/bcftools_stats/plots/hwe.0.dat",
        "results/qc/bcftools_stats/plots/tstv_by_sample.0.png",
        "results/qc/bcftools_stats/plots/snps_by_sample.0.png",
        "results/qc/bcftools_stats/plots/hwe.0.png",
    params:
        d="results/qc/bcftools_stats/plots/",
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
        sites="results/qc/relatedness/sites.hg38.vcf.gz",
        s="results/qc/relatedness/somalier",
        o1="results/qc/relatedness/somalier.html",
        o2="results/qc/relatedness/somalier.pairs.tsv",
        o3="results/qc/relatedness/somalier.samples.tsv",
        #o4="results/qc/relatedness/somalier.groups.tsv",
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

rule subset_for_contam_check:
    input:
        vcf="results/HaplotypeCaller/filtered/snps.hardfiltered.vcf.gz",
        i="results/HaplotypeCaller/filtered/snps.hardfiltered.vcf.gz.tbi",
    output:
        vcf=temp("results/qc/contamination_check/chr5_and_10.snps.hardfiltered.vcf.gz"),
        i=temp("results/qc/contamination_check/chr5_and_10.snps.hardfiltered.vcf.gz.tbi"),
    benchmark:
        "results/results/performance_benchmarks/subset_for_contam_check/subset.tsv"
    threads: 8
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "bcftools view --threads {threads} -r chr5,chr10 -Oz -o {output.vcf} {input.vcf} && "
        "tabix -p vcf {output.vcf}"


rule contamination_check:
    input:
        vcf="results/qc/contamination_check/chr5_and_10.snps.hardfiltered.vcf.gz",
        i1="results/qc/contamination_check/chr5_and_10.snps.hardfiltered.vcf.gz.tbi",
        bam="results/bqsr/{sample}.bam",
        i2="results/bqsr/{sample}.bam.bai",
    output:
        "results/qc/contamination_check/{sample}.selfSM",
        "results/qc/contamination_check/{sample}.selfRG",
        "results/qc/contamination_check/{sample}.log",
        "results/qc/contamination_check/{sample}.depthSM",
        "results/qc/contamination_check/{sample}.depthRG",
    params:
        d="results/qc/contamination_check/{sample}",
    benchmark:
        "results/performance_benchmarks/contamination_check/{sample}.tsv"
    conda:
        "../envs/verifybamid.yaml"
    resources:
        mem_mb=4000,
    shell:
        "verifyBamID "
        "--best "
        "--vcf {input.vcf} "
        "--bam {input.bam} "
        "--out {params.d}"
