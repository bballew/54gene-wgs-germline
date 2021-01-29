rule split_multiallelics:
    """Normalize, left-align, and split mulit-allelics in the HC vcf.
    """
    input:
        vcf="results/HaplotypeCaller/genotyped/HC_variants.vcf.gz",
        index="results/HaplotypeCaller/genotyped/HC_variants.vcf.gz.tbi",
        ref="resources/Homo_sapiens_assembly38.fasta",
    output:
        vcf=temp("results/HaplotypeCaller/genotyped/HC_variants_split_multiallelics.vcf.gz"),
        index=temp("results/HaplotypeCaller/genotyped/HC_variants_split_multiallelics.vcf.gz.tbi"),
    threads: nt
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "bcftools norm "
        "-f {input.ref} "
        "-m - "
        "--threads {threads} "
        "{input.vcf} "
        "-Oz -o {output.vcf} && "
        "tabix -p vcf {output.vcf}"


rule split_snps:
    """
    """
    input:
        vcf="results/HaplotypeCaller/genotyped/HC_variants_split_multiallelics.vcf.gz",
        index="results/HaplotypeCaller/genotyped/HC_variants_split_multiallelics.vcf.gz.tbi",
    output:
        temp("results/HaplotypeCaller/filtered/snps.all.vcf.gz"),
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk SelectVariants "
        "-V {input.vcf} "
        "-select-type SNP "
        "-O {output}"


rule split_indels:
    """
    """
    input:
        vcf="results/HaplotypeCaller/genotyped/HC_variants_split_multiallelics.vcf.gz",
        index="results/HaplotypeCaller/genotyped/HC_variants_split_multiallelics.vcf.gz.tbi",
    output:
        temp("results/HaplotypeCaller/filtered/indels.all.vcf.gz"),
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk SelectVariants "
        "-V {input.vcf} "
        "-select-type INDEL "
        "-O {output}"


rule hard_filter_snps:
    """
    """
    input:
        "results/HaplotypeCaller/filtered/snps.all.vcf.gz",
    output:
        temp("results/HaplotypeCaller/filtered/snps.hardfiltered.vcf.gz"),
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk VariantFiltration "
        "-V {input} "
        '-filter "QD < 2.0" --filter-name "QD2" '
        '-filter "QUAL < 30.0" --filter-name "QUAL30" '
        '-filter "SOR > 3.0" --filter-name "SOR3" '
        '-filter "FS > 60.0" --filter-name "FS60" '
        '-filter "MQ < 40.0" --filter-name "MQ40" '
        '-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" '
        '-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" '
        "-O {output}"


rule hard_filter_indels:
    """
    """
    input:
        "results/HaplotypeCaller/filtered/indels.all.vcf.gz",
    output:
        temp("results/HaplotypeCaller/filtered/indels.hardfiltered.vcf.gz"),
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk VariantFiltration "
        "-V {input} "
        '-filter "QD < 2.0" --filter-name "QD2" '
        '-filter "QUAL < 30.0" --filter-name "QUAL30" '
        '-filter "FS > 200.0" --filter-name "FS200" '
        '-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" '
        "-O {output}"


# rule recalibrate_calls:
#     """VQSR
#     """
#     input:
#     output:
#     conda: "../envs/gatk.yaml"
#     shell:


rule merge_calls:
    """
    """
    input:
        snps="results/HaplotypeCaller/filtered/snps.hardfiltered.vcf.gz",
        indels="results/HaplotypeCaller/filtered/indels.hardfiltered.vcf.gz",
    output:
        vcf=protected("results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz"),
        index=protected("results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz.tbi"),
    params:
        tempDir,
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "bcftools concat -a -Ov {input} | "
        "bcftools sort -T {params} -Oz -o {output.vcf} && "
        "tabix -p vcf {output.vcf}" # 'gatk --java-options "-Xmx4G" GatherVcfs -I {input.snps} -I {input.indels} -O {output}'


rule picard_metrics:
    """
    """
    input:
        calls="results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz",
        index="results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz.tbi",
        dbsnp="resources/Homo_sapiens_assembly38.dbsnp138.vcf",
        d="resources/Homo_sapiens_assembly38.dict",
    output:
        "results/HaplotypeCaller/filtered/HC.variant_calling_detail_metrics",
        "results/HaplotypeCaller/filtered/HC.variant_calling_summary_metrics",
    params:
        "results/HaplotypeCaller/filtered/HC",
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk CollectVariantCallingMetrics "
        "-I {input.calls} "
        "--DBSNP {input.dbsnp} "
        "-SD {input.d} "
        "-O {params}"


# rule :
#     """
#     """
#     input:
#     output:
#     conda: "../envs/gatk.yaml"
#     shell:
