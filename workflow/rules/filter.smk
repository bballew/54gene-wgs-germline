rule split_multiallelics:
    """Normalize, left-align, and split mulit-allelics in the HC vcf.
    Splitting multi-allelics means we won't have each row in the vcf
    representing a unique genomic location.  But, this is allowed in
    the vcf spec, and most downstream tools are ok with it.  Plus,
    applying gatk hard filtering occurs via first separating out indels
    and snps.  Multi-type loci don't get split into either category,
    so by splitting the multi-allelics, you can apply the appropriate
    filter to all alt alleles.
    """
    input:
        vcf="results/HaplotypeCaller/genotyped/HC_variants.vcf.gz",
        index="results/HaplotypeCaller/genotyped/HC_variants.vcf.gz.tbi",
        ref="resources/Homo_sapiens_assembly38.fasta",
    output:
        vcf=temp("results/HaplotypeCaller/genotyped/HC_variants_split_multiallelics.vcf.gz"),
        index=temp("results/HaplotypeCaller/genotyped/HC_variants_split_multiallelics.vcf.gz.tbi"),
    benchmark:
        "results/performance_benchmarks/split_multiallelics/benchmarks.tsv"
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
        vcf=temp("results/HaplotypeCaller/filtered/snps.all.vcf.gz"),
        index=temp("results/HaplotypeCaller/filtered/snps.all.vcf.gz.tbi"),
    benchmark:
        "results/performance_benchmarks/split_snps/benchmarks.tsv"
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk SelectVariants "
        "-V {input.vcf} "
        "-select-type SNP "
        "-O {output.vcf}"


rule split_indels:
    """
    """
    input:
        vcf="results/HaplotypeCaller/genotyped/HC_variants_split_multiallelics.vcf.gz",
        index="results/HaplotypeCaller/genotyped/HC_variants_split_multiallelics.vcf.gz.tbi",
    output:
        vcf=temp("results/HaplotypeCaller/filtered/indels.all.vcf.gz"),
        index=temp("results/HaplotypeCaller/filtered/indels.all.vcf.gz.tbi"),
    benchmark:
        "results/performance_benchmarks/split_indels/benchmarks.tsv"
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk SelectVariants "
        "-V {input.vcf} "
        "-select-type INDEL "
        "-O {output.vcf}"


rule hard_filter_snps:
    """
    """
    input:
        vcf="results/HaplotypeCaller/filtered/snps.all.vcf.gz",
        index="results/HaplotypeCaller/filtered/snps.all.vcf.gz.tbi",
    output:
        vcf=temp("results/HaplotypeCaller/filtered/snps.hardfiltered.vcf.gz"),
        index=temp("results/HaplotypeCaller/filtered/snps.hardfiltered.vcf.gz.tbi"),
    benchmark:
        "results/performance_benchmarks/hard_filter_snps/benchmarks.tsv"
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk VariantFiltration "
        "-V {input.vcf} "
        '-filter "QD < 2.0" --filter-name "QD2" '
        '-filter "QUAL < 30.0" --filter-name "QUAL30" '
        '-filter "SOR > 3.0" --filter-name "SOR3" '
        '-filter "FS > 60.0" --filter-name "FS60" '
        '-filter "MQ < 40.0" --filter-name "MQ40" '
        '-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" '
        '-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" '
        "-O {output.vcf}"


rule hard_filter_indels:
    """
    """
    input:
        vcf="results/HaplotypeCaller/filtered/indels.all.vcf.gz",
        index="results/HaplotypeCaller/filtered/indels.all.vcf.gz.tbi",
    output:
        vcf=temp("results/HaplotypeCaller/filtered/indels.hardfiltered.vcf.gz"),
        index=temp("results/HaplotypeCaller/filtered/indels.hardfiltered.vcf.gz.tbi"),
    benchmark:
        "results/performance_benchmarks/hard_filter_indels/benchmarks.tsv"
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk VariantFiltration "
        "-V {input.vcf} "
        '-filter "QD < 2.0" --filter-name "QD2" '
        '-filter "QUAL < 30.0" --filter-name "QUAL30" '
        '-filter "FS > 200.0" --filter-name "FS200" '
        '-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" '
        "-O {output.vcf}"


# rule recalibrate_calls:
#     """VQSR
#     """
#     input:
#     output:
#     conda: "../envs/gatk.yaml"
#     shell:


rule merge_calls:
    """
    Opted for bcftools here, just because it's easy to multithread,
    sort, zip, etc. all in one data stream.
    """
    input:
        snps="results/HaplotypeCaller/filtered/snps.hardfiltered.vcf.gz",
        s_index="results/HaplotypeCaller/filtered/snps.hardfiltered.vcf.gz.tbi",
        indels="results/HaplotypeCaller/filtered/indels.hardfiltered.vcf.gz",
        i_index="results/HaplotypeCaller/filtered/indels.hardfiltered.vcf.gz.tbi",
    output:
        vcf=protected("results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz"),
        index=protected("results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz.tbi"),
    benchmark:
        "results/performance_benchmarks/merge_calls/benchmarks.tsv"
    params:
        tempDir,
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "bcftools concat -a -Ov {input.snps} {input.indels} | "
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
    benchmark:
        "results/performance_benchmarks/picard_metrics/benchmarks.tsv"
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
