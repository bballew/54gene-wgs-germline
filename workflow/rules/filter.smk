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
        vcf="results/HaplotypeCaller/genotyped/chrom_vcfs/{chrom}.vcf.gz",
        index="results/HaplotypeCaller/genotyped/chrom_vcfs/{chrom}.vcf.gz.tbi",
        ref="resources/Homo_sapiens_assembly38.fasta",
        fai="resources/Homo_sapiens_assembly38.fasta.fai",
    output:
        vcf=temp(
            "results/HaplotypeCaller/genotyped/chrom_vcfs/{chrom}_split_multiallelics.vcf.gz"
        ),
        index=temp(
            "results/HaplotypeCaller/genotyped/chrom_vcfs/{chrom}_split_multiallelics.vcf.gz.tbi"
        ),
    benchmark:
        "results/performance_benchmarks/split_multiallelics/{chrom}/benchmarks.tsv"
    threads: config["bcftools"]["threads"]
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
        vcf="results/HaplotypeCaller/genotyped/chrom_vcfs/{chrom}_split_multiallelics.vcf.gz",
        index=(
            "results/HaplotypeCaller/genotyped/chrom_vcfs/{chrom}_split_multiallelics.vcf.gz.tbi"
        ),
    output:
        vcf=temp("results/HaplotypeCaller/filtered/{chrom}_snps.vcf.gz"),
        index=temp("results/HaplotypeCaller/filtered/{chrom}_snps.vcf.gz.tbi"),
    params:
        t=tempDir,
        java_opts=utils.allow_blanks(config["selectVariants"]["java_opts"]),
    benchmark:
        "results/performance_benchmarks/split_snps/{chrom}/benchmarks.tsv"
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["selectVariants"]["memory"],
        xmx=lambda wildcards, attempt: attempt * config["selectVariants"]["xmx"],
    shell:
        'export _JAVA_OPTIONS="" && '
        'gatk --java-options "-Xmx{resources.xmx}m {params.java_opts}" SelectVariants '
        "--tmp-dir {params.t} "
        "-V {input.vcf} "
        "--select-type SNP "
        "-O {output.vcf}"


rule split_indels:
    """
    """
    input:
        vcf="results/HaplotypeCaller/genotyped/chrom_vcfs/{chrom}_split_multiallelics.vcf.gz",
        index=(
            "results/HaplotypeCaller/genotyped/chrom_vcfs/{chrom}_split_multiallelics.vcf.gz.tbi"
        ),
    output:
        vcf=temp("results/HaplotypeCaller/filtered/{chrom}_indels.vcf.gz"),
        index=temp("results/HaplotypeCaller/filtered/{chrom}_indels.vcf.gz.tbi"),
    params:
        t=tempDir,
        java_opts=utils.allow_blanks(config["selectVariants"]["java_opts"]),
    benchmark:
        "results/performance_benchmarks/split_indels/{chrom}/benchmarks.tsv"
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["selectVariants"]["memory"],
        xmx=lambda wildcards, attempt: attempt * config["selectVariants"]["xmx"],
    shell:
        'export _JAVA_OPTIONS="" && '
        'gatk --java-options "-Xmx{resources.xmx}m {params.java_opts}" SelectVariants '
        "--tmp-dir {params.t} "
        "-V {input.vcf} "
        "--select-type INDEL "
        "-O {output.vcf}"


rule hard_filter_snps:
    """
    """
    input:
        vcf="results/HaplotypeCaller/filtered/{chrom}_snps.vcf.gz",
        index="results/HaplotypeCaller/filtered/{chrom}_snps.vcf.gz.tbi",
    output:
        vcf=temp("results/HaplotypeCaller/filtered/{chrom}_snps.hardfiltered.vcf.gz"),
        index=temp("results/HaplotypeCaller/filtered/{chrom}_snps.hardfiltered.vcf.gz.tbi"),
    params:
        t=tempDir,
        java_opts=utils.allow_blanks(config["variantFiltration"]["java_opts"]),
    benchmark:
        "results/performance_benchmarks/hard_filter_snps/{chrom}/benchmarks.tsv"
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["variantFiltration"]["memory"],
        xmx=lambda wildcards, attempt: attempt * config["variantFiltration"]["xmx"],
    shell:
        'export _JAVA_OPTIONS="" && '
        'gatk --java-options "-Xmx{resources.xmx}m {params.java_opts}" VariantFiltration '
        "--tmp-dir {params.t} "
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
        vcf="results/HaplotypeCaller/filtered/{chrom}_indels.vcf.gz",
        index="results/HaplotypeCaller/filtered/{chrom}_indels.vcf.gz.tbi",
    output:
        vcf=temp("results/HaplotypeCaller/filtered/{chrom}_indels.hardfiltered.vcf.gz"),
        index=temp("results/HaplotypeCaller/filtered/{chrom}_indels.hardfiltered.vcf.gz.tbi"),
    params:
        t=tempDir,
        java_opts=utils.allow_blanks(config["variantFiltration"]["java_opts"]),
    benchmark:
        "results/performance_benchmarks/hard_filter_indels/{chrom}/benchmarks.tsv"
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["variantFiltration"]["memory"],
        xmx=lambda wildcards, attempt: attempt * config["variantFiltration"]["xmx"],
    shell:
        'export _JAVA_OPTIONS="" && '
        'gatk --java-options "-Xmx{resources.xmx}m {params.java_opts}" VariantFiltration '
        "--tmp-dir {params.t} "
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

if full:

    rule contamination_check:
        input:
            vcf=expand(
                "results/HaplotypeCaller/filtered/{region}_snps.hardfiltered.vcf.gz",
                region=config["verifyBamID"]["region"],
            ),
            i1=expand(
                "results/HaplotypeCaller/filtered/{region}_snps.hardfiltered.vcf.gz.tbi",
                region=config["verifyBamID"]["region"],
            ),
            bam="results/bqsr/{sample}.bam",
            i2="results/bqsr/{sample}.bai",
        output:
            "results/qc/contamination_check/{region}/{sample}.selfSM",
            "results/qc/contamination_check/{region}/{sample}.selfRG",
            "results/qc/contamination_check/{region}/{sample}.log",
            "results/qc/contamination_check/{region}/{sample}.depthSM",
            "results/qc/contamination_check/{region}/{sample}.depthRG",
            "results/qc/contamination_check/{region}/{sample}.bestSM",
            "results/qc/contamination_check/{region}/{sample}.bestRG",
        params:
            d="results/qc/contamination_check/{region}/{sample}",
        benchmark:
            "results/performance_benchmarks/contamination_check/{region}/{sample}.tsv"
        conda:
            "../envs/verifybamid.yaml"
        resources:
            mem_mb=lambda wildcards, attempt: attempt * config["verifyBamID"]["memory"],
        shell:
            "verifyBamID "
            "--best "
            "--vcf {input.vcf} "
            "--bam {input.bam} "
            "--out {params.d}"

    rule summarize_contam_check:
        """This pulls the FREEMIX and depth estimate from the verifyBamID output"""
        input:
            expand(
                "results/qc/contamination_check/{region}/{sample}.selfSM",
                sample=SAMPLES,
                region=config["verifyBamID"]["region"],
            ),
        output:
            "results/qc/contamination_check/summary.txt",
        params:
            "results/qc/contamination_check/{region}",
        shell:
            'grep -v "^#" {params}*selfSM  > {output}'


rule merge_calls:
    """
    Opted for bcftools here, just because it's easy to multithread,
    sort, zip, etc. all in one data stream.
    """
    input:
        snps=expand(
            "results/HaplotypeCaller/filtered/{chrom}_snps.hardfiltered.vcf.gz", chrom=chromList
        ),
        s_index=expand(
            "results/HaplotypeCaller/filtered/{chrom}_snps.hardfiltered.vcf.gz.tbi", chrom=chromList
        ),
        indels=expand(
            "results/HaplotypeCaller/filtered/{chrom}_indels.hardfiltered.vcf.gz", chrom=chromList
        ),
        i_index=expand(
            "results/HaplotypeCaller/filtered/{chrom}_indels.hardfiltered.vcf.gz.tbi",
            chrom=chromList,
        ),
    output:
        vcf="results/HaplotypeCaller/filtered/{chrom}.hardfiltered.vcf.gz",
        index="results/HaplotypeCaller/filtered/{chrom}.hardfiltered.vcf.gz.tbi",
    params:
        t=tempDir,
    benchmark:
        "results/performance_benchmarks/merge_calls/{chrom}/benchmarks.tsv"
    conda:
        "../envs/bcftools_tabix.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["bcftools"]["memory"],
    shell:
        "bcftools concat -a -Ov {input.snps} {input.indels} | "
        "bcftools sort -T {params.t} -Oz -o {output.vcf} && "
        "tabix -p vcf {output.vcf}"


rule picard_metrics:
    """
    """
    input:
        calls=expand(
            "results/HaplotypeCaller/filtered/{chrom}.hardfiltered.vcf.gz", chrom=chromList
        ),
        index=expand(
            "results/HaplotypeCaller/filtered/{chrom}.hardfiltered.vcf.gz.tbi", chrom=chromList
        ),
        dbsnp="resources/Homo_sapiens_assembly38.dbsnp138.vcf",
        d="resources/Homo_sapiens_assembly38.dict",
    output:
        "results/HaplotypeCaller/filtered/{chrom}.variant_calling_detail_metrics",
        "results/HaplotypeCaller/filtered/{chrom}.variant_calling_summary_metrics",
    benchmark:
        "results/performance_benchmarks/picard_metrics/{chrom}/benchmarks.tsv"
    params:
        d="results/HaplotypeCaller/filtered/{chrom}",
        t=tempDir,
        java_opts=utils.allow_blanks(config["picardCollectVariantCallingMetrics"]["java_opts"]),
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=(
            lambda wildcards, attempt: attempt
            * config["picardCollectVariantCallingMetrics"]["memory"]
        ),
        xmx=(
            lambda wildcards, attempt: attempt * config["picardCollectVariantCallingMetrics"]["xmx"]
        ),
    shell:
        'export _JAVA_OPTIONS="" && '
        'gatk --java-options "-Xmx{resources.xmx}m {params.java_opts}" CollectVariantCallingMetrics '
        "--TMP_DIR {params.t} "
        "-I {input.calls} "
        "--DBSNP {input.dbsnp} "
        "-SD {input.d} "
        "-O {params.d}"


if config["combine_vcfs"]:

    rule HC_concat_vcfs_bcftools:
        """If specified as 'True' in the config, this rule will combine the per-chromosome joint-called and filtered VCFs into one multi-sample all-region project-level VCF. Default value is 'False'."""
        input:
            vcfList=expand(
                "results/HaplotypeCaller/filtered/{chrom}.hardfiltered.vcf.gz", chrom=chromList
            ),
            indexList=expand(
                "results/HaplotypeCaller/filtered/{chrom}.hardfiltered.vcf.gz.tbi", chrom=chromList
            ),
        output:
            projectVCF="results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz",
            projectIndex="results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz.tbi",
        params:
            t=tempDir,
        benchmark:
            "results/performance_benchmarks/HC_combine_chrom_vcfs/benchmarks.tsv"
        conda:
            "../envs/bcftools_tabix.yaml"
        resources:
            mem_mb=lambda wildcards, attempt: attempt * config["bcftools"]["memory"],
            queue=config["compute_queue"],
        shell:
            "bcftools concat -a -Ou {input.vcfList} | "
            "bcftools sort -T {params.t} -Oz -o {output.projectVCF} && "
            "tabix -p vcf {output.projectVCF}"
