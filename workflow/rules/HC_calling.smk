
rule HC_call_variants:
    """Call gVCFs with GATK4.
    Runs over each interval in a list of 50 (obtained from resources provided by GATK
    in parallel. This is based on the GRCH38 build. 
    """
    input:
        r="resources/Homo_sapiens_assembly38.fasta",
        d="resources/Homo_sapiens_assembly38.dict",
        fai="resources/Homo_sapiens_assembly38.fasta.fai",
        amb="resources/Homo_sapiens_assembly38.fasta.64.amb",
        ann="resources/Homo_sapiens_assembly38.fasta.64.ann",
        bwt="resources/Homo_sapiens_assembly38.fasta.64.bwt",
        pac="resources/Homo_sapiens_assembly38.fasta.64.pac",
        sa="resources/Homo_sapiens_assembly38.fasta.64.sa",
        # bed="results/split_regions/{chrom}.bed",
        bam="results/bqsr/{sample}.bam",
        bai="results/bqsr/{sample}.bai",
        intervals="resources/scattered_calling_intervals/{directory}_of_50/scattered.interval_list"
    output:
        gvcf=temp("results/HaplotypeCaller/called/interval_{directory}/{sample}.g.vcf"),
        idx=temp("results/HaplotypeCaller/called/interval_{directory}/{sample}.g.vcf.idx"),
    params:
        t=tempDir,
        java_opts=utils.allow_blanks(config["haplotypeCaller"]["java_opts"]),
    benchmark:
        "results/performance_benchmarks/HC_call_variants/{sample}_interval_{directory}.tsv"
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["haplotypeCaller"]["memory"],
        xmx=lambda wildcards, attempt: attempt * config["haplotypeCaller"]["xmx"],
        queue=config["memory_queue"],
    shell:
        'gatk --java-options "-Xmx{resources.xmx}m {params.java_opts}" HaplotypeCaller '
        "--tmp-dir {params.t} "
        "-R {input.r} "
        "-I {input.bam} "
        "-ERC GVCF "
        "-L {input.intervals} "
        "-O {output.gvcf} "
        "-G StandardAnnotation "
        "-G StandardHCAnnotation"


rule HC_compress_gvcfs:
    """Zip and index gVCFs."""
    input:
        gvcf="results/HaplotypeCaller/called/interval_{directory}/{sample}.g.vcf",
        idx="results/HaplotypeCaller/called/interval_{directory}/{sample}.g.vcf.idx",
    output:
        gvcf=temp("results/HaplotypeCaller/called/interval_{directory}/{sample}.g.vcf.gz"),
        tbi=temp("results/HaplotypeCaller/called/interval_{directory}/{sample}.g.vcf.gz.tbi"),
    benchmark:
        "results/performance_benchmarks/HC_compress_gvcfs/{sample}_interval_{directory}.tsv"
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "bgzip -c {input.gvcf} > {output.gvcf} && "
        "tabix -p vcf {input.gvcf}.gz"


rule HC_concat_gvcfs:
    """Combine per-interval (instead of per-chromosome) gvcfs; maintain one per sample.
    Here I just rename it to "all_regions" instead.

    Prev: Not clear whether it would be fastest to concat per-chrom gvcfs and
    then genotype, or genotype and then concat.
    """
    input:
        vcfList=expand(
            "results/HaplotypeCaller/called/interval_{directory}/{{sample}}.g.vcf.gz", directory=INTERVALS
        ),
        indexList=expand(
            "results/HaplotypeCaller/called/interval_{directory}/{{sample}}.g.vcf.gz.tbi", directory=INTERVALS
        ),
    output:
        "results/HaplotypeCaller/called/{sample}_all_regions.g.vcf.gz",
    benchmark:
        "results/performance_benchmarks/HC_concat_gvcfs/{sample}.tsv"
    params:
        l=lambda wildcards, input: " -I ".join(input.vcfList),
        t=tempDir,
        java_opts=utils.allow_blanks(config["gatherVcfs"]["java_opts"]),
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["gatherVcfs"]["memory"],
        xmx=lambda wildcards, attempt: attempt * config["gatherVcfs"]["xmx"],
        queue=config["compute_queue"],
    shell:
        'gatk --java-options "-Xmx{resources.xmx}m {params.java_opts}" GatherVcfs '
        "-I {params.l} "
        "-O {output}"


rule HC_index_gvcf:
    """index per-sample gvcfs."""
    input:
        "results/HaplotypeCaller/called/{sample}_all_regions.g.vcf.gz",
    output:
        "results/HaplotypeCaller/called/{sample}_all_regions.g.vcf.gz.tbi",
    benchmark:
        "results/performance_benchmarks/HC_index_gvcf/{sample}.tsv"
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "tabix -p vcf {input}"
