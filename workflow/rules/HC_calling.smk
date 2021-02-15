rule split_bed_file:
    """Separates bed regions by chromosome.

    For DV:
    If you're not assigning a number of shards by which to divide
    and parallelize, then the pipeline will parallelize by chrom.
    To do this, we take the bed file (e.g. exome capture region)
    and split the regions by chromosome.  Subsequent steps are run
    concurrently on each of the single-chromosome bed files.

    For GATK:
    HaplotypeCaller can't be parallelized per task (e.g. threads),
    so must be run over sub-regions if you want parallelization.
    **Do we want to use the old 4000-region bed file, or is by-chrom
    sufficient?

    Note that grep exits with 0 if a match is found, 1 if no match,
    and 2 if error.  Snakemake looks for exit codes of 0 to determine
    that a job finished successfully.  No match is an acceptable outcome
    here, so the shell command below should allow match or no match.
    """
    input:
        bed,
    output:
        temp("results/split_regions/{chrom}.bed"),
    benchmark:
        "results/performance_benchmarks/split_bed_file/{chrom}.tsv"
    shell:
        'grep "^{wildcards.chrom}[[:space:]]" {input} > {output};'
        "if [ $? -le 1 ]; then exit 0; else exit 1; fi"


rule HC_call_variants:
    """Call gVCFs with GATK4.
    Runs over each chrom in parallel.
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
        bed="results/split_regions/{chrom}.bed",
        bam="results/bqsr/{sample}.bam",
        bai="results/bqsr/{sample}.bam.bai",
    output:
        gvcf=temp("results/HaplotypeCaller/called/{chrom}/{sample}.g.vcf"),
        idx=temp("results/HaplotypeCaller/called/{chrom}/{sample}.g.vcf.idx"),
    params:
        t=tempDir,
    benchmark:
        "results/performance_benchmarks/HC_call_variants/{sample}_{chrom}.tsv"
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=16000,
    shell:
        'gatk --java-options "-Xmx4g" HaplotypeCaller '
        "--tmp-dir {params.t} "
        "-R {input.r} "
        "-I {input.bam} "
        "-ERC GVCF "
        "-L {input.bed} "
        "-O {output.gvcf} "
        "-G StandardAnnotation "
        "-G StandardHCAnnotation"


rule HC_compress_gvcfs:
    """Zip and index gVCFs."""
    input:
        gvcf="results/HaplotypeCaller/called/{chrom}/{sample}.g.vcf",
        idx="results/HaplotypeCaller/called/{chrom}/{sample}.g.vcf.idx",
    output:
        temp("results/HaplotypeCaller/called/{chrom}/{sample}.g.vcf.gz"),
        temp("results/HaplotypeCaller/called/{chrom}/{sample}.g.vcf.gz.tbi"),
    benchmark:
        "results/performance_benchmarks/HC_compress_gvcfs/{sample}_{chrom}.tsv"
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "bgzip {input.gvcf} && "
        "tabix -p vcf {input.gvcf}.gz"


rule HC_concat_gvcfs:
    """Combine per-chromosome gvcfs; maintain one per sample.
    Not clear whether it would be fastest to concat per-chrom gvcfs and
    then genotype, or genotype and then concat.
    """
    input:
        vcfList=expand(
            "results/HaplotypeCaller/called/{chrom}/{{sample}}.g.vcf.gz", chrom=chromList
        ),
        indexList=expand(
            "results/HaplotypeCaller/called/{chrom}/{{sample}}.g.vcf.gz.tbi", chrom=chromList
        ),
    output:
        "results/HaplotypeCaller/called/{sample}_all_chroms.g.vcf.gz",
    benchmark:
        "results/performance_benchmarks/HC_concat_gvcfs/{sample}.tsv"
    params:
        l=lambda wildcards, input: " -I ".join(input.vcfList),
        t=tempDir,
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=8000,
    shell:
        'gatk --java-options "-Xmx4g" GatherVcfs '
        "--tmp-dir {params.t} "
        "-I {params.l} "
        "-O {output}"


rule HC_index_gvcf:
    """index per-sample gvcfs."""
    input:
        "results/HaplotypeCaller/called/{sample}_all_chroms.g.vcf.gz",
    output:
        "results/HaplotypeCaller/called/{sample}_all_chroms.g.vcf.gz.tbi",
    benchmark:
        "results/performance_benchmarks/HC_index_gvcf/{sample}.tsv"
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "tabix -p vcf {input}"
