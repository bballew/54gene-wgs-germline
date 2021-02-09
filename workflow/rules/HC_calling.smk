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
    benchmark:
        "results/performance_benchmarks/HC_call_variants/{sample}_{chrom}.tsv"
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=8000
    shell:
        'gatk HaplotypeCaller '
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
        lambda wildcards, input: " -I ".join(input.vcfList),
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=8000
    shell:
        'gatk GatherVcfs -I {params} -O {output}'


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


rule HC_create_each_sample_map_file:
    """Create sample file to be read by GenomicsDBImport.

    Can pass all vcfs on the command line, each with a -V before it,
    but with enough samples you will hit the character limit.  Instead,
    create a sample map file and read from it.

    sample map is formatted as follows: sample_name--tab--path_to_sample_vcf per line

    See GATK's documentation for additional details.
    """
    input:
        gvcf="results/HaplotypeCaller/called/{sample}_all_chroms.g.vcf.gz",
        index="results/HaplotypeCaller/called/{sample}_all_chroms.g.vcf.gz.tbi",
    output:
        temp("results/HaplotypeCaller/DBImport/{sample}.map"),
    benchmark:
        "results/performance_benchmarks/HC_create_each_sample_map_file/{sample}.tsv"
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "n=$(bcftools query -l {input.gvcf});"
        'echo "${{n}}\t{input.gvcf}" > {output}'


rule HC_create_cohort_map_file:
    """Combine per-sample map file into one cohort map file."""
    input:
        expand("results/HaplotypeCaller/DBImport/{sample}.map", sample=SAMPLES),
    output:
        temp("results/HaplotypeCaller/DBImport/cohort.sample_map"),
    benchmark:
        "results/performance_benchmarks/HC_create_cohort_map_file/benchmarks.tsv"
    params:
        mapDir="results/HaplotypeCaller/DBImport/",
    shell:
        "for i in {params.mapDir}*.map; do cat $i >> {output}; done"


rule HC_consolidate_gvcfs:
    """Split DBimport databases by chromosome.

    The output of this step includes some files that are in subdirectories
    with unpredictable names.  I've attempted to include them in the input
    of the next rule via glob in the functions below (o5 and o6 in rule
    HC_genotypeGVCFs) but this is still being tested.

    Note that DBImport requires a new or empty directory for
    --genomicsdb-workspace-path.  Possible issue when resuming a pipeline.
    Snakemake's implicit directory management results in an error when
    DBImport finds that the workspace path already exists (even though it
    seems empty?).  Removing the directory just prior to running DBImport
    seems to solve this problem, but will be problematic on resuming the
    pipeline.  Added a datetime stamp (dt) to the dir to help address this.
    However, this datetime stamp will need to be overridden if the pipeline
    is stopped and then resumed prior to completion of rule HC_genotypeGVCFs.

    What exactly is GenomicsDB for?  GenotypeGVCFs can only take one single
    input.  GenomicsDB is a method of consolidating gvcfs across samples,
    to provide the input for genotyping.  The alternative is to use
    CombineGVCFs.  See
    https://software.broadinstitute.org/gatk/documentation/article?id=11813
    for more details.

    As documented here (https://software.broadinstitute.org/gatk/documentation/
    tooldocs/4.0.11.0/org_broadinstitute_hellbender_tools_genomicsdb_
    GenomicsDBImport.php), DBImport can take up a lot of /tmp space.  GATK
    recommends using --tmp-dir to redirect to a large temp space.  Annoyingly,
    this tool appears unable to implicitly create temp subdirectories, so you
    have to explicitly include temp dir creation in the shell section.
    Note that GenotypeGVCFs below does not suffer from this directory creation
    issue.
    """
    input:
        sampleMap="results/HaplotypeCaller/DBImport/cohort.sample_map",
        gvcfList=expand(
            "results/HaplotypeCaller/called/{sample}_all_chroms.g.vcf.gz", sample=SAMPLES
        ),
        indexList=expand(
            "results/HaplotypeCaller/called/{sample}_all_chroms.g.vcf.gz.tbi", sample=SAMPLES
        ),
    output:
        t=temp(directory(tempDir + "HC_DBImport/{chrom}/")),
        o1="results/HaplotypeCaller/DBImport/{chrom}/vcfheader.vcf",
        o2="results/HaplotypeCaller/DBImport/{chrom}/vidmap.json",
        o3="results/HaplotypeCaller/DBImport/{chrom}/callset.json",
        o4="results/HaplotypeCaller/DBImport/{chrom}/__tiledb_workspace.tdb",
    benchmark:
        "results/performance_benchmarks/HC_consolidate_gvcfs/{chrom}.tsv"
    params:
        interval="{chrom}",
        db="results/HaplotypeCaller/DBImport/{chrom}",
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=23000
    shell:
        "mkdir -p {output.t} && "
        "rm -r {params.db} && "
        'gatk --java-options "-Xmx20G" GenomicsDBImport '
        "--sample-name-map {input.sampleMap} "
        "--genomicsdb-workspace-path {params.db} "
        "-L {params.interval} "
        "--tmp-dir {output.t}"


rule HC_genotype_gvcfs:
    """Joint genotyping."""
    input:
        r="resources/Homo_sapiens_assembly38.fasta",
        bed=bed,
        o1="results/HaplotypeCaller/DBImport/{chrom}/vcfheader.vcf",
        o2="results/HaplotypeCaller/DBImport/{chrom}/vidmap.json",
        o3="results/HaplotypeCaller/DBImport/{chrom}/callset.json",
        o4="results/HaplotypeCaller/DBImport/{chrom}/__tiledb_workspace.tdb",
        o5=get_DBImport_path1,
        o6=get_DBImport_path2,
    output:
        vcf=temp("results/HaplotypeCaller/genotyped/{chrom}.vcf.gz"),
        idx=temp("results/HaplotypeCaller/genotyped/{chrom}.vcf.gz.tbi"),
    benchmark:
        "results/performance_benchmarks/HC_genotype_gvcfs/{chrom}.tsv"
    params:
        db="results/HaplotypeCaller/DBImport/{chrom}",
        t=tempDir,
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=8000
    shell:
        'gatk GenotypeGVCFs '
        "-R {input.r} "
        "-V gendb://{params.db} "
        "-O {output.vcf} "
        "--tmp-dir {params.t} "
        "-stand-call-conf 30 "
        "-G StandardAnnotation "
        "-G StandardHCAnnotation"


rule HC_concat_vcfs_bcftools:
    """Combine per-chromosome joint-called multi-sample VCFs."""
    input:
        vcfList=expand("results/HaplotypeCaller/genotyped/{chrom}.vcf.gz", chrom=chromList),
        indexList=expand("results/HaplotypeCaller/genotyped/{chrom}.vcf.gz.tbi", chrom=chromList),
    output:
        projectVCF=protected("results/HaplotypeCaller/genotyped/HC_variants.vcf.gz"),
        idx=protected("results/HaplotypeCaller/genotyped/HC_variants.vcf.gz.tbi"),
    params:
        t=tempDir,
    benchmark:
        "results/performance_benchmarks/HC_concat_vcfs_bcftools/benchmarks.tsv"
    conda:
        "../envs/bcftools_tabix.yaml"
    shell:
        "bcftools concat -a {input.vcfList} -Ou | "
        "bcftools sort -T {params.t} -Oz -o {output.projectVCF} && "
        "tabix -p vcf {output.projectVCF}"

