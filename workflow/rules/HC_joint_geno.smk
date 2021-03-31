if jointgeno:

    rule symlink_gvcfs:
        """"""
        input:
            gvcf=utils.get_gvcf,
            index=utils.get_gvcf_index,
        output:
            gvcf="results/HaplotypeCaller/called/{sample}_all_chroms.g.vcf.gz",
            index="results/HaplotypeCaller/called/{sample}_all_chroms.g.vcf.gz.tbi",
        shell:
            "ln -s {input.gvcf} {output.gvcf} && "
            "ln -s {input.index} {output.index}"


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
    seems to solve this problem.  Adding the directory as an input ensures
    that a failed job will clean up after itself, so resumption of the
    pipeline shouldn't be an issue.

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

    For the batch and cache options - see DBImport docs.  May want to enable
    caching for exome, due to many intervals.  May want to tweak batch size
    based on sample size.
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
        o1="results/HaplotypeCaller/DBImport/{chrom}/vcfheader.vcf",
        o2="results/HaplotypeCaller/DBImport/{chrom}/vidmap.json",
        o3="results/HaplotypeCaller/DBImport/{chrom}/callset.json",
        o4="results/HaplotypeCaller/DBImport/{chrom}/__tiledb_workspace.tdb",
        d=directory("results/HaplotypeCaller/DBImport/{chrom}"),
    benchmark:
        "results/performance_benchmarks/HC_consolidate_gvcfs/{chrom}.tsv"
    params:
        interval="{chrom}",
        db="results/HaplotypeCaller/DBImport/{chrom}",
        t=tempDir,
        java_opts=config["genomicsDBImport"]["java_opts"],
        batch_size=config["genomicsDBImport"]["batch_size"],
        reader_threads=config["genomicsDBImport"]["reader_threads"],
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=config["genomicsDBImport"]["memory"],
    shell:
        'export _JAVA_OPTIONS="" && '
        "rm -r {params.db} && "
        'gatk --java-options "{params.java_opts}" GenomicsDBImport '
        "--batch-size {params.batch_size} "
        "--disable-bam-index-caching "
        "--sample-name-map {input.sampleMap} "
        "--genomicsdb-workspace-path {params.db} "
        "-L {params.interval} "
        "--tmp-dir {params.t} "
        "--reader-threads {params.reader_threads} "
        "--genomicsdb-shared-posixfs-optimizations"


rule HC_genotype_gvcfs:
    """Joint genotyping."""
    input:
        r="resources/Homo_sapiens_assembly38.fasta",
        bed=bed,
        o1="results/HaplotypeCaller/DBImport/{chrom}/vcfheader.vcf",
        o2="results/HaplotypeCaller/DBImport/{chrom}/vidmap.json",
        o3="results/HaplotypeCaller/DBImport/{chrom}/callset.json",
        o4="results/HaplotypeCaller/DBImport/{chrom}/__tiledb_workspace.tdb",
        o5=utils.get_DBImport_path1,
        o6=utils.get_DBImport_path2,
    output:
        vcf=temp("results/HaplotypeCaller/genotyped/{chrom}.vcf.gz"),
        idx=temp("results/HaplotypeCaller/genotyped/{chrom}.vcf.gz.tbi"),
    benchmark:
        "results/performance_benchmarks/HC_genotype_gvcfs/{chrom}.tsv"
    params:
        db="results/HaplotypeCaller/DBImport/{chrom}",
        t=tempDir,
        java_opts=config["genotypeGVCFs"]["java_opts"],
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=config["genotypeGVCFs"]["memory"],
    shell:
        'export _JAVA_OPTIONS="" && '
        'gatk --java-options "{params.java_opts}" GenotypeGVCFs '
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
    resources:
        mem_mb=config["bcftools"]["memory"],
    shell:
        "bcftools concat -a {input.vcfList} -Ou | "
        "bcftools sort -T {params.t} -Oz -o {output.projectVCF} && "
        "tabix -p vcf {output.projectVCF}"
