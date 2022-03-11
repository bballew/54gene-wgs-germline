if jointgeno:

    localrules:
        symlink_gvcfs,

    rule symlink_gvcfs:
        """"""
        input:
            gvcf=utils.get_gvcf,
            index=utils.get_gvcf_index,
        output:
            gvcf="results/HaplotypeCaller/called/{sample}_all_regions.g.vcf.gz",
            index="results/HaplotypeCaller/called/{sample}_all_regions.g.vcf.gz.tbi",
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
        gvcf="results/HaplotypeCaller/called/{sample}_all_regions.g.vcf.gz",
        index="results/HaplotypeCaller/called/{sample}_all_regions.g.vcf.gz.tbi",
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
    """Split DBimport databases by intervals for faster import.
    A database will be created for each interval which will allow for faster
    iteration over the intervals during the joint genotyping step that proceeds.

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
            "results/HaplotypeCaller/called/{sample}_all_regions.g.vcf.gz", sample=SAMPLES
        ),
        indexList=expand(
            "results/HaplotypeCaller/called/{sample}_all_regions.g.vcf.gz.tbi", sample=SAMPLES
        ),
        intervalfile=lambda wildcards: intervals_df.loc[wildcards.interval, "file_path"],
    output:
        o1="results/HaplotypeCaller/DBImport/{interval}/vcfheader.vcf",
        o2="results/HaplotypeCaller/DBImport/{interval}/vidmap.json",
        o3="results/HaplotypeCaller/DBImport/{interval}/callset.json",
        o4="results/HaplotypeCaller/DBImport/{interval}/__tiledb_workspace.tdb",
        d=directory("results/HaplotypeCaller/DBImport/{interval}"),
    benchmark:
        "results/performance_benchmarks/HC_consolidate_gvcfs/{interval}.tsv"
    params:
        db="results/HaplotypeCaller/DBImport/{interval}",
        t=tempDir,
        java_opts=utils.allow_blanks(config["genomicsDBImport"]["java_opts"]),
        batch_size=config["genomicsDBImport"]["batch_size"],
        reader_threads=config["genomicsDBImport"]["reader_threads"],
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["genomicsDBImport"]["memory"],
        xmx=lambda wildcards, attempt: attempt * config["genomicsDBImport"]["xmx"],
        queue=config["memory_queue"],
    shell:
        'export _JAVA_OPTIONS="" && '
        "rm -r {params.db} && "
        'gatk --java-options "-Xmx{resources.xmx}m {params.java_opts}" GenomicsDBImport '
        "--batch-size {params.batch_size} "
        "--disable-bam-index-caching "
        "--sample-name-map {input.sampleMap} "
        "--genomicsdb-workspace-path {params.db} "
        "-L {input.intervalfile} "
        "--tmp-dir {params.t} "
        "--reader-threads {params.reader_threads} "
        "--genomicsdb-shared-posixfs-optimizations"


rule HC_genotype_gvcfs:
    """Joint genotyping."""
    input:
        r="resources/Homo_sapiens_assembly38.fasta",
        f="resources/Homo_sapiens_assembly38.fasta.fai",
        d="resources/Homo_sapiens_assembly38.dict",
        o1="results/HaplotypeCaller/DBImport/{interval}/vcfheader.vcf",
        o2="results/HaplotypeCaller/DBImport/{interval}/vidmap.json",
        o3="results/HaplotypeCaller/DBImport/{interval}/callset.json",
        o4="results/HaplotypeCaller/DBImport/{interval}/__tiledb_workspace.tdb",
        o5=utils.get_DBImport_path1,
        o6=utils.get_DBImport_path2,
    output:
        vcf=temp("results/HaplotypeCaller/genotyped/{interval}.vcf.gz"),
        idx=temp("results/HaplotypeCaller/genotyped/{interval}.vcf.gz.tbi"),
    benchmark:
        "results/performance_benchmarks/HC_genotype_gvcfs/{interval}.tsv"
    params:
        db="results/HaplotypeCaller/DBImport/{interval}",
        t=tempDir,
        m=config["genotypeGVCFs"]["max_alt_alleles"],
        g=config["genotypeGVCFs"]["genomicsdb_max_alt_alleles"],
        c=config["genotypeGVCFs"]["max_genotype_count"],
        java_opts=utils.allow_blanks(config["genotypeGVCFs"]["java_opts"]),
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["genotypeGVCFs"]["memory"],
        xmx=lambda wildcards, attempt: attempt * config["genotypeGVCFs"]["xmx"],
        queue=config["memory_queue"],
    shell:
        'export _JAVA_OPTIONS="" && '
        'gatk --java-options "-Xmx{resources.xmx}m {params.java_opts}" GenotypeGVCFs '
        "-R {input.r} "
        "-V gendb://{params.db} "
        "-O {output.vcf} "
        "--tmp-dir {params.t} "
        "--max-alternate-alleles {params.m} "
        "--genomicsdb-max-alternate-alleles {params.g} "
        "--max-genotype-count {params.c} "
        "-stand-call-conf 30 "
        "-G StandardAnnotation "
        "-G StandardHCAnnotation"


rule HC_split_intervals_chrom:
    """Separate the VCFs joint-called over intervals by chromosome."""
    input:
        vcfList=expand("results/HaplotypeCaller/genotyped/{interval}.vcf.gz", interval=intervalList),
        indexList=expand(
            "results/HaplotypeCaller/genotyped/{interval}.vcf.gz.tbi", interval=intervalList
        ),
    output:
        vcf=temp("results/HaplotypeCaller/split_vcfs/{interval}/{chrom}.vcf.gz"),
        index=temp("results/HaplotypeCaller/split_vcfs/{interval}/{chrom}.vcf.gz.tbi"),
    benchmark:
        "results/performance_benchmarks/HC_split_intervals_chrom/{interval}/{chrom}/benchmarks.tsv"
    conda:
        "../envs/bcftools_tabix.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["bcftools"]["memory"],
        queue=config["compute_queue"],
    shell:
        "bcftools view -r {wildcards.chrom} -Oz -o {output.vcf} {input.vcfList} && "
        "tabix -p vcf {output.vcf}"


rule HC_combine_chrom_vcfs:
    """Combine the VCFs for each chromosome from each interval to generate multi-sample per-chromosome VCFs."""
    input:
        vcfList=expand(
            "results/HaplotypeCaller/split_vcfs/{interval}/{{chrom}}.vcf.gz",
            interval=intervalList,
            chrom=chromList,
        ),
        indexList=expand(
            "results/HaplotypeCaller/split_vcfs/{interval}/{{chrom}}.vcf.gz.tbi",
            interval=intervalList,
            chrom=chromList,
        ),
    output:
        vcf="results/HaplotypeCaller/genotyped/chrom_vcfs/{chrom}.vcf.gz",
        index="results/HaplotypeCaller/genotyped/chrom_vcfs/{chrom}.vcf.gz.tbi",
    conda:
        "../envs/bcftools_tabix.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["bcftools"]["memory"],
        queue=config["compute_queue"],
    shell:
        "bcftools concat -a -Ou {input.vcfList} | "
        "bcftools sort -Oz -o {output.vcf} && "
        "tabix -p vcf {output.vcf}"
