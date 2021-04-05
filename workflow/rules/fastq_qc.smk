rule symlink_fastqs:
    """Create symbolic links with a different naming scheme.
    FastQC doesn't allow you to change the output file names, so
    this is a way to get the sample and readgroup info into the reports
    and use the filenaming convention followed in the rest of the pipeline.

    Note that we're assuming paired end data throughout.
    """
    input:
        r1=utils.get_read1_fastq,
        r2=utils.get_read2_fastq,
    output:
        r1="results/input/{rg}_r1.fastq.gz",
        r2="results/input/{rg}_r2.fastq.gz",
    group:
        "symlink_group"
    benchmark:
        "results/performance_benchmarks/symlink_fastqs/{rg}.tsv"
    shell:
        "ln -s {input.r1} {output.r1} && "
        "ln -s {input.r2} {output.r2}"


rule group_symlinks:
    """
    This is a dummy rule to facilitate grouping of the upstream symlink
    rule.  The group will ensure that all symlink jobs will be submitted
    at once to one node, to save queuing and execution time.
    """
    input:
        r1=expand("results/input/{rg}_r1.fastq.gz", rg=sampleDict.keys()),
        r2=expand("results/input/{rg}_r2.fastq.gz", rg=sampleDict.keys()),
    output:
        temp("results/input/grouped.out"),
    group:
        "symlink_group"
    shell:
        "sleep 10 && touch {output}"


rule fastqc:
    """Generate FastQC reports for all input fastqs."""
    input:
        t="results/input/grouped.out",
        r1="results/input/{rg}_r1.fastq.gz",
        r2="results/input/{rg}_r2.fastq.gz",
    output:
        html1="results/fastqc/{rg}_r1_fastqc.html",
        zip1="results/fastqc/{rg}_r1_fastqc.zip",
        html2="results/fastqc/{rg}_r2_fastqc.html",
        zip2="results/fastqc/{rg}_r2_fastqc.zip",
    benchmark:
        "results/performance_benchmarks/fastqc/{rg}.tsv"
    params:
        t=tempDir,
    threads: config["fastqc"]["threads"]
    conda:
        "../envs/fastqc_multiqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["fastqc"]["memory"],
        batch=concurrent_limit,
    shell:
        "fastqc {input.r1} -d {params.t} --quiet -t {threads} --outdir=results/fastqc/ && "
        "fastqc {input.r2} -d {params.t} --quiet -t {threads} --outdir=results/fastqc/"


rule multiqc:
    """Generate one multiQC report for all input fastqs.
    Should add samtools stats output and possibly others eventually,
    dedup metrics, ...
    """
    input:
        expand("results/fastqc/{rg}_r1_fastqc.zip", rg=sampleDict.keys()),
        expand("results/fastqc/{rg}_r2_fastqc.zip", rg=sampleDict.keys()),
    output:
        "results/multiqc/multiqc.html",
    benchmark:
        "results/performance_benchmarks/multiqc/benchmarks.tsv"
    params:
        outDir="results/multiqc/",
        outName="multiqc.html",
    conda:
        "../envs/fastqc_multiqc.yaml"
    shell:
        "multiqc --force -o {params.outDir} -n {params.outName} {input}"


rule quality_trimming:
    """Quality trimming of read ends.
    May want to tweak params; possibly put in config.  Could remove
    adapter sequences here if we want to move that downstream at
    some point.

    Adapter sequences are from Illumina's TruSeq adapters.
    """
    input:
        r1=utils.get_read1_fastq,
        r2=utils.get_read2_fastq,
    output:
        r1_paired=temp("results/paired_trimmed_reads/{rg}_r1.fastq.gz"),
        r2_paired=temp("results/paired_trimmed_reads/{rg}_r2.fastq.gz"),
        h="results/paired_trimmed_reads/{rg}_fastp.html",
        j="results/paired_trimmed_reads/{rg}_fastp.json",
    benchmark:
        "results/performance_benchmarks/quality_trimming/{rg}.tsv"
    threads: config["fastp"]["threads"]
    conda:
        "../envs/fastp.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["fastp"]["memory"],
        batch=concurrent_limit,
        queue=config["compute_queue"],
    shell:
        "fastp -i {input.r1} -I {input.r2} "
        "-w {threads} "
        "-o {output.r1_paired} -O {output.r2_paired} "
        "-h {output.h} -j {output.j} "
        "--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA "
        "--adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT "
        "-l 36"
