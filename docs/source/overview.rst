Overview
===============

This workflow was designed by the Genomics & Data Science team (GDS) at 54gene and is used to analyze paired-end short-read germline whole-genome sequencing data.  This pipeline is designed to first be deployed in small batches (e.g. per flow cell), starting with FASTQs and resulting in gVCFs and a small batch joint-called VCF.  A second run of the pipeline can receive a larger batch of gVCFs (e.g. gVCFs derived from many flow cells), and generates a large batch joint-called VCF.  The workflow, which is designed to support reproducible bioinformatics, is written in `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ and is platform-agnostic.  All dependencies are installed by the pipeline as-needed using `conda <https://docs.conda.io/en/latest/>`_.  Development and testing has been predominantly on AWS' `ParallelCluster <https://aws.amazon.com/hpc/parallelcluster/>`_ using `Amazon Linux <https://aws.amazon.com/amazon-linux-2/>`_ using Snakemake version 7.8.2.

Features:

- Read filtering and trimming
- Read alignment, deduplication, and BQSR
- Variant calling and filtering
- Joint-genotyping
- Sex discordance and relatedness assessment
- Generate MultiQC reports

To install the latest release, type::

    git clone https://gitlab.com/data-analysis5/dna-sequencing/54gene-wgs-germline.git

Inputs
------

The pipeline requires the following inputs:

- A headerless, whitespace delimited ``manifest.txt`` file with sample names and paths (columns dependent on the run-mode)
- Config file with the run-mode specified and other pipeline parameters configured (see default config provided in ``config/config.yaml``)
- A tab-delimited ``intervals.tsv`` file with names of intervals and paths to region (BED) files of the genome you want to parallelize the variant calling and joint-calling steps by (i.e. 50 BED files each with a small region of the genome to parallelize by)
- A tab-delimited ``sex_linker.tsv`` file with the sample names in one column and sex in the other to identify discordances in reported vs. inferred sex
- A ``multiqc.yaml`` config file for generating MultiQC reports (provided for you)

Outputs
-------

Depending on which run-mode you have set, you will be able to generate:

- A hard-filtered, multi-sample joint-called VCF in ``full`` and ``joint_genotyping`` mode
- Per-sample gVCFs for all regions of the genome for future joint-calling in ``full`` mode
- Deduplicated and post-BQSR BAM files in ``full`` mode
- Various QC metrics (e.g. FastQC, MultiQC, bcftools stats) in all three modes


See the :doc:`installation`, :doc:`execution`, and :doc:`configuration` for details on setting up and running the pipeline.
