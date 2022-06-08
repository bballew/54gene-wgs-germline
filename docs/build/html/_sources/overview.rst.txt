Overview
===============

This workflow was designed by the Genomics & Data Science team at 54gene and is used to analyze germline whole-genome sequencing data in the form of either FASTQs or gVCFs. This pipeline emits a joint-called multi-sample VCF. It is currently optimized to be run on HPC infrastructure and was developed and tested on AWS' `ParallelCluster <https://aws.amazon.com/hpc/parallelcluster/>`_ and Snakemake version 6.15.5.

Features:

- Read filtering and trimming
- Read alignment
- Variant calling and filtering
- Joint-genotyping
- Check sex and relatedness
- Generate a multiQC report

To install the latest release, type::

    git clone https://gitlab.com/data-analysis5/dna-sequencing/54gene-wgs-germline.git

Inputs
------

The pipeline requires the following inputs:

- A headerless, whitespace delimited ``manifest.txt`` file with sample names and paths (columns dependent on the run-mode)
- Config file with the run-mode specified and other pipeline parameters configured (see default config provided in ``config/config.yaml``)
- A tab-delimited ``intervals.tsv`` file with names of intervals and paths to region (BED) files of the genome you want to parallelize the variant calling and joint-calling steps by (i.e. 50 BED files of 50 intervals to parallelize by)
- A tab-delimited ``sex_linker.tsv`` file with the sample names in one column and sex in the other to identify discordances

Outputs
-------

Depending on which run-mode you have set, you will be able to generate:
- A hard-filtered, multi-sample joint-called VCF in ``full`` and ``joint_genotyping`` mode
- Per-sample gVCFs for all regions of the genome for future joint-calling in ``full`` mode
- Deduplicated and post-BQSR BAM files in ``full`` mode
- Various QC metrics (i.e. fastQC, multiQC, bcftools stats) in all three modes


See the :doc:`installation` and :doc:`usage` for details on setting up and running the pipeline.
