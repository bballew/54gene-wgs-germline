Installation and Setup
======================

This workflow was designed to use `conda`` for dependency management and utilizes `Snakemake` as the workflow management system, ensuring reproducible and scalable analyses.

This installation guide is designed for Unix/Linux environments.

Step 1. Obtain a copy of this workflow
--------------------------------------

If you haven't already, clone this repository to your local system, into the place where you want to perform the data analysis::

    git clone git@gitlab.com:data-analysis5/54gene-wgs-germline.git

Step 2. Install the run-time conda environment
----------------------------------------------
If you don't already have conda installed, follow this guide `here <https://docs.conda.io/en/latest/miniconda.html#installing>`_ to install Miniconda.

Once installed, create the run-time conda environment with minimal dependencies defined using the following command::
    
    conda env create -f environment.yaml

Step 3. Configure the workflow
------------------------------
The workflow needs to configured to perform the analysis of your choice by editing the following files in the `config/` folder:

A. Config file
^^^^^^^^^^^^^^
As of v1.0, the pipeline offers three run modes. Please specify the run mode in `config.yaml`.

- **full**: This mode starts with fastqs and emits a joint-called, filtered, multi-sample VCF
- **joint_genotyping**: This mode starts with gVCFs and runs joint-calling and filtering, emitting a multi-sample VCF. In the event you have analyzed batches of samples in the full-run mode, these batches can then jointly re-genotyped with this run mode.
- **fastqc_only**: This mode starts with fastqs and emits trimmed fastqs as well as a multiQC report for raw and trimmed reads. This run mode is meant for performing QC on fastq data before further downstream analysis.

A. Manifest file
^^^^^^^^^^^^^^^^
You will need to provide a headerless, white-space delimited manifest file to run the pipeline for all three run-modes. 

For **full** and **fastqc_only** mode:

The `manifest.txt` would include the following columns:

+------------+-----------+----------------+-----------------+
| readgroup  | sample_ID |path/to/r1.fastq| path/to/r2.fastq|
+------------+-----------+----------------+-----------------+

Where `readgroup` values are unique (i.e. sampleID_barcode) and `sample_ID` values are the same for all fastq pairs from a single sample and can be different from the fastq filenames themselves.

For example:

+--------------------+-----------+-----------------------------------+---+
| Sample001_S1_L001  | Sample001 | fastqs/Sample_001_S1_L004_R1.fastq|...|
+--------------------+-----------+-----------------+-----------------+---+

For **joint_genotyping** mode:

The `manifest.txt` would include the following two columns:

+-------------+-----------------------------+
| sample_ID   |  path/to/sample_ID.g.vcf.gz |
+-------------+-----------------------------+

For example:

+---------------+-----------------------------+
| Sample_001    |  vcfs/Sample_001.g.vcf.gz   |
+---------------+-----------------------------+

*Note*: The gVCFs should be zipped and indexed. 