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

B. Manifest file
^^^^^^^^^^^^^^^^
You will need to provide a headerless, white-space delimited manifest file to run the pipeline for all three run-modes. 

For **full** and **fastqc_only** mode:

The `manifest.txt` would include the following columns:

- First column with the readgroup for each sample (contains the full sample ID, barcode, and lane)
- Second column with only the sample ID 
- Third column with the path to the read 1 FASTQ file
- Fourth column with the path to the read 2 FASTQ file 
  
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

- First column with the sample IDs for each gVCF
- Second column with the paths to the gVCFs

+-------------+-----------------------------+
| sample_ID   |  path/to/sample_ID.g.vcf.gz |
+-------------+-----------------------------+

For example:

+---------------+-----------------------------+
| Sample_001    |  vcfs/Sample_001.g.vcf.gz   |
+---------------+-----------------------------+

*Note*: The gVCFs should be zipped and indexed. 

C. Intervals file
^^^^^^^^^^^^^^^^^

For **full** and **joint_genotyping** modes only.

Joint-calling for a large number of samples can be computationally expensive and very time-consuming. This pipeline was designed to mitigate these issues by parallelizing joint-calling over multiple intervals of the genome. To specify the number of intervals, and which regions to parallelize over, a 2-column tab-delmited ``intervals.tsv`` file can be specified. 

This file contains two columns:

- ``interval_name`` for the name of the particular interval or region 
- ``file_path`` full path to the interval/region BED file, Picard-style ``.interval_list``, VCF file, or GATK-style ``.list`` or ``.intervals`` file (see further details on these formats `here <https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists>`_)

For example:

+---------------+-------------------------------------------------------+
| interval_name |  file_path                                            |
+---------------+-------------------------------------------------------+
| interval_1    | /resources/scattered_calling_intervals/interval_1.bed |
+---------------+-------------------------------------------------------+


The pipeline will supply these interval files to the GATK ``HaplotypeCaller``, ``GenomicsDBImport`` and ``GenotypeGVCFs`` steps to run concurrent instances of these rules at each specified interval(s), reducing overall execution time.

We recommend specifying regions of equal size for parallelization.

D. Sex linker file
^^^^^^^^^^^^^^^^^^

The pipeline provides an option to check the relatdness amongst the samples using Somalier in the ``config.yaml`` (see ``check_relatedness`` parameter in :doc:`usage`). This requires a 2-column, tab-delimited ``sex_linker.tsv`` file provided and specified in the ``config.yaml``. This file will have:

- First column with the header ``Sample`` with all sample names 
- Second column with the header ``Sex`` containing one-letter formattted sex of all samples 

For example:

+---------+-----+
| Sample  | Sex |
+---------+-----+
| NA12878 | F   |
+---------+-----+


E. MultiQC yaml
^^^^^^^^^^^^^^^

A configuration file for MultiQC can be found in ``config/multiqc.yaml`` and is used for generating and specifying the order of the various modules in the multiQC report from the pipeline. We **do not** recommend modifying this file unless you understand how this configuration file is setup or how multiQC works. 