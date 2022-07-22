Usage and Execution
==============================

Configure the workflow
------------------------------

The workflow needs to configured to perform the analysis of your choice by editing the following files in the ``config/`` folder:

Configuration file
^^^^^^^^^^^^^^^^^^

As of v1.0, the pipeline offers three run modes. Please specify the run mode in ``config.yaml``.

- **full**: This mode starts with fastqs and emits a joint-called, filtered, multi-sample VCF.
- **joint_genotyping**: This mode starts with gVCFs and runs joint-calling and filtering, emitting a multi-sample VCF. In the event you have analyzed batches of samples in the full-run mode, these batches can then jointly re-genotyped with this run mode.
- **fastqc_only**: This mode starts with fastqs and emits trimmed fastqs as well as a multiQC report for raw and trimmed reads. This run mode is meant for performing QC on fastq data before further downstream analysis.

B. Manifest file
^^^^^^^^^^^^^^^^
You will need to provide a headerless, white-space delimited manifest file to run the pipeline for all three run-modes. 

For **full** and **fastqc_only** mode, the ``manifest.txt`` requires the following columns:

- First column with the readgroup for each sample
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
Config Parameters
-----------------

Below are descriptions and usage options for the various config parameters specified in ``config.yaml``.

+---------------+-----------+----------------+---------------------------+
| Parameter     |  Required |  Default value |         Description       |
+===============+===========+================+===========================+
| sampleFile    |     Y     |      NA        |  Manifest file with IDs   |
+---------------+-----------+----------------+---------------------------+
| intervalsFile |     Y     |      NA        | File with interval names  |
|               |           |                | and file paths            |
+---------------+-----------+----------------+---------------------------+
| jobs          |     Y     |     1000       | Jobs to submit to cluster |
+---------------+-----------+----------------+---------------------------+
| sexLinker     |     Y     |      NA        | File with reported sex of |
|               |           |                | sexes of each sample ID   |
+---------------+-----------+----------------+---------------------------+
| tempDir       |     Y     |     temp/      | Location of temp directory|
+---------------+-----------+----------------+---------------------------+
| runType       |     Y     |    full: yes   | Specify run mode to use   |
+---------------+-----------+----------------+---------------------------+
| global_vars   |     N     |  As specified  | Set global java options   |
+---------------+-----------+----------------+---------------------------+
| cluster_mode  |     N     |  As specified  | Used to submit jobs to a  |
|               |           |                | cluster. See              |
|               |           |                | :ref:`execution`          |
+---------------+-----------+----------------+---------------------------+
| default_queue |     Y     |     "big"      | Name of your default      |
|               |           |                | cluster partition/queue   |
+---------------+-----------+----------------+---------------------------+
| compute_queue |     Y     |     "big"      | Name of queue/partition   |
|               |           |                | with compute-heavy nodes  |
+---------------+-----------+----------------+---------------------------+
| memory_queue  |     Y     |     "big"      | Name of queue/partition   |
|               |           |                | with nodes for memory-    |
|               |           |                | intensive jobs            |
+---------------+-----------+----------------+---------------------------+
| bed           |     N     | `Homo_sapiens_ | Split genome by chromosome|
|               |           | assembly38.bed`| when calling variants     |
+---------------+-----------+----------------+---------------------------+
| max_concurent |     Y     |      40        |Max concurrent running jobs|
|               |           |                |to run at once. Allows you |
|               |           |                |to throttle heavy rules    |
|               |           |                |depending on your env      |
+---------------+-----------+----------------+---------------------------+
| max_het_ratio |     Y     |      2.5       | Max het/hom ratio to allow|
|               |           |                |                           |
+---------------+-----------+----------------+---------------------------+
| min_avg_depth |     Y     |      20        | Minimum depth required for|
|               |           |                | sample                    |
+---------------+-----------+----------------+---------------------------+
| max_contam    |     Y     |      0.03      | Max % of contam allowed   |
+---------------+-----------+----------------+---------------------------+
|time_threshold |     Y     |    5 (mins)    | Exclude rules from the    |
|               |           |                | benchmarking report if    |
|               |           |                | elapsed time is below this|
+---------------+-----------+----------------+---------------------------+
| somalier      |     Y     |      True      | Check relatedness and sex |
|               |           |                | discordance with Somalier |
|               |           |                | (requires sex_linker.tsv) |
|               |           |                | only available in full    |
|               |           |                | run mode                  |
+---------------+-----------+----------------+---------------------------+

Resource Allocation
-------------------

This pipeline was originally developed to be run on an Unix-based HPC system for scalable and efficient analyses. As a result, there are several configuration options available for allocating resources and running jobs within the workflow.


Within the ``config/config.yaml`` there are several options for allocating memory, threads, and setting Java-specific variables (for GATK suite of tools) for the various tools used in the workflow.

- To modify the number of maximum jobs to run, you can specify a number for ``jobs``
- Set a maximum number of jobs to run concurrently if you have bandwidth constraints using the ``max_concurrent`` variable
- Specify which paritions/queues to submit your jobs to, depending on your HPC using the ``default_queue``, ``compute_queue`` and ``memory_queue`` variables (some jobs such as alignment, require more memory and can use nodes with more memory available when specified with these variables)
- Specify the number of threads and memory in MB for each tool, where available using the ``threads`` and ``memory`` variables
- Specify the space to allocate for Java class metadata using the ``global_vars`` variable

