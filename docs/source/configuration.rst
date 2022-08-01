Configuration
==============================

The workflow needs to be configured to perform the analysis of your choice by editing the following files in the ``config/`` folder.  Each file is described in more detail below.

- :ref:`Configuration file`
- :ref:`Manifest file`
- :ref:`Intervals file`
- :ref:`Sex linker file`
- :ref:`MultiQC configuration`

.. _Configuration file:

Configuration file
------------------

The pipeline offers three run modes. Please specify the run mode in ``config.yaml``.  The name of the file defaults to ``config.yaml`` but you can use other filenames in conjunction with Snakemake's ``--configfile`` command line flag.

- **full**: This mode starts with FASTQs and emits a joint-called, filtered, multi-sample VCF.
- **joint_genotyping**: This mode starts with gVCFs and runs joint-calling and filtering, emitting a multi-sample VCF. In the event you have analyzed batches of samples in the full-run mode, these batches can then jointly re-genotyped with this run mode.
- **fastqc_only**: This mode starts with FASTQs and emits trimmed FASTQs as well as a multiQC report for raw and trimmed reads. This run mode is meant for performing QC on FASTQ data before further downstream analysis.

.. _Manifest file:

Manifest file
-------------

You will need to provide a headerless, white-space delimited manifest file to run the pipeline for all three run-modes.  The default name for the file is ``manifest.txt`` but this is user configurable in the config file under ``sampleFile``.

For **full** and **fastqc_only** mode, the ``manifest.txt`` requires the following columns:

- Columns: ``readgroup  sample_ID   path/to/r1.fastq    path/to/r2.fastq``
- ``readgroup`` values should be unique, e.g. sampleID_flowcellID
- ``sample_ID`` should be the same for all FASTQ pairs from a single sample, and can be different from the FASTQ filenames

For example::

    Sample1_S1_L001 Sample1 input/Sample_001_S1_L001_R1.fastq   input/Sample_001_S1_L001_R2.fastq
    Sample1_S1_L002 Sample1 input/Sample_001_S1_L002_R1.fastq   input/Sample_001_S1_L002_R2.fastq

For **joint_genotyping** mode:

- Columns: ``sample_ID   path/to/file.g.vcf.gz``
- ``sample_ID`` values should be unique, and should correspond to the sample IDs in the gVCFs
- gVCFs should be bgzipped and indexed

For example::

    Sample1 vcfs/Sample1.g.vcf.gz
    Sample2 vcfs/Sample2.g.vcf.gz

.. _Intervals file:

Intervals file
--------------

For **full** and **joint_genotyping** modes only.

Joint-calling for a large number of samples is computationally expensive and time-consuming. This pipeline was designed to mitigate these issues by parallelizing joint-calling over multiple intervals of the genome.  To specify the number of intervals, and which regions to parallelize over, a 2-column tab-delmited ``intervals.tsv`` file can be specified.  The filename can be customized and edited in the config file under ``intervalsFile``.

This file contains two columns with headers:

- ``interval_name`` for the name of the particular interval or region
- ``file_path`` full path to the interval/region BED file, Picard-style ``.interval_list``, VCF file, or GATK-style ``.list`` or ``.intervals`` file (see further details on these formats `here <https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists>`_)

For example::

    interval_name   file_path
    interval_1  resources/scattered_calling_intervals/interval_1.bed
    interval_2  resources/scattered_calling_intervals/interval_2.bed


The pipeline will supply these interval files to the GATK ``HaplotypeCaller``, ``GenomicsDBImport``, and ``GenotypeGVCFs`` steps to run concurrent instances of these rules at each specified interval(s), reducing overall execution time.

We recommend specifying regions of equal size for parallelization.

.. _Sex linker file:

Sex linker file
---------------

The pipeline provides a boolean option ``somalier`` to estimate relatedness amongst the samples using `Somalier <https://github.com/brentp/somalier>`_ in the ``config.yaml`` (see ``check_relatedness`` parameter in :doc:`configuration`).  This requires a 2-column, tab-delimited file.  The filename defaults to ``sex_linker.tsv`` and is specified in the ``config.yaml`` under ``sexLinker``.  This file requires:

- First column with the header ``Sample`` with all sample names
- Second column with the header ``Sex`` containing case-insensitive encodings of sex in either m/f or male/female format

For example::

    Sample  Sex
    NA12878 F
    Subject1    female
    Subject2    m

.. _MultiQC configuration:

MultiQC configuration
---------------------

A configuration file for MultiQC can be found in ``config/multiqc.yaml`` and is used for generating and specifying the order of the various modules in the MultiQC report from the pipeline. We **do not** recommend modifying this file unless you understand how this configuration file is setup or how MultiQC works.

.. _Config parameters:

Config parameters
-----------------

Below are descriptions and usage options for the various config parameters specified in ``config.yaml``.

+-------------------------+-----------+------------------------------------------------------------+
| Parameter               | Required  | Description                                                |
+=========================+===========+============================================================+
| ``sampleFile``          |     Y     | Manifest file with IDs                                     |
+-------------------------+-----------+------------------------------------------------------------+
| ``intervalsFile``       |     Y     | File with interval names and file paths                    |
+-------------------------+-----------+------------------------------------------------------------+
| ``jobs``                |     Y     | Max jobs to run concurrently                               |
+-------------------------+-----------+------------------------------------------------------------+
| ``sexLinker``           |     Y     | File with reported sex of each sample ID                   |
+-------------------------+-----------+------------------------------------------------------------+
| ``tempDir``             |     Y     | Location of temp directory; does not have to exist prior   |
|                         |           | to pipeline execution                                      |
+-------------------------+-----------+------------------------------------------------------------+
| ``runType``             |     Y     | Specify run mode to use (see below)                        |
+-------------------------+-----------+------------------------------------------------------------+
|    ``full``             |     Y     | ``[yes|no]`` Set to yes for full run mode                  |
+-------------------------+-----------+------------------------------------------------------------+
|    ``joint_genotyping`` |     Y     | ``[yes|no]`` Set to yes for joint calling from gVCFs       |
+-------------------------+-----------+------------------------------------------------------------+
|    ``fastq_qc_only``    |     Y     | ``[yes|no]`` Set to yes for FASTQ QC and trimming          |
+-------------------------+-----------+------------------------------------------------------------+
| ``global_vars``         |     N     | Set global java options                                    |
+-------------------------+-----------+------------------------------------------------------------+
| ``cluster_mode``        |     N     | Used to submit jobs to a cluster only if you are using     |
|                         |           | the optional wrapper script.  See :doc:`execution`         |
+-------------------------+-----------+------------------------------------------------------------+
| ``default_queue``       |     Y     | Name of your default cluster partition/queue; can be ``~`` |
+-------------------------+-----------+------------------------------------------------------------+
| ``compute_queue``       |     Y     | Name of queue/partition best suited for compute-           |
|                         |           | intensive jobs; can be ``~``                               |
+-------------------------+-----------+------------------------------------------------------------+
| ``memory_queue``        |     Y     | Name of queue/partition best suited for memory-intensive   |
|                         |           | jobs; can be ``~``                                         |
+-------------------------+-----------+------------------------------------------------------------+
| ``center_id``           |     Y     | Name of sequencing center for use in ``@RG`` tag in bams   |
+-------------------------+-----------+------------------------------------------------------------+
| ``max_concurrent``      |     Y     | Max concurrent jobs for specific high-bandwidth rules,     |
|                         |           | to avoid potentially hitting bandwidth caps if deployed    |
|                         |           | in a cloud environment; see wrapper script for an example  |
|                         |           | of how to pass this in to snakemake.  Set to the same      |
|                         |           | number as ``jobs`` if you don't want to limit concurrent   |
|                         |           | rules in this way                                          |
+-------------------------+-----------+------------------------------------------------------------+
| ``max_het_ratio``       |     Y     | Max het/hom ratio to allow through post-calling QC         |
+-------------------------+-----------+------------------------------------------------------------+
| ``min_avg_depth``       |     Y     | Minimum depth required for sample to pass post-calling QC  |
+-------------------------+-----------+------------------------------------------------------------+
| ``max_contam``          |     Y     | Max % contamination to allow through post-calling QC       |
+-------------------------+-----------+------------------------------------------------------------+
| ``time_threshold``      |     Y     | (minutes) Exclude rules from the benchmarking report if    |
|                         |           | elapsed time is below this threshold                       |
+-------------------------+-----------+------------------------------------------------------------+
| ``somalier``            |     Y     | Check relatedness and sex discordance with Somalier        |
|                         |           | (requires sex_linker.tsv) only available in full run       |
|                         |           | mode.  Support of Mac OSX is experimental, so you may      |
|                         |           | want to set this to False on a Mac                         |
+-------------------------+-----------+------------------------------------------------------------+

The remainder of the ``config.yaml`` file contains a selected set of exposed per-tool parameters.  For the most part, this allows tuning of resource allocation on a per-tool basis (i.e. ``threads`` and ``memory`` in MB).  Java-based tools also allow for arbitrary java options to be passed through via ``java_opts``.  Additional exposed parameters include:

- ``genomicsDBImport`` and ``genotypeGVCFs``: We have exposed some useful parameters that have been helpful to adjust as scale increases.  Please see GATK documentation for the relevant tools to learn more.
- ``verifyBamID``: A ``region`` field allows the user to specify chromosomes over which to run contamination analysis, in an attempt to mitigate large memory requirements.
