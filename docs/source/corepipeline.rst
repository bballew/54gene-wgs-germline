Core pipeline
=============

About
-----

This page describe details of the various run-modes available in this pipeline, the rules used within them and further specifications on the tools used. This page provides information on the parameters used for certain tools and behaviours between the run-modes. 

Workflow
--------

Pulling in resources
^^^^^^^^^^^^^^^^^^^^

The pipeline first pulls in reference fasta, annotated VCF and index files (GRCh38 build) from the Broad Institute's public AWS S3 bucket in the ``rule get_resources``. The source code for this rule can be modified to pull from their GCP bucket instead. 

FastQC and Read Trimming
^^^^^^^^^^^^^^^^^^^^^^^^
**In full and fastqc_only run modes:**

The pipeline will create symbolic links for the supplied input fastqs from the ``manifest.txt`` file in ``rule symlink_fastqs`` using a module from ``utilities.py`` to obtain the readgroup portion of the sample names for ease of file naming of FastQC outputs which cannot be renamed and of outputs later on as well. Please bare in mind these symlinks when managing data on your infrastructure. 


``rule fastqc`` will then generate fastQC reports for all input fastqs followed by read trimming and adapter removal (Illumina's TruSeq adapters) using ``fastp``. If you choose to use alternate adapters, you will need to modify the source code for this particular rule. There also exist many additional parameters and functionality for fastp which we have not exposed in the config or included here.

Why fastp over other standard tools such as Trimmomatic? 
We chose fastp over Trimmomatic after encountering run errors and periodic stalling of the tool in our initial development of the pipeline. 


Post-trimming, fastQC will be run again on the trimmed reads. This results in fastQC results for the raw input reads in ``results/fastqc``, and post-trimmed reads in ``results/post_trimming_fastqc``.

Read Alignment
^^^^^^^^^^^^^^

**In full mode only:**

Trimmed reads are aligned to the reference genome using BWA in ``rule align_reads``. The read1 and read2 data are combined into a single output BAM per fastq pair. If samples were run over several lanes (i.e. 4 lanes), the read1 and read2 fastq pair will be aligned separately then combined during a deduplication step later on (see ``rule mark_duplicates``). This helps speed up the alignment process and helps parallelize this step for efficiency. The output of this rule will be a sorted BAM for each sample within ``results/mapped/{rg}.bam``.

  
Duplicates for each sorted sample BAM are the marked and removed. Multi-lane sample BAMs will be concatenated in this step. The pipeline will then use GATK's BaseRecalibrator to generate a recalibration table for base quality score using known sites VCF pulled from the Broad's resources bucket. This recalibration is then applied to each BAM in ``rule apply_bqsr``. Samtools stats are then generated for all BAMs. See the :doc:`investigateresults` page for more information on where to find the stats. 

We are currently using the dbSNP VCFs of known SNPs and the known indels VCF from the GATK resource bundle, available on `AWS S3 <s3://broad-references/hg38/v0/>`_ and `GCP <https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0>`_. 

  
Variant Calling
^^^^^^^^^^^^^^^

Single-sample gVCFs are generated using GATK's HaplotypeCaller and parallelized over a user-specified number of intervals and/or regions of the genome using the ``intervals.tsv`` file. A gVCF for each sample will be generated for each specified interval/region. This method allows for parallelization and reduction in overall execution time for variant calling. Following the `GVCF workflow <https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller>`_, these are to be used for joint genotyping of multiple samples later in the pipeline for scalable analyses.

Each gVCF will be compressed and then concatenated across each interval/region to generate a single gVCF for each sample along with its index, with the naming convention of ``/HaplotypeCaller/called/{sample}_all_regions.g.vcf.gz``. 

Joint genotyping
^^^^^^^^^^^^^^^^

**In joint_genotyping mode only:**

It is at this step in the workflow that a secondary entry point is provided when the run-mode in the ``config.yaml`` is set to ``joint_genotyping``. In this run mode, the gVCFs provided in the ``manifest.txt`` file and their indexes will be symlinked to a subdirectory within ``/results/HaplotypeCaller/called/`` prior to continuing on to the rest of the workflow.
  
**In joint_genotyping and full run modes:**

In order to perform joint-genotyping over multiple samples using GATK's GenotypeGVCFs, the input gVCFs must be consolidated across samples as the genotyping step can only take one single input. To circumvent this issue, we use GATK's GenomicsDBImport in ``rule HC_consolidate_gvcfs`` to generate database stores for each sample, parallelized again across intervals/regions, to then pass into GenotypeGVCFs. DBImport can potentially take up alot of /temp space so it is recommended that ``--tmp-dir`` be used to redirect to a larger temp space. The ``--batch-size`` and ``--reader-threads`` parameters can be tweaked in the ``config.yaml`` to read in more data stores concurrently or in larger batch sizes but the default settings are those suggested by GATK developers. 

Joint genotyping using the various database stores created is performed in ``rule HC_genotype_gvcfs`` to emit a genotyped gVCF for each interval/region in ``results/HaplotypeCaller/genotyped/{interval}.vcf.gz``. The ``--max_alt_alleles`` to genotype and ``--max_genotype_count`` for each site can be tweaked in the ``config.yaml``. 

We exposed these and other parameters for GenomicsDBImport after encountering `recent issues <https://github.com/broadinstitute/gatk/issues/7542>`_ where the ``--max-alternate-alleles`` flag for GenotypeGVCFs was set at a default of 6 but was not actually being applied as a threshold. A fix in GATK v4.2.4.1 attempted to apply this threshold instead resulted in a bug where the relevant locus would result in the program crashing. Subsequently, an `update in v4.2.5.0 <https://github.com/broadinstitute/gatk/pull/7655>`_ introduced a new parameter for GenotypeGVCFs called ``--genomicsdb-max-alternate-alleles`` and required to be minimum one greater than ``--max-alternate-alleles`` to account for the NON_REF allele.   


The per-interval/region, genotyped gVCFs will be concatenated into one sorted, indexed, project-level multi-sample gVCF for downstream analysis in ``results/HaplotypeCaller/genotyped/HC_variants.vcf.gz``.  
  
*Note*: While GenomicsDBImport supports adding N+1 samples to the datastores, however, our pipeline does not utilize this functionality and instead creates the databases every time from scratch. This was a development choice made to avoid issues with potential failures with maintaining the datastores and revisitng them in future analyses. 

.. _variant_filtration:

Variant Filtration
^^^^^^^^^^^^^^^^^^

The project-level VCF is normalized and split into a separate VCF with multiallelics to allow for further hard-filtering into SNPs and indels later. Hard-filtering using GATK's VariantFiltration tool is performed separately on the SNP and indel-specific project-level VCFs in ``rule hard_filter_snps`` and ``rule_hard_filter_indels``. For more information on how we perform hard-filtering, see GATK's documentation on hard-filtering recommendations `here <https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants>`_.

*Note*: We currently do not remove the filtered sites themselves from the VCF but instead just update the filter field. 

Post-Calling QC
^^^^^^^^^^^^^^^

Contamination Check
*******************

**In full mode only:**

As an added QC measure, we perform a contamination check on the BAM files using a tool called `VerifyBamID <https://genome.sph.umich.edu/wiki/VerifyBamID>`_. This tool will verify whether reads in a particular BAM are concordant with previous genotypes found for a specific sample. This can help identify potential sample mixture. 

The tool normally takes the entire BAM file as an input but to reduce the computational burden of performing this check, we opted to only subset particular chromosomes (ideally one or two) from the BAM files to perform the check. We have found that is this sufficient for initial flagging of contamination for further ni-depth investigation of troublesome samples. We allow the ability to select these chromosomes within the ``config.yaml``. 

This step in ``rule contamination_check`` will output various contamination metrics for each sample BAM file that are combined in a summary file. This summary file will be later used for automated filtering of samples out of the project-level VCF based on thresholds defined in the ``config.yaml``. See :ref:`exclude_samples` section for more information.

Merging calls
*************

The hard-filtered SNP and indels generated in :ref:`variant_filtration` are merged into one sorted, and indexed project-level hard-filtered multi-sample VCF in ``results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz``. We subsequnently run GATK's CollectVariantCallMetrics on this VCF file. See :doc:`investigateresults` for more info. 

Checking Relatedness with Somalier
**********************************

If specified as 'Yes' for ``check_relatedness`` in the ``config.yaml``, the pipeline will run Somalier to check for relatedness amongst the samples. `Somalier <https://github.com/brentp/somalier>`_ is a tool that can be used to check any number of samples from joint-called VCFs for identity and to infer relationships. The tool ideally requires a jointly-called cohort VCF and PED file of expected sexes and relationships. An example of the Somalier output can be found `here <https://brentp.github.io/somalier/ex.html>`_. 

In our pipeline, we create a PED-style file from the ``sex_linker.tsv`` specified in the ``config.yaml`` using a module that will convert the linker file to the appropriate format for Somalier in ``utilities/create_ped.py``. 

How do we use this tool?

We use this tool as a rough estimate of relatedness and to highlight potential concerns within our data for putative relatedness or with unexpected duplicates. If identified, we then perform more thorough analyses using other tools such as KING, plink, graf and etc. 

This is largely because Somalier uses the following equation to determine relatedness::

    (shared-hets)(i,j)-2*ibs0(i,j)/min (hets(i),hets(j))


This assumes, as noted in their `publication <https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00761-2>`_, that the sites they've selected on which to assess relatedness are "high-quality, unlinked sites with a population allele frequency of around 0.5."  We suspect this will not hold true across all populations. It is unclear how much this will degrade across multiple populations with some degree of shared ancestry.  Note that the relatedness value will always be depressed when comparing, for example, NA12878 with Nigerian subjects.

*Note*: We currently run this tool using our organization's conda channel. You will likely have to modify the ``envs/somalier.yaml`` to point elsewhere to a different package location.

Sex Check
*********

In addition to Somalier above, we also use bcftools' guess-ploidy to determine sample sex from genotype likelihoods. These results are included in the MultiQC report generated at the end of the post-calling QC stage. See :ref:`multiqc` for more information. 

.. _exclude_samples:

Sample Exclusions
*****************

We exclude samples from the project-level hard-filtered VCF using a module found in ``script/create_exclude_list.py`` in ``rule create_exclude_list`` that will use the metrics and information generated from the contamination check and bcftools stats to exclude samples based on the following default thresholds:

- Max het/hom ratio of 2.5 
- Minimum average depth of 20
- Maximum contamination estimate of 0.03 (only used if run in full run mode)

These thresholds can be tweaked in the ``config.yaml``. A list of samples to exclude and another list with these samples and annotations for why they were excluded will be generated in ``results/post_qc_exclusions/``.

Post sample exclusion, another sorted and indexed, project-level, hard-filtered VCF will emitted in ``results/post_qc_exclusions/samples_excluded.HC_variants.hardfiltered.vcf.gz``. This is to avoid altering the original project-level, hard-filtered VCF. Note that the ID column here will also be updated to ``%CHROM:%POS:%REF:%ALT`` using bcftools annotate. If you want this field to updated to something else, you may choose to modify the source code in ``rule exclude_samples``. 

.. _multiqc:

MultiQC
*******

A multiQC report is generated for all three run-modes and will differ in content depending on which post-calling QC checks were performed. 

For **fastqc_only** run mode, the multiQC report will include:

- Pre- and post-read-trimming fastQC results 

For the **full** run mode, the multiQC report will include:

-  Pre- and post-read-trimming fastQC results
-  Bcftool stats on joint-called variants 
-  Deduplication metrics for BAM files
-  Sex check results from bcftools guess-ploidy 
-  Contamination check results from verifyBamID
-  If specified in config, relatedness check results from Somalier
-  Variant calling metrics

For **joint_genotyping** mode, the multiQC report will include:

- Variant calling metrics
- Sex check results from bcftools guess-ploidy
- Bcftool stats on joint-called variants 
- If specified in config, relatedness check results from Somalier 