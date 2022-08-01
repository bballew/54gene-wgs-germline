Core pipeline
=============

This page describe details of the various run-modes available in this pipeline, the rules used within them and further specifications on the tools used. This page provides information on the parameters used for certain tools and behaviours between the run-modes.

Pulling in resources
--------------------

There is no config option to point to a reference genome.  The pipeline automatically pulls in the required GRCh38 reference files, including the reference genome and all requisite indices as well as the known SNP and indel files for running BQSR, from the Broad Institute's public AWS S3 bucket (``s3://broad-references/hg38/v0/``) in the ``rule get_resources``.  We have not provided an option for hg19/GRCh37.

FastQC and read trimming
------------------------

**In full and fastqc_only run modes:**

The pipeline will create symbolic links for the supplied input FASTQs from the manifest file in ``rule symlink_fastqs``.  FastQC generates reports based on filename, and has no option to change the output filenames.  Symlinking allows harmonization of filenames to the convention followed throughout the pipeline; for FASTQs, that convention is ``<readgroup>_r[12].fastq.gz``.  Please bear in mind these symlinks when managing data on your infrastructure.


``rule fastqc`` will then generate FastQC reports for all input FASTQs.  Note that this is one of the rules governed by the ``max_concurrent`` config argument (see :ref:`Config parameters`).  On filesystems where IO bandwidth is capped, you may want to control the number of concurrent rules running at this stage.

Next, we perform read trimming and adapter removal (currently hard-coded to use Illumina's TruSeq adapters) using `fastp <https://github.com/OpenGene/fastp>`_.  If you need to use alternate adapters or adjust other ``fastp`` parameters, please submit a feature request to expose these as parameters in config space.

Post-trimming, FastQC will be run again on the trimmed reads.  This results in FastQC results for the raw input reads in ``results/fastqc``, and post-trimmed reads in ``results/post_trimming_fastqc``.  Review these reports to monitor read quality and effective adapter removal.

Note that **fastqc_only** run mode will stop here, allowing a quick turnaround in sharing read quality information with lab, and assessing whether there are any samples to drop before performing a more computationally costly full run.

Read alignment, deduplication, and BQSR
---------------------------------------

**In full mode only:**

Trimmed reads are aligned to the reference genome using `BWA <https://github.com/lh3/bwa>`_ in ``rule align_reads``. The read1 and read2 data are combined into a single output BAM per FASTQ pair. If samples were run over several lanes (e.g. 4 lanes), each per-lane read1 and read2 FASTQ pair will be aligned individually, then combined during the subsequent deduplication step (see ``rule mark_duplicates``). This helps with efficient alignment by running multiple smaller alignments in parallel. The read group IDs of the BAM files will include the sequencing center ID specified in the ``config.yaml`` under the ``center_id`` parameter.

The alignment step outputs aligned and sorted BAMs, one per sample-and-lane, at ``results/mapped/<readgroup>.bam``. These BAMs are flagged as temp, so they are automatically removed unless run with the ``--notemp`` Snakemake flag (see `Snakemake documentation <https://snakemake.readthedocs.io/en/stable/executing/cli.html#BEHAVIOR>`_).

After alignment, duplicate reads in the sorted BAMs generated for each readgroup are then marked and removed in ``rule mark_duplicates``.  It is at this step that samples split over multiple lanes will be merged, and subsequently named with the sample ID provided in the manifest. This generates one BAM file for each sample, found as ``results/dedup/<sample ID>.bam``. Subsequently, the pipeline will use GATK's BaseRecalibrator to generate a recalibration table for base quality scores using known sites VCF pulled from the Broad's resources bucket. This recalibration is then applied to each BAM in ``rule apply_bqsr``.  Samtools stats are then generated for all BAMs.  See the :doc:`investigateresults` page for more information on where to find the stats. Upon recalibration, the per-sample, sorted BAM files and their indexes can be found in ``results/bqsr/<sample ID>.bam``.

Variant calling
---------------

Per-sample gVCFs are generated in ``rule HC_call_variants`` using GATK's HaplotypeCaller.  Calling is parallelized over a user-specified number of intervals and/or regions of the genome using the interval file listed in the config.  A temp-flagged gVCF for each sample will be generated for each specified interval/region; these are automatically cleaned up once they are all present and have been successfully combined into a single per-sample gVCF using GATK's GatherVcfs.  This method allows for parallelization and reduction in overall execution time for variant calling.  Following the `GVCF workflow <https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller>`_, these are to be used for joint genotyping of multiple samples later in the pipeline for scalable analyses.  The resulting per-sample gVCF is compressed and indexed, and can be found at ``results/HaplotypeCaller/called/<sample>_all_regions.g.vcf.gz``.

Joint genotyping
----------------

**In joint_genotyping mode only:**

It is at this step in the workflow that a second entry point is provided when the run mode in the ``config.yaml`` is set to ``joint_genotyping``.  In this run mode, the gVCFs provided in the manifest file and their indices will be symlinked to a subdirectory within ``/results/HaplotypeCaller/called/`` prior to continuing on to the rest of the workflow.

**In joint_genotyping and full run modes:**

In order to perform joint-genotyping over multiple samples using GATK's GenotypeGVCFs, the input gVCFs must be consolidated across samples as the genotyping step can only take one single input. To circumvent this issue, we use GATK's GenomicsDBImport in ``rule HC_consolidate_gvcfs`` to generate database stores for each sample, parallelized again across intervals/regions, to then pass into GenotypeGVCFs.  DBImport can potentially take up a lot of temp space so it is recommended that ``--tmp-dir`` be used to redirect to a larger temp space.  The ``--batch-size`` and ``--reader-threads`` parameters can be tweaked in the ``config.yaml`` to read in more data stores concurrently or in larger batch sizes but the default settings are those suggested by GATK developers.

Joint genotyping using the various database stores created is performed in ``rule HC_genotype_gvcfs`` to emit a genotyped gVCF for each interval/region in ``results/HaplotypeCaller/genotyped/{interval}.vcf.gz``. The ``--max_alt_alleles`` to genotype and ``--max_genotype_count`` for each site can be tweaked in the ``config.yaml``.

We exposed these and other parameters for GenomicsDBImport after encountering `recent issues <https://github.com/broadinstitute/gatk/issues/7542>`_ where the ``--max-alternate-alleles`` flag for GenotypeGVCFs was set at a default of 6 but was not actually being applied as a threshold.  A fix in GATK v4.2.4.1 attempted to apply this threshold, but instead resulted in a bug where the tool would crash upon reaching a locus exceeding this threshold. Subsequently, an `update in v4.2.5.0 <https://github.com/broadinstitute/gatk/pull/7655>`_ introduced a new parameter for GenotypeGVCFs called ``--genomicsdb-max-alternate-alleles``, which is required to be minimum one greater than ``--max-alternate-alleles`` to account for the NON_REF allele.

The per-interval/region, genotyped gVCFs will be concatenated into one sorted, indexed, project-level multi-sample gVCF for downstream analysis in ``results/HaplotypeCaller/genotyped/HC_variants.vcf.gz``.

*Note*: While GenomicsDBImport supports adding N+1 samples to the datastores, our pipeline does not utilize this functionality and instead creates the databases every time from scratch.  This was a development choice made to avoid issues with potential failures with maintaining the datastores and revisiting them in future analyses.

Variant filtering
-----------------

The project-level VCF is normalized and multiallelics are split using ``bcftools norm`` in ``rule split_multiallelics``.  This means that the resulting VCF may have multiple lines representing the same genomic position.  This is conformant with VCF specifications, and may not be expected as input by all downstream tools.  We have elected to split multiallelics for several reasons, including:

- Inability to apply hard filtering to multi-type loci.  GATK's hard filters require first splitting indels and SNPs; multi-type loci don't get split into either category.  So, by splitting multiallelics, you can apply the appropriate filter to all alt alleles
- Difficulty in parsing which annotations refer to which allele after using a tool like VEP or SNPeff

Hard-filtering using GATK's VariantFiltration tool is performed separately on the SNP and indel-specific project-level VCFs in ``rule hard_filter_snps`` and ``rule_hard_filter_indels``.  After variants are flagged in the FILTER column based on hard filters, indels and snps are recombined and can be found at ``results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz``.  For more information on how we perform hard-filtering, see GATK's `documentation <https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants>`_ on hard-filtering recommendations.

*Note*: We currently do not remove the filtered sites themselves from the VCF but instead just update the filter field.  You will want to do a pass with GATK or bcftools to filter out non-PASS variants.

Post-calling QC
---------------

Contamination Check
^^^^^^^^^^^^^^^^^^^

**In full mode only:**

As an added QC measure, we perform a contamination check on the BAM files using a tool called `VerifyBamID <https://genome.sph.umich.edu/wiki/VerifyBamID>`_. This tool estimates the most likely proportion of contaminant DNA present in a sample given phred likelihoods of actual basecalls, assuming HWE.

The tool normally takes the entire BAM file as an input but to reduce the computational burden of performing this check, we opted to only subset particular chromosomes (ideally one or two) from the BAM files to perform the check.  We have found that is this sufficient for initial flagging of contamination for further in-depth investigation of troublesome samples.  We allow the ability to select these chromosomes within the ``config.yaml``.

This step in ``rule contamination_check`` will output various contamination metrics for each sample BAM file that are combined in a summary file.  This summary file will be later used for automated filtering of samples out of the project-level VCF based on thresholds defined in the ``config.yaml``.  See the :ref:`Sample exclusions` section for more information.

Checking relatedness with Somalier
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If ``check_relatedness`` is set to ``yes`` in the ``config.yaml``, the pipeline will run Somalier to check for relatedness amongst the samples. `Somalier <https://github.com/brentp/somalier>`_ is a tool that can be used to check any number of samples from joint-called VCFs for identity and to infer relationships.  The tool takes as input a jointly-called cohort VCF and PED file of expected sexes and relationships.  Our pipeline requires a simple sex linker file described in :doc:`configuration` and creates the PED file for you.  An example of the Somalier output can be found `here <https://brentp.github.io/somalier/ex.html>`_.

This tool provides a rough estimate of relatedness which we mainly use to identify unexpected genetic duplicates.  To confirm specific relationships, we perform a second pass evaluation of the relevant samples using more specialized software, e.g. KING, graf, etc.  Somalier uses the following equation to determine relatedness::

    (shared-hets)(i,j)-2*ibs0(i,j)/min (hets(i),hets(j))

This assumes, as noted in their `publication <https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00761-2>`_, that the sites they've selected on which to assess relatedness are "high-quality, unlinked sites with a population allele frequency of around 0.5."  We suspect this will not hold true across all populations, and we are currently working in a relatively underrepresented ancestry group.  It is unclear how much this will degrade across multiple populations with some degree of shared ancestry.  Note that the relatedness value will always be depressed when comparing samples from disparate ancestries, for example, NA12878 with continental African subjects.

Sex Check
^^^^^^^^^

Somalier also provides functionality to assess sex discordance.  The HTML report provided by Somalier, and in the MultiQC report that ingests this data, includes a plot of scaled mean depth on X vs. self-reported sex.  This plot allows quick identification of disagreement between reported and genetic sex.

In addition to Somalier, we also use bcftools' guess-ploidy plugin to determine sample sex from genotype likelihoods.  These results are also included in the MultiQC report generated at the end of the post-calling QC stage. See :ref:`multiqc` for more information.

.. _Sample exclusions:

Sample exclusions
^^^^^^^^^^^^^^^^^

We exclude samples from the project-level hard-filtered VCF in ``rule create_exclude_list`` based on metrics and information generated from the contamination check and bcftools stats.  Samples are excluded based on the following default thresholds:

- Max het/hom ratio of 2.5
- Minimum average depth of 20
- Maximum contamination estimate of 0.03 (only used if run in full run mode)

These thresholds can be tweaked in the ``config.yaml``.  A list of samples to exclude and another list with these samples and annotations for why they were excluded will be generated in ``results/post_qc_exclusions/``.

Post sample exclusion, another sorted and indexed, project-level, hard-filtered VCF will emitted in ``results/post_qc_exclusions/samples_excluded.HC_variants.hardfiltered.vcf.gz``.  Note that the ID column here will also be updated to ``CHROM:POS:REF:ALT`` using bcftools annotate.

.. _multiqc:

MultiQC
^^^^^^^

A MultiQC report is generated for all three run-modes and will differ in content depending on which post-calling QC checks were performed.

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
