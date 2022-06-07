Core pipeline
=============

About
#####

This page describe details of the various run-modes available in this pipeline, the rules used within them and further specifications on the tools used. This page provides information on the parameters used for certain tools and behaviours between the run-modes. 

Workflow
********

Pulling in resources
^^^^^^^^^^^^^^^^^^^^

The pipeline first pulls in reference fasta, annotated VCF and index files (GRCh38 build) from the Broad Institute's public AWS S3 bucket in the ``rule get_resources``. The source code for this rule can be modified to pull from their GCP bucket instead. 

FastQC and Read Trimming
^^^^^^^^^^^^^^^^^^^^^^^^
**In full and fastqc_only run modes:**
The pipeline will create symbolic links for the supplied input fastqs from the ``manifest.txt`` file in ``rule symlink_fastqs`` using a module from ``utilities.py`` to obtain the readgroup portion of the sample names for ease of naming files thereon. 

``rule fastqc`` will then generate fastQC reports for all input fastqs followed by read trimming and adapter removal (Illumina's TruSeq adapters) using ``fastp``. 

Post-trimming, fastQC will be run again on the trimmed reads. This results in fastQC results for the raw input reads in ``results/fastqc``, and post-trimmed reads in ``results/post_trimming_fastqc``.

Read Alignment
^^^^^^^^^^^^^^

**In full mode only:**

Trimmed reads are aligned to the reference genome using BWA in ``rule align_reads``. The read1 and read2 data are combined into a single output BAM per fastq pair. If samples were run over several lanes (i.e. 4 lanes), the read1 and read2 fastq pair will be aligned separately then combined during a deduplication step later on (see ``rule mark_duplicates``). This helps speed up the alignment process and helps parallelize this step for efficiency.

Duplicates for each sample BAM are the marked and removed. Multi-lane sample BAMs will be concatenated in this step. The pipeline will then use GATK's BaseRecalibrator to generate a recalibration table for base quality score using known sites VCF pulled from the Broad's resources bucket. This recalibration is then applied to each BAM in ``rule apply_bqsr``. Samtools stats are then generated for all BAMs. 

Variant Calling
^^^^^^^^^^^^^^^

Single-sample gVCFs are generated using GATK's HaplotypeCaller and parallelized over a user-specified number of intervals and/or regions of the genome using the ``intervals.tsv`` file. A gVCF for each sample will be generated for each specified interval/region. This method allows for parallelization and reduction in overall execution time for variant calling. Following the GVCF workflow, these are to be used for joint genotyping of multiple samples later in the pipeline for scalable analyses.

Each gVCF will be compressed, indexed and then concatenated across each interval/region to generate a single gVCF for each sample with the naming convention of ``/HaplotypeCaller/called/{sample}_all_regions.g.vcf.gz``. 

Joint genotyping
^^^^^^^^^^^^^^^^

**In joint_genotyping mode only:**
It is at this step in the workflow that a secondary entry point is provided when the run-mode in the ``config.yaml`` is set to ``joint_genotyping``. In this run mode, the gVCFs provided in the ``manifest.txt`` file and their indexes will be symlinked prior to continuing on to the rest of the workflow.

**In joint_genotyping and full run modes:**
In order to perform joint-genotyping over multiple samples using GATK's GenotypeGVCFs, the input gVCFs must be consolidated across samples as the genotyping step can only take one single input. To circumvent this issue, we use GATK's GenomicsDBImport in ``rule HC_consolidate_gvcfs`` to generate database stores for each sample, parallelized again across intervals/regions, to then pass into GenotypeGVCFs. DBImport can potentially take up alot of /temp space so it is recommended that ``--tmp-dir`` be used to redirect to a larger temp space. The ``--batch-size`` and ``--reader-threads`` parameters can be tweaked in the ``config.yaml`` to read in more data stores concurrently or in larger batch sizes but the default settings are those suggested by GATK developers. 

Joint genotyping using the various database stores created is performed in ``rule HC_genotype_gvcfs`` to emit a genotyped gVCF for each interval/region in ``results/HaplotypeCaller/genotyped/{interval}.vcf.gz``. The ``--max_alt_alleles`` to genotype and ``--max_genotype_count`` for each site can be tweaked in the ``config.yaml``. 

The subsequent per-interval/region, genotyped gVCFs are concatenated into one project-level multi-sample gVCF for downstream analysis in ``results/HaplotypeCaller/genotyped/HC_variants.vcf.gz``.

Variant Filtration
^^^^^^^^^^^^^^^^^^

The project-level VCF is normalized and split into a separate VCF with multiallelics to allow for further hard-filtered into SNPs and indels later. Hard-filtering using GATK's VariantFiltration tool is performed separately on the SNP and indel-specific project-level VCFs in ``rule hard_filter_snps`` and ``rule_hard_filter_indels``. 

