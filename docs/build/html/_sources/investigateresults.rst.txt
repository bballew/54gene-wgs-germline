Investigate the Results
=======================

Assessing completion
--------------------

Upon pipeline completion, verify that all steps have completed without error by checking the top-level log ``WGS_<datestamp>.out``.  The bottom few lines of the file should contain something like nnn of nnn steps (100%) done.  Additional job logs (when run on a high-performance computing cluster) are stored in the ``logs/`` sub-directory.

Outputs and Results
-------------------

All pipeline results are stored in the ``results/`` directory.

- The hard-filtered, joint-called VCF can be found in ``results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz``
- For future joint-calling, the gVCFs are located at ``results/HaplotypeCaller/called/<sample>_all_chroms.g.vcf.gz``
- Deduplicated and post-BQSR bams are found at ``results/bqsr/<sample>.bam``


Samples that fail the following thresholds are automatically removed from the above file, and the output is placed in ``results/post_qc_exclusions/samples_excluded.HC_variants.hardfiltered.vcf.gz``:

1. Average depth of coverage < 20x
2. Contamination > 3%
3. Het/Hom ratio > 2.5

The record of sample exclusions, along with reasons for exclusion, is found at ``results/post_qc_exclusions/exclude_list_with_annotation.tsv``.  Samples are excluded for at least one of the following reasons.  Values listed are defaults, but can be changed in the ``config.yaml``.

QC
---

The following QC metrics are available:

- fastqc at ``results/fastqc/``
- Trimming stats via fastp at ``results/paired_trimmed_reads/``
- Alignment stats via samtools at ``results/alignment_stats/``
- Recalibration stats from bqsr at ``results/bqsr/``
- Relatedness via somalier at ``results/qc/relatedness/``
- Sample contamination via verifyBamID at ``results/qc/contamination_check/`` (for full runs; **not** included in joint-genotyping only)
- Inferred sex via bcftools +guess-ploidy at ``results/qc/sex_check/``
- Picard metrics at ``results/HaplotypeCaller/filtered/``
- bcftools stats at ``results/qc/bcftools_stats/``
- multiqc report at ``results/multiqc/``


