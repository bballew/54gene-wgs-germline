# 54gene workflow: 54gene WGS germline

This workflow is designed to take either fastqs or gvcfs as input, and emit a joint-called multi-sample VCF.  The pipeline analyzes germline samples only, and does not currently support multiplexed data.  It is optimized for deployment on AWS parallelcluster, which can be set up as described [here](https://gitlab.com/data-analysis5/gds-docs/parallelcluster_docs), though it should run without issue on any HPC system or local machine.

## Authors

* Bari Jane Ballew (@bballew)

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

### Step 1: Obtain a copy of this workflow

1. Clone this repository to your local system, into the place where you want to perform the data analysis.
```
    git clone git@gitlab.com:data-analysis5/54gene-wgs-germline.git
```

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `manifest.txt` to specify your sample setup.

#### A.  Config file

This pipeline currently offers two modes.  Please specify the run mode in `config.yaml`.
- __full__: This mode starts with fastq pairs and emits a joint-called, filtered, multi-sample VCF.
- __joint_genotyping__: This mode starts with gVCFs and runs joint-calling and filtering, emitting a multi-sample VCF.  The idea here is that samples can be analyzed in batches in the full run mode, and then the batches can be jointly re-genotyped with this mode.
- __fastq_qc_only__: This mode starts with fastq pairs and performs trimming, then runs FastQC on both pre- and post-trimming fastq pairs, and finally
  emits a MultiQC report.  This is designed to be run initially before a full pipeline run, so that if any lanes need to be dropped or trimming
  parameters adjusted, it can be identified before spending time and compute on alignment, joint calling, etc.

After joint-calling the samples listed in the manifest, the pipeline will create a new multi-sample VCF where samples have been automatically removed based on the following thresholds sourced from the config:
- max_het_ratio (defaults to 2.5): excludes samples with het/hom-alt ratios, as calculated by bcftools stats, that are above the threshold
- min_avg_depth (defaults to 20.0): excludes samples with average depth of coverage, as calclated by bcftools stats, that are below the threshold
- max_contam (defaults to 0.03): only applied in __full__ runs; excludes samples with contamination estimates, as calculated by verifyBamID, that are above the threshold

#### B.  Manifest file

You will need to provide a headerless, whitespace-delimited manifest file to run the pipeline.  For __full__ mode and for __fastq_qc_only__ mode:
- Columns: `readgroup (should be unique)  sample_ID   path/to/r1.fastq    path/to/r2.fastq`
- `readgroup` values should be unique, e.g. sampleID_flowcellID
- `sample_ID` should be the same for all fastq pairs from a single sample, and can be different from the fastq filenames

For __joint_genotyping__ mode:
- Columns: `sample_ID   path/to/file.g.vcf.gz`
- `sample_ID` values should be unique, and should correspond to the sample IDs in the gvcfs
- gvcfs should be zipped and indexed


### Step 3: Install the run-time environment

If needed, install miniconda by following the steps [here](https://docs.conda.io/en/latest/miniconda.html).

- Create a conda environment with, minimally, the dependencies defined in `environment.yaml`.

```
# create the env
conda env create -f environment.yaml
```

### Step 4: Execute workflow

Activate the conda environment:

    conda activate 54gene-wgs-germline

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

To run the pipeline in a cluster environment, edit `wrapper.sh` as needed for your system, and then run via

    bash run.sh

Alternatively, you may run snakemake pipelines on a cluster via something like this

    snakemake --use-conda --cluster qsub --jobs 100


### Step 5: Investigate results

Upon pipeline completion, verify that all steps have completed without error by checking the top-level log `WGS_<datestamp>.out`.  The bottom few lines of the file should contain something like `nnn of nnn steps (100%) done`.  Additional job logs (when run on a cluster) are stored in the `logs/` directory.

All pipeline results are stored in the `results/` directory.

The hard-filtered, joint-called VCF can be found in `results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz`.

For future joint-calling, the gVCFs are located at `results/HaplotypeCaller/called/<sample>_all_chroms.g.vcf.gz`.

Deduplicated and post-BQSR bams are found at `results/bqsr/<sample>.bam`.

Samples that fail the following thresholds are automatically removed from the above file, and the output is placed in `results/post_qc_exclusions/samples_excluded.HC_variants.hardfiltered.vcf.gz`.  The record of sample exclusions, along with reasons for exclusion, is found at `results/post_qc_exclusions/exclude_list_with_annotation.tsv`.  Samples are excluded for at least one of the following reasons.  Values listed are defaults, but can be changed in the `config.yaml`.
- Average depth of coverage <20x
- Contamination > 3%
- Het/Hom ratio > 2.5

The following QC metrics are available:
- fastqc at `results/fastqc/`
- Trimming stats via fastp at `results/paired_trimmed_reads/`
- Alignment stats via samtools at `results/alignment_stats/`
- Recalibration stats from bqsr at `results/bqsr/`
- Relatedness via somalier at `results/qc/relatedness/`
- Sample contamination via verifyBamID at `results/qc/contamination_check/` - for full runs; not included in joint-genotyping only
- Inferred sex via bcftools +guess-ploidy at `results/qc/sex_check/`
- Picard metrics at `results/HaplotypeCaller/filtered/`
- bcftools stats at `results/qc/bcftools_stats/`
- multiqc report at `results/multiqc/`

A final summary report for the run is available at `results/run_summary/`.  This report contains overview information, such as samples starting vs.
samples emitted into the final VCF, a table of key QC stats, and lists of samples flagged for unexpected relatedness, duplication, or sex discordance.


### Step 6: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push


### Step 7: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@github.com:snakemake-workflows/54gene-wgs-germline.git` or `git remote add -f upstream https://github.com/snakemake-workflows/54gene-wgs-germline.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/master workflow > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/master config`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.


### Step 8: Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local system, to a different place than where you ran your analysis.
3. Copy the modified files from your analysis to the clone of your fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not** accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary.
4. Commit and push your changes to your fork.
5. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against the original repository.


## Testing

Test cases are in the subfolder `tests`. They are automatically executed via continuous integration (TBD).

Unit tests for scripts are found in `workflow/scripts/`.  They can be executed with `pytest` or `testthat` for python and R scripts, respectively.

Test input data, configs, and example manifests for both pipeline modes can be found [here](https://gitlab.com/data-analysis5/54gene-wgs-test-data).  Note that there are a few important caveats when testing.
- Somalier doesn't seem to be functional on Mac, so make sure you're in a linux environment or comment out that target from the Snakefile (ugh).  (Using a container for the OS would solve this problem....)
- Make sure you update the manifests to point to wherever you put the input data
- You can run tests locally using the script in that same repo, via `bash run_local_test.sh`.
