# 54gene workflow: 54gene WGS germline

This workflow is designed to take either fastqs or gvcfs as input, and emit a joint-called multi-sample VCF.  Please see [Read the Docs](https://54gene-wgs-germline.readthedocs.io/en/latest/) for additional documentation.

You can find a small test dataset and pre-configured files for this pipeline [here](https://gitlab.com/data-analysis5/dna-sequencing/54gene-wgs-test-data).

## Authors

* Esha Joshi
* Cameron Palmer
* Bari Jane Ballew (@bballew)

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository.

### Step 1: Obtain a copy of this workflow

Clone this repository to your local system, into the place where you want to perform the data analysis.
```
    git clone git@gitlab.com:data-analysis5/54gene-wgs-germline.git
```

### Step 2: Configure workflow

The pipeline inputs include:
- A configuration file
- A manifest file
- A list of intervals
- A sex linker file
- A MultiQC config file (provided)

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

    snakemake --use-conda --cluster sbatch --jobs 100


### Step 5: Investigate results

Upon pipeline completion, verify that all steps have completed without error by checking the top-level Snakemake log.  The bottom few lines of should contain something like `nnn of nnn steps (100%) done`.  Additional job logs (when run on a cluster) are stored in the `logs/` directory.

All pipeline results are stored in the `results/` directory.

The hard-filtered, joint-called VCF can be found in `results/HaplotypeCaller/filtered/HC_variants.hardfiltered.vcf.gz`.

For future joint-calling, the gVCFs are located at `results/HaplotypeCaller/called/<sample>_all_chroms.g.vcf.gz`.

Deduplicated and post-BQSR bams are found at `results/bqsr/<sample>.bam`.
