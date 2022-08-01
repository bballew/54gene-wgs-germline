# 54gene workflow: 54gene WGS germline

This workflow is designed to take either fastqs or gvcfs as input, and emit a joint-called multi-sample VCF.  Please see <insert url here> for additional documentation.

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

Please see the [documentation]() for more details.

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

To run an end-to-end test of the pipeline, please clone this repository and follow the instructions there: []()

Test cases are in the subfolder `tests`. They are automatically executed via continuous integration (TBD).

Unit tests for scripts are found in `workflow/scripts/`.  They can be executed with `pytest` or `testthat` for python and R scripts, respectively.

Test input data, configs, and example manifests for both pipeline modes can be found [here](https://gitlab.com/data-analysis5/54gene-wgs-test-data).  Note that there are a few important caveats when testing.
- Somalier doesn't seem to be functional on Mac, so make sure you're in a linux environment or comment out that target from the Snakefile (ugh).  (Using a container for the OS would solve this problem....)
- Make sure you update the manifests to point to wherever you put the input data
- You can run tests locally using the script in that same repo, via `bash run_local_test.sh`.
