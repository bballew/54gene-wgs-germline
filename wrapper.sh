#!/bin/sh

set -euo pipefail

DATE=$(date +"%Y%m%d%H%M")

f=$(awk '($0~/^sampleFile/){print $2}' config/config.yaml | sed "s/['\"]//g")
j=$(awk '($0~/^jobs/){print $2}' config/config.yaml)

snakemake -p \
	--use-conda \
	--conda-frontend mamba \
	--rerun-incomplete \
	--default-resources mem_mb=1024 batch=1 \
	--cluster "qsub -S /bin/bash -V -j y -o logs/ -cwd -pe mpi {threads} -l h_vmem={resources.mem_mb}M" \
	--resources batch=${j} \
	--jobs ${j} \
	--latency-wait 300 \
	&> WGS_${DATE}.out
