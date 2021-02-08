#!/bin/sh

set -euo pipefail

DATE=$(date +"%Y%m%d%H%M")

t=$(awk '($0~/^tempDir/){print $2}' config/config.yaml | sed "s/['\"]//g")
mkdir -p "${t}"
export TMPDIR="${t}"
export JDK_JAVA_OPTIONS=-Djava.io.tmpdir="${t}"
export _JAVA_OPTIONS=-Djava.io.tmpdir="${t}"
export IBM_JAVA_OPTIONS=-Djava.io.tmpdir="${t}"
mkdir -p logs/

snakemake -p \
	--use-conda \
	--conda-frontend mamba \
	--rerun-incomplete \
	--cluster "qsub -S /bin/bash -V -j y -o logs/ -cwd -pe mpi {threads}" \
	--jobs 40 \
	--latency-wait 300 \
	&> WGS_${DATE}.out
