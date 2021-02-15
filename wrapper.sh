#!/bin/sh

set -euo pipefail

DATE=$(date +"%Y%m%d%H%M")

t=$(awk '($0~/^tempDir/){print $2}' config/config.yaml | sed "s/['\"]//g")
f=$(awk '($0~/^sampleFile/){print $2}' config/config.yaml | sed "s/['\"]//g")
mkdir -p "${t}"
export TMPDIR="${t}"
#export JDK_JAVA_OPTIONS="-Xmx4g -Djava.io.tmpdir=${t} -XX:CompressedClassSpaceSize=200m"
export _JAVA_OPTIONS="-XX:CompressedClassSpaceSize=200m"
#export JAVA_OPTS="-Xmx4g -Djava.io.tmpdir=${t} -XX:CompressedClassSpaceSize=200m"
mkdir -p logs/

b=$(wc -l "${f}" | cut -d" " -f1) 

snakemake -p \
	--use-conda \
	--conda-frontend mamba \
	--rerun-incomplete \
	--default-resources mem_mb=1024 batch=1 \
	--cluster "qsub -S /bin/bash -V -j y -o logs/ -cwd -pe mpi {threads} -l h_vmem={resources.mem_mb}M" \
	--resources batch=${b} \
	--jobs 40 \
	--latency-wait 300 \
	&> WGS_${DATE}.out
