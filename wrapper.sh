#!/bin/sh

set -euo pipefail

DATE=$(date +"%Y%m%d%H%M")

j=$(awk '($0~/^jobs/){print $2}' config/config.yaml)
q=$(awk '($0~/^default_queue/){print $2}' config/config.yaml | sed "s/\"//g")
cluster_mode=$(awk '($0~/^cluster_mode/){print $0}' config/config.yaml | sed "s/\"/'/g" | awk -F\' '($0~/^cluster_mode/){print $2}')
t=$(awk '($0~/^tempDir/){print $2}' config/config.yaml | sed "s/\"//g")
echo "Parameters from config:"
echo "jobs: $j"
echo "default queue: $q"
echo "cluster submission string: $cluster_mode"
echo "temp directory: $t"
echo ""
v=$(snakemake --version)
echo "Snakemake version ${v}"
echo ""

CMD="snakemake -p \
	--use-conda \
	--conda-frontend mamba \
	--rerun-incomplete \
	--restart-times 3 \
	--default-resources mem_mb=1024 batch=1 queue=${q} \"tmpdir='$PWD/${t}'\" \
	--cluster \"${cluster_mode}\" \
	--resources batch=${j} \
	--jobs ${j} \
	--latency-wait 300 \
	&> WGS_${DATE}.out"

echo "Command run:"
echo $CMD
eval $CMD
