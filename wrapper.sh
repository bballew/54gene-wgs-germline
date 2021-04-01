#!/bin/sh

set -euo pipefail

DATE=$(date +"%Y%m%d%H%M")

j=$(awk '($0~/^jobs/){print $2}' config/config.yaml)
q=$(awk '($0~/^default_queue/){print $2}' config/config.yaml)
cluster_mode=$(awk '($0~/^cluster_mode/){print $0}' config/config.yaml | sed "s/\"/'/g" | awk -F\' '($0~/^cluster_mode/){print $2}')

snakemake -p \
	--use-conda \
	--conda-frontend mamba \
	--rerun-incomplete \
	--restart-times 3 \
	--default-resources mem_mb=1024 batch=1 queue="${q}" \
	--cluster "${cluster_mode}" \
	--resources batch=${j} \
	--jobs ${j} \
	--latency-wait 300 \
	&> WGS_${DATE}.out
