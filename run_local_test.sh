#!/bin/bash

set -euo pipefail
snakemake -p \
        --notemp \
	--use-conda \
	--conda-frontend mamba \
	--rerun-incomplete \
	--jobs 1
