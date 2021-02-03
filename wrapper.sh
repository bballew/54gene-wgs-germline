#!/bin/sh

set -euo pipefail

DATE=$(date +"%Y%m%d%H%M")

snakemake -p --rerun-incomplete --cluster "qsub -S /bin/bash -V -j y -o logs -cwd -pe mpi {threads}" --jobs 60 --latency-wait 300 &> WGS_${DATE}.out"
