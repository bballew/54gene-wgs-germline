#!/bin/sh

set -euo pipefail

sbatch -D ./ -o %x.o%j wrapper.sh
# qsub -S /bin/bash -V -j y -cwd wrapper.sh
