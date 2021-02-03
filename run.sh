#!/bin/sh

set -euo pipefail

qsub -S /bin/bash -V -j y -cwd wrapper.sh
