#!/bin/bash
#SBATCH --job-name=pull_sci
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/pull_sci_%j.log
set -euo pipefail
export APPTAINER_CACHEDIR=/scratch/adas/apptainer/cache
export APPTAINER_TMPDIR=/scratch/adas/apptainer/tmp
cd /scratch/adas/apptainer
apptainer pull sci_container.sif docker://yuz12012/sci_container:latest
ls -lh sci_container.sif
