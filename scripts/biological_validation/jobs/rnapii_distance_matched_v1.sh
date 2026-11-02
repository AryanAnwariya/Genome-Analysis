#!/bin/bash
#SBATCH --job-name=rnapii_distmatch
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/rnapii_distmatch_%j.log

set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

/home/adas/.conda/envs/genome/bin/python \
  /home/adas/sci-hic-validation/scripts/rnapii_distance_matched_v1.py
