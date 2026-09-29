#!/bin/bash
#SBATCH --job-name=hope_interkr
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/hope_interkr_%j.log
set -euo pipefail
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
R=/scratch/adas/genome-analysis/results/gm12878_interkr
/home/adas/.conda/envs/genome/bin/python /scratch/adas/genome-analysis/scripts/hope_katz.py \
  --csr $R/gm12878_INTERKR_csr.npz \
  --out $R/gm12878_INTERKR_hope_katz.embedding
