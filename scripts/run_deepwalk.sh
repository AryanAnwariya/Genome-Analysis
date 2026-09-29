#!/bin/bash
#SBATCH --job-name=dw_interkr
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/dw_interkr_%j.log
set -euo pipefail
export PYTHONHASHSEED=0   # needed for reproducible gensim results
R=/scratch/adas/genome-analysis/results/gm12878_interkr
/home/adas/.conda/envs/genome/bin/python /scratch/adas/genome-analysis/scripts/weighted_deepwalk_csr.py \
  --graph $R/gm12878_INTERKR_HiC_graph.txt \
  --csr-cache $R/gm12878_INTERKR_csr.npz \
  --out $R/gm12878_INTERKR_weighted_deepwalk.embedding
