#!/bin/bash
#SBATCH --job-name=line_gwkr
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/line_gwkr_%j.log
set -euo pipefail
LINE=/scratch/adas/genome-analysis/sci/LINE/line_hpc
GRAPH=/scratch/adas/genome-analysis/results/gm12878_gwkr/gm12878_GWKR_HiC_graph.txt
OUT=/scratch/adas/genome-analysis/results/gm12878_gwkr
for ORDER in 1 2 3; do
  echo "=== order $ORDER start: $(date)"
  $LINE -train $GRAPH -order $ORDER -samples 25 -size 100 -negative 5 -rho 0.025 -threads 1 \
        -output $OUT/gm12878_GWKR_order_${ORDER}_samples_25M.embedding
  echo "=== order $ORDER end: $(date)"
done
