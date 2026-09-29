#!/bin/bash
#SBATCH --job-name=line_interkr
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/line_interkr_%j.log

# LINE embedding of the GM12878 INTERKR graph, SCI mode "both":
# order 1 and order 2 trained separately (paper hyperparameters: size 100, negative 5, samples 25M)
set -euo pipefail

LINE=/scratch/adas/genome-analysis/sci/LINE/line_hpc
GRAPH=/scratch/adas/genome-analysis/results/gm12878_interkr/gm12878_INTERKR_HiC_graph.txt
OUT=/scratch/adas/genome-analysis/results/gm12878_interkr

for ORDER in 1 2; do
  echo "=== order $ORDER start: $(date)"
  $LINE -train $GRAPH -order $ORDER -samples 25 -size 100 -negative 5 -rho 0.025 -threads 1 \
        -output $OUT/gm12878_INTERKR_order_${ORDER}_samples_25M.embedding
  echo "=== order $ORDER end: $(date)"
done
echo "=== ALL DONE: $(date)"
