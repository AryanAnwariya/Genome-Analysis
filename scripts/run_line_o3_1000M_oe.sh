#!/bin/bash
#SBATCH --job-name=line_o3_1000M_oe
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/line_order3_1000M_oe_%j.log
set -euo pipefail
LINE=/scratch/adas/genome-analysis/sci/LINE/line_hpc
GRAPH=/scratch/adas/genome-analysis/results/gm12878_interkr_oe/gm12878_INTERKR_OE_HiC_graph.txt
OUT=/scratch/adas/genome-analysis/results/gm12878_interkr_oe
echo "=== order 3 start: $(date)"
$LINE -train $GRAPH -order 3 -samples 1000 -size 100 -negative 5 -rho 0.025 -threads 1 \
      -output $OUT/gm12878_INTERKR_OE_order_3_samples_1000M.embedding
echo "=== order 3 end: $(date)"
