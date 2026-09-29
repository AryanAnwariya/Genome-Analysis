#!/bin/bash
#SBATCH --job-name=score_o3_1000M
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/score_o3_1000M_%j.log
set -euo pipefail
cd /scratch/adas/genome-analysis
RAO=reference/GSE63525_GM12878_subcompartments.bed.gz
SCI=sci/predictions/GM12878_SCI_sub_compartments.bed
set -- results/gm12878_interkr_oe/*_nodes.tsv
if [ "$#" -ne 1 ] || [ ! -e "$1" ]; then echo "need exactly one OE nodes file, found: $*"; exit 1; fi
NODES_O=$1
NODES_I=results/gm12878_interkr/gm12878_INTERKR_nodes.tsv
score() {
  mkdir -p $1/order3_1000M
  [ -e $2 ] || { echo "embedding missing: $2"; exit 1; }
  for s in 1 2 3; do
    /home/adas/.conda/envs/genome/bin/python scripts/evaluate_embedding.py --emb $2 --nodes $3 --rao $RAO --sci $SCI --k 5 --seed $s --n-init 100 --out-prefix $1/order3_1000M/LINE_o3_1000M_s$s > /dev/null
  done
}
score results/gm12878_interkr results/gm12878_interkr/gm12878_INTERKR_order_3_samples_1000M.embedding $NODES_I
score results/gm12878_interkr_oe results/gm12878_interkr_oe/gm12878_INTERKR_OE_order_3_samples_1000M.embedding $NODES_O
grep -H -E "kmeans inertia|cluster sizes|vs Rao|vs SCI" results/gm12878_interkr/order3_1000M/*_metrics.txt results/gm12878_interkr_oe/order3_1000M/*_metrics.txt
echo "=== done ==="
