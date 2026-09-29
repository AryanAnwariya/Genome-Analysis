#!/bin/bash
#SBATCH --job-name=score_o3_oe
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/score_o3_oe_%j.log
set -euo pipefail
cd /scratch/adas/genome-analysis
EMB=results/gm12878_interkr_oe/gm12878_INTERKR_OE_order_3_samples_25M.embedding
RAO=reference/GSE63525_GM12878_subcompartments.bed.gz
SCI=sci/predictions/GM12878_SCI_sub_compartments.bed
set -- results/gm12878_interkr_oe/*_nodes.tsv
if [ "$#" -ne 1 ] || [ ! -e "$1" ]; then echo "need exactly one nodes file, found: $*"; exit 1; fi
NODES=$1
[ -e "$EMB" ] || { echo "embedding missing: $EMB"; exit 1; }
echo "nodes file: $NODES"
mkdir -p results/gm12878_interkr_oe/order3
for s in 1 2 3; do
  /home/adas/.conda/envs/genome/bin/python scripts/evaluate_embedding.py \
    --emb $EMB --nodes $NODES --rao $RAO --sci $SCI \
    --k 5 --seed $s --n-init 100 \
    --out-prefix results/gm12878_interkr_oe/order3/LINE_o3_25M_s$s > /dev/null
done
grep -H -E "kmeans inertia|cluster sizes|vs Rao|vs SCI" results/gm12878_interkr_oe/order3/*_metrics.txt
echo "=== done ==="
