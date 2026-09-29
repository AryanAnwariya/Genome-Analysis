#!/bin/bash
#SBATCH --job-name=km_seeds
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/kmeans_seeds_%j.log
set -euo pipefail
cd /scratch/adas/genome-analysis
R=results/gm12878_interkr
mkdir -p $R/seeds
EV="/home/adas/.conda/envs/genome/bin/python scripts/evaluate_embedding.py --nodes $R/gm12878_INTERKR_nodes.tsv --rao reference/GSE63525_GM12878_subcompartments.bed.gz --sci sci/predictions/GM12878_SCI_sub_compartments.bed"
for SEED in 1 2 3 4 5; do
  $EV --seed $SEED --emb $R/gm12878_INTERKR_order_1_samples_25M.embedding $R/gm12878_INTERKR_order_2_samples_25M.embedding --out-prefix $R/seeds/LINE25M_scimerge_s$SEED
  $EV --seed $SEED --l2-each --emb $R/gm12878_INTERKR_order_1_samples_1000M.embedding $R/gm12878_INTERKR_order_2_samples_1000M.embedding --out-prefix $R/seeds/LINE1000M_l2each_s$SEED
  $EV --seed $SEED --emb $R/gm12878_INTERKR_weighted_deepwalk.embedding --out-prefix $R/seeds/DeepWalk_s$SEED
  $EV --seed $SEED --emb $R/gm12878_INTERKR_hope_katz.embedding --out-prefix $R/seeds/HOPE_s$SEED
done
echo "ALL DONE"
