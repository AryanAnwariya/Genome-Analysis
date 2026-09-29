#!/bin/bash
#SBATCH --job-name=ninit100
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/ninit100_%j.log
set -euo pipefail
cd /scratch/adas/genome-analysis
R=results/gm12878_interkr
mkdir -p $R/ninit100
EV="/home/adas/.conda/envs/genome/bin/python scripts/evaluate_embedding.py --nodes $R/gm12878_INTERKR_nodes.tsv --rao reference/GSE63525_GM12878_subcompartments.bed.gz --sci sci/predictions/GM12878_SCI_sub_compartments.bed --n-init 100"
DW=$R/gm12878_INTERKR_weighted_deepwalk.embedding
HOPE=$R/gm12878_INTERKR_hope_katz.embedding
for K in 1 2 3; do
  $EV --seed $K --emb $HOPE --out-prefix $R/ninit100/HOPE_s$K
  $EV --seed $K --l2-each --emb $DW $HOPE --out-prefix $R/ninit100/DW_HOPE_s$K
  $EV --seed $K --emb $DW --out-prefix $R/ninit100/DeepWalk_s$K
  $EV --seed $K --l2-each --emb $R/gm12878_INTERKR_order_1_samples_1000M.embedding $R/gm12878_INTERKR_order_2_samples_1000M.embedding --out-prefix $R/ninit100/LINE1000M_l2each_s$K
  $EV --seed $K --emb $R/gm12878_INTERKR_order_1_samples_25M.embedding $R/gm12878_INTERKR_order_2_samples_25M.embedding --out-prefix $R/ninit100/LINE25M_scimerge_s$K
done
echo "ALL DONE"
