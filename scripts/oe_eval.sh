#!/bin/bash
#SBATCH --job-name=oe_eval
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/oe_eval_%j.log
set -uo pipefail
cd /scratch/adas/genome-analysis
O=results/gm12878_interkr_oe; mkdir -p $O/eval
E=$O/gm12878_INTERKR_OE
EV="/home/adas/.conda/envs/genome/bin/python scripts/evaluate_embedding.py --nodes ${E}_nodes.tsv --rao reference/GSE63525_GM12878_subcompartments.bed.gz --sci sci/predictions/GM12878_SCI_sub_compartments.bed --n-init 100"
for K in 1 2 3; do
  $EV --seed $K --emb ${E}_order_1_samples_25M.embedding ${E}_order_2_samples_25M.embedding --out-prefix $O/eval/LINE25M_scimerge_s$K
  $EV --seed $K --l2-each --emb ${E}_order_1_samples_1000M.embedding ${E}_order_2_samples_1000M.embedding --out-prefix $O/eval/LINE1000M_l2each_s$K
  $EV --seed $K --emb ${E}_order_1_samples_1000M.embedding --out-prefix $O/eval/LINE_o1_1000M_s$K
  $EV --seed $K --emb ${E}_weighted_deepwalk.embedding --out-prefix $O/eval/DeepWalk_s$K
  $EV --seed $K --emb ${E}_hope_katz.embedding --out-prefix $O/eval/HOPE_s$K
  $EV --seed $K --l2-each --emb ${E}_weighted_deepwalk.embedding ${E}_hope_katz.embedding --out-prefix $O/eval/DW_HOPE_s$K
done
echo "ALL DONE"
