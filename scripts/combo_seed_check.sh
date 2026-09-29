#!/bin/bash
#SBATCH --job-name=combo_seeds
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/combo_seeds_%j.log
set -euo pipefail
cd /scratch/adas/genome-analysis
R=results/gm12878_interkr
mkdir -p $R/seeds
EV="/home/adas/.conda/envs/genome/bin/python scripts/evaluate_embedding.py --nodes $R/gm12878_INTERKR_nodes.tsv --rao reference/GSE63525_GM12878_subcompartments.bed.gz --sci sci/predictions/GM12878_SCI_sub_compartments.bed"
HOPE=$R/gm12878_INTERKR_hope_katz.embedding
# DeepWalk embedding seeds (new walks + Word2Vec), K-means seed fixed at 12345
for S in 1 2 3; do
  DW=$R/gm12878_INTERKR_weighted_deepwalk_seed$S.embedding
  $EV --emb $DW --out-prefix $R/seeds/DeepWalk_embseed$S
  $EV --l2-each --emb $DW $HOPE --out-prefix $R/seeds/DW_HOPE_embseed$S
done
# K-means seeds on the original combination
for K in 1 2 3 4 5; do
  $EV --seed $K --l2-each --emb $R/gm12878_INTERKR_weighted_deepwalk.embedding $HOPE --out-prefix $R/seeds/DW_HOPE_kmseed$K
done
echo "ALL DONE"
