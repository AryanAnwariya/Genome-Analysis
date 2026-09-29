#!/bin/bash
#SBATCH --job-name=score_gwkr
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/score_gwkr_%j.log
set -euo pipefail
cd /scratch/adas/genome-analysis
D=results/gm12878_gwkr
NODES=$D/gm12878_GWKR_nodes.tsv
RAO=reference/GSE63525_GM12878_subcompartments.bed.gz
SCI=sci/predictions/GM12878_SCI_sub_compartments.bed
E1=$D/gm12878_GWKR_order_1_samples_25M.embedding
E2=$D/gm12878_GWKR_order_2_samples_25M.embedding
E3=$D/gm12878_GWKR_order_3_samples_25M.embedding
PY=/home/adas/.conda/envs/genome/bin/python
for s in 1 2 3; do
  $PY scripts/evaluate_embedding.py --emb $E3 --nodes $NODES --rao $RAO --sci $SCI --k 5 --seed $s --n-init 100 --out-prefix $D/eval/LINE_o3_s$s > /dev/null
  $PY scripts/evaluate_embedding.py --emb $E1 $E2 --nodes $NODES --rao $RAO --sci $SCI --k 5 --seed $s --n-init 100 --out-prefix $D/eval/LINE_concat_s$s > /dev/null
  $PY scripts/evaluate_embedding.py --emb $E1 $E2 --l2-each --nodes $NODES --rao $RAO --sci $SCI --k 5 --seed $s --n-init 100 --out-prefix $D/eval/LINE_l2each_s$s > /dev/null
done
grep -H -E "nodes clustered|cluster sizes|kmeans inertia|vs Rao|vs SCI" $D/eval/*_metrics.txt
echo "=== done ==="
