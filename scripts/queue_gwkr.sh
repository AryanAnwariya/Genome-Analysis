#!/bin/bash
set -euo pipefail
cd /scratch/adas/genome-analysis
RAO_DIR=/scratch/adas/genome-analysis/data/extracted/GM12878_combined_interchromosomal/100kb_resolution_interchromosomal
CHROM=/scratch/adas/genome-analysis/sci/chromosome_sizes/hg19.chrom.sizes
for f in scripts/build_sci_graph.py scripts/evaluate_embedding.py sci/LINE/line_hpc $CHROM; do
  test -e $f || { echo "MISSING $f"; exit 1; }
done
test -d $RAO_DIR || { echo "MISSING $RAO_DIR"; exit 1; }
mkdir -p results/gm12878_gwkr/eval logs

cat >> ~/genome-analysis-backup/DECISIONS.txt << 'EOF'

[GWKR input test, 2026-09-28, BEFORE result] Archive contains GWKRnorm vectors; authors' hic2sci.sh uses Juicer 'dump observed GW_KR' from the .hic. ASSUMED equivalent (not verified). Test: build graph with build_sci_graph.py --norm GWKRnorm (same rules: value>=1, ln weight), LINE at 25M, size 100, negative 5, threads 1: order 3, and orders 1+2 (concatenated SCI merge, and per-order L2). Score with 3 seeds, n_init=100, lowest inertia, against PUBLISHED SCI labels. Rule: ARI vs published >=0.7 matches in spirit, 0.4-0.7 partial, <0.4 no. Also record edge count (INTERKR: 227.5M undirected; old demo graph 8.35M, origin unknown) and nodes with edges. One graph, one pass, no tuning, reported whatever happens.
EOF

cat > scripts/build_gwkr.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=build_gwkr
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/build_gwkr_%j.log
set -euo pipefail
cd /scratch/adas/genome-analysis
/home/adas/.conda/envs/genome/bin/python scripts/build_sci_graph.py \
  --rao-dir /scratch/adas/genome-analysis/data/extracted/GM12878_combined_interchromosomal/100kb_resolution_interchromosomal \
  --chrom-sizes /scratch/adas/genome-analysis/sci/chromosome_sizes/hg19.chrom.sizes \
  --norm GWKRnorm --out-prefix results/gm12878_gwkr/gm12878_GWKR
echo "=== build done: $(date)"
EOF

cat > scripts/run_line_gwkr.sh << 'EOF'
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
EOF

cat > scripts/score_gwkr.sh << 'EOF'
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
EOF

J1=$(sbatch --parsable scripts/build_gwkr.sh)
J2=$(sbatch --parsable --dependency=afterok:${J1} scripts/run_line_gwkr.sh)
J3=$(sbatch --parsable --dependency=afterok:${J2} scripts/score_gwkr.sh)
echo "queued: build=$J1  line=$J2  score=$J3"
