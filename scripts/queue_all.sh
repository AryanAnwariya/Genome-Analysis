#!/bin/bash
set -euo pipefail
cd /scratch/adas/genome-analysis
for f in scripts/compare_vs_o3.py scripts/run_line_order3.sh scripts/run_line_order3_oe.sh; do
  test -f $f || { echo "MISSING $f - run the previous batch that creates it first"; exit 1; }
done

cat >> ~/genome-analysis-backup/DECISIONS.txt << 'EOF'

[Queued chain, 2026-09-28, BEFORE results] Last LINE variant: -order 3, 1000M samples, size 100, both graphs, -threads 4 (all earlier LINE runs used 1 thread; 4 used to fill the 4-CPU minimum and save time; multithreaded LINE is slightly non-deterministic). Scored with the usual protocol (3 seeds, n_init=100), then paired chromosome bootstrap against order 3 at 25M and 1000M. No further LINE variants will be run. Every LINE variant is reported and none dropped; the claim 'alternatives beat LINE' is stated per variant.
EOF

sed 's/-samples 25 /-samples 1000 /; s/-threads 1/-threads 4/; s/samples_25M/samples_1000M/; s/line_o3/line_o3_1000M/; s/line_order3_/line_order3_1000M_/' scripts/run_line_order3.sh > scripts/run_line_o3_1000M.sh
sed 's/-samples 25 /-samples 1000 /; s/-threads 1/-threads 4/; s/samples_25M/samples_1000M/; s/line_o3/line_o3_1000M/; s/line_order3_/line_order3_1000M_/' scripts/run_line_order3_oe.sh > scripts/run_line_o3_1000M_oe.sh
grep -n -E "job-name|--output|GRAPH=|samples|threads" scripts/run_line_o3_1000M.sh scripts/run_line_o3_1000M_oe.sh

python3 - << 'EOF'
s = open("scripts/compare_vs_o3.py").read()
s = s.replace("B = 1000", 'TAG = os.environ.get("O3TAG", "LINE_o3_25M"); LAB = os.environ.get("O3LAB", "LINE o3 25M"); DIRN = os.environ.get("O3DIR", "order3")\nB = 1000', 1)
s = s.replace('"results/gm12878_interkr/order3"', '"results/gm12878_interkr/" + DIRN').replace('"results/gm12878_interkr_oe/order3"', '"results/gm12878_interkr_oe/" + DIRN')
s = s.replace('"LINE_o3_25M"', "TAG").replace('"LINE o3 25M"', "LAB").replace("'LINE o3 25M'", "LAB")
s = s.replace('TAG = os.environ.get(TAG', 'TAG = os.environ.get("O3TAG"').replace('LAB = os.environ.get(LAB', 'LAB = os.environ.get("O3LAB"')
open("scripts/compare_vs_o3_param.py", "w").write(s)
EOF
/home/adas/.conda/envs/genome/bin/python -c "import ast; ast.parse(open('scripts/compare_vs_o3_param.py').read()); print('param compare script parses OK')"
grep -n -E "^TAG|DIRN" scripts/compare_vs_o3_param.py | head

cat > scripts/score_o3_1000M.sh << 'EOF'
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
EOF

cat > scripts/compare_final.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=cmp_final
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/cmp_final_%j.log
cd /scratch/adas/genome-analysis
PY=/home/adas/.conda/envs/genome/bin/python
O3TAG=LINE_o3_25M O3LAB="LINE o3 25M" O3DIR=order3 $PY scripts/compare_vs_o3_param.py > results/bio/compare_vs_o3_25M.txt 2>&1
O3TAG=LINE_o3_1000M O3LAB="LINE o3 1000M" O3DIR=order3_1000M $PY scripts/compare_vs_o3_param.py > results/bio/compare_vs_o3_1000M.txt 2>&1
echo "=== done ==="
EOF

J1=$(sbatch --parsable scripts/run_line_o3_1000M.sh)
J2=$(sbatch --parsable scripts/run_line_o3_1000M_oe.sh)
J3=$(sbatch --parsable --dependency=afterok:${J1}:${J2} scripts/score_o3_1000M.sh)
J4=$(sbatch --parsable --dependency=afterok:${J3} scripts/compare_final.sh)
echo "queued: LINE INTERKR=$J1  LINE OE=$J2  score=$J3  compare=$J4"
