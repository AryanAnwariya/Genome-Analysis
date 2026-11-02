#!/bin/bash
#SBATCH --job-name=ctcf_validation
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/ctcf_validation_%j.log

set -euo pipefail

PYTHON=/home/adas/.conda/envs/genome/bin/python
SCRIPT=/home/adas/sci-hic-validation/scripts/validate_chia_pet.py

BASE=/scratch/adas/genome-analysis
RESULTS=$BASE/results/chia_pet
CLUSTERS=$BASE/results/gm12878_interkr_oe

OUT=$RESULTS/full_ctcf_v1

mkdir -p "$OUT"

echo "=== FULL CTCF BIOLOGICAL VALIDATION ==="

"$PYTHON" "$SCRIPT" \
    --clusters-dir "$CLUSTERS" \
    --loops \
    "$RESULTS/4DNFIA3H6KXY.pet5.hg19.bedpe" \
    "$RESULTS/4DNFIS9CCN6R.pet5.hg19.bedpe" \
    --output "$OUT"

echo
echo "=== FINAL RESULTS ==="
cat "$OUT/summary.csv"

echo
echo "=== VALIDATION COMPLETED ==="
