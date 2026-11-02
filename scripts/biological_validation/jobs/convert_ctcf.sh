#!/bin/bash
#SBATCH --job-name=chia_ctcf_lift
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/chia_ctcf_lift_%j.log

set -euo pipefail

PYTHON=/home/adas/.conda/envs/genome/bin/python
SCRIPT=/home/adas/sci-hic-validation/scripts/liftover_bedpe.py
LIFTOVER=/home/adas/genome-tools/liftover/liftOver_v479
CHAIN=/home/adas/genome-tools/liftover/hg38ToHg19.over.chain.gz

BASE=/scratch/adas/genome-analysis/data/chia_pet
OUT=/scratch/adas/genome-analysis/results/chia_pet

mkdir -p "$OUT"

for ID in 4DNFIA3H6KXY 4DNFIS9CCN6R; do

    echo "======================================"
    echo "Processing $ID"
    echo "======================================"

    "$PYTHON" "$SCRIPT" \
        --input "$BASE/prepared/${ID}.pet5.bedpe" \
        --output "$OUT/${ID}.pet5.hg19.bedpe" \
        --audit "$OUT/${ID}.liftover_audit.json" \
        --liftover "$LIFTOVER" \
        --chain "$CHAIN"

done

echo "=== FINAL OUTPUT FILES ==="
wc -l "$OUT"/*.pet5.hg19.bedpe

echo "=== COMPLETED ==="
