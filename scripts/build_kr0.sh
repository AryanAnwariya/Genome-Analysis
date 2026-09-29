#!/bin/bash
#SBATCH --job-name=build_kr0
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/build_kr0_%j.log
set -euo pipefail
cd /scratch/adas/genome-analysis
/home/adas/.conda/envs/genome/bin/python scripts/build_sci_graph.py \
  --rao-dir /scratch/adas/genome-analysis/data/extracted_g0/GM12878_combined_interchromosomal/100kb_resolution_interchromosomal \
  --chrom-sizes /scratch/adas/genome-analysis/sci/chromosome_sizes/hg19.chrom.sizes \
  --norm KRnorm --mapq MAPQG0 --out-prefix results/gm12878_kr0/gm12878_KR0
echo "=== build done: $(date)"
