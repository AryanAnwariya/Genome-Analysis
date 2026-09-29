#!/bin/bash
#SBATCH --job-name=build_g0
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/build_g0_%j.log
set -euo pipefail
cd /scratch/adas/genome-analysis
/home/adas/.conda/envs/genome/bin/python scripts/build_sci_graph.py \
  --rao-dir /scratch/adas/genome-analysis/data/extracted_g0/GM12878_combined_interchromosomal/100kb_resolution_interchromosomal \
  --chrom-sizes /scratch/adas/genome-analysis/sci/chromosome_sizes/hg19.chrom.sizes \
  --norm GWKRnorm --mapq MAPQG0 --write-7col --out-prefix results/gm12878_g0/gm12878_G0
echo "=== build done: $(date)"
