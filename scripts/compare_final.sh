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
