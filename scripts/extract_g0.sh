#!/bin/bash
#SBATCH --job-name=extract_g0
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/extract_g0_%j.log
set -euo pipefail
cd /scratch/adas/genome-analysis
LIST=reference/extract_list_100kb_autosomes_MAPQG0.txt
grep -E "100kb_resolution_interchromosomal/chr([1-9]|1[0-9]|2[0-2])_chr([1-9]|1[0-9]|2[0-2])/MAPQG0/[^/]+$" reference/archive_listing.txt | awk '{print $NF}' > $LIST
n=$(wc -l < $LIST)
echo "list lines: $n"
if [ "$n" -ne 3465 ]; then echo "expected 3465 lines, stopping"; head -n 3 $LIST; exit 1; fi
if grep -q " " $LIST; then echo "list lines contain spaces, stopping"; head -n 2 $LIST; exit 1; fi
bad=$(grep -v -c "^GM12878_combined_interchromosomal/100kb_resolution_interchromosomal/" $LIST || true)
if [ "$bad" -ne 0 ]; then echo "$bad lines have an unexpected start, stopping"; grep -v "^GM12878_combined_interchromosomal/100kb_resolution_interchromosomal/" $LIST | head -n 3; exit 1; fi
head -n 2 $LIST
echo "=== extract start: $(date)"
tar -xzf data/GSE63525_GM12878_combined_interchromosomal_contact_matrices.tar.gz -C data/extracted_g0 -T $LIST
echo "=== extract end: $(date)"
find data/extracted_g0 -type f | wc -l
