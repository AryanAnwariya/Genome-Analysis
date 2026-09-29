#!/bin/bash
set -euo pipefail
cd /scratch/adas/genome-analysis
for f in scripts/build_sci_graph.py scripts/evaluate_embedding.py scripts/run_line_gwkr.sh scripts/score_gwkr.sh sci/LINE/line_hpc reference/archive_listing.txt data/GSE63525_GM12878_combined_interchromosomal_contact_matrices.tar.gz; do
  test -e $f || { echo "MISSING $f"; exit 1; }
done
mkdir -p results/gm12878_g0/eval logs data/extracted_g0 /scratch/adas/apptainer/cache /scratch/adas/apptainer/tmp

cat >> ~/genome-analysis-backup/DECISIONS.txt << 'EOF'

[MAPQ>=0 test + container, 2026-09-28, BEFORE results] Hypothesis (unverified, from memory): the paper's 4.9 billion reads may mean unfiltered data, so the authors may have used MAPQ>=0, not our MAPQ>=30. Archive has MAPQG0 (4048 listing entries = 253 pairs x 16). Test: extract 100kb autosomal MAPQG0 (expect 3465 files), build graph with --norm GWKRnorm --mapq MAPQG0 --write-7col, LINE 25M size 100 negative 5 threads 1 in orders 1, 2, 3, score against PUBLISHED labels (3 seeds, n_init=100). Rule: ARI vs published >=0.7 matches in spirit, 0.4-0.7 partial, <0.4 no. If no, stop reproduction attempts. Separately: pull authors' container (yuz12012/sci_container:latest) with apptainer to run their code on our 7-col input (not yet run). Observation: preliminary Docker LINE silhouette 0.006 / DB 9.09 resembles our collapsed concatenated LINE (0.005 / 8.3-8.8). README test run: -o both -s 1 on demo_data.txt; README default order is 1; utils.py single-order command string lacks a space (effect unknown).
EOF

cat > scripts/extract_g0.sh << 'EOF'
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
grep -E "100kb_resolution_interchromosomal/chr([1-9]|1[0-9]|2[0-2])_chr([1-9]|1[0-9]|2[0-2])/MAPQG0/[^/]+$" reference/archive_listing.txt > $LIST || true
n=$(wc -l < $LIST)
echo "list lines: $n"
if [ "$n" -ne 3465 ]; then echo "expected 3465 lines, stopping"; head -n 3 $LIST; exit 1; fi
if grep -q " " $LIST; then echo "list lines contain spaces, stopping"; head -n 2 $LIST; exit 1; fi
echo "=== extract start: $(date)"
tar -xzf data/GSE63525_GM12878_combined_interchromosomal_contact_matrices.tar.gz -C data/extracted_g0 -T $LIST
echo "=== extract end: $(date)"
find data/extracted_g0 -type f | wc -l
EOF

cat > scripts/build_g0.sh << 'EOF'
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
EOF

sed 's/gwkr/g0/g; s/GWKR/G0/g' scripts/run_line_gwkr.sh > scripts/run_line_g0.sh
sed 's/gwkr/g0/g; s/GWKR/G0/g' scripts/score_gwkr.sh > scripts/score_g0.sh

cat > scripts/pull_sci_container.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=pull_sci
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/pull_sci_%j.log
set -euo pipefail
export APPTAINER_CACHEDIR=/scratch/adas/apptainer/cache
export APPTAINER_TMPDIR=/scratch/adas/apptainer/tmp
cd /scratch/adas/apptainer
apptainer pull sci_container.sif docker://yuz12012/sci_container:latest
ls -lh sci_container.sif
EOF

grep -n -E "job-name|GRAPH=|OUT=|NODES=|^D=" scripts/run_line_g0.sh scripts/score_g0.sh
J1=$(sbatch --parsable scripts/extract_g0.sh)
J2=$(sbatch --parsable --dependency=afterok:${J1} scripts/build_g0.sh)
J3=$(sbatch --parsable --dependency=afterok:${J2} scripts/run_line_g0.sh)
J4=$(sbatch --parsable --dependency=afterok:${J3} scripts/score_g0.sh)
J5=$(sbatch --parsable scripts/pull_sci_container.sh)
echo "queued: extract=$J1 build=$J2 line=$J3 score=$J4 | container pull=$J5"
