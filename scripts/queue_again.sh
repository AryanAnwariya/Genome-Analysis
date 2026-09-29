#!/bin/bash
set -euo pipefail
cd /scratch/adas/genome-analysis
J1=$(sbatch --parsable scripts/run_line_o3_1000M.sh)
J2=$(sbatch --parsable scripts/run_line_o3_1000M_oe.sh)
J3=$(sbatch --parsable --dependency=afterok:${J1}:${J2} scripts/score_o3_1000M.sh)
J4=$(sbatch --parsable --dependency=afterok:${J3} scripts/compare_final.sh)
echo "queued: LINE INTERKR=$J1  LINE OE=$J2  score=$J3  compare=$J4"
