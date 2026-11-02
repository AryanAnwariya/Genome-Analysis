#!/bin/bash
#SBATCH --job-name=chia_test
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/chia_test_%j.log

set -euo pipefail

PYTHON=/home/adas/.conda/envs/genome/bin/python
SCRIPT=/home/adas/sci-hic-validation/scripts/validate_chia_pet.py

BASE=/scratch/adas/genome-analysis
OUT=$BASE/results/chia_pet/smoke_test

mkdir -p "$OUT"

echo "=== CREATE REPRESENTATIVE TEST SAMPLE ==="

awk 'NR % 14 == 1' \
  "$BASE/results/chia_pet/4DNFIS9CCN6R.pet5.hg19.bedpe" \
  > "$OUT/ctcf_smoke.bedpe"

wc -l "$OUT/ctcf_smoke.bedpe"

echo "=== RUN VALIDATION ==="

"$PYTHON" "$SCRIPT" \
  --clusters-dir "$BASE/results/gm12878_interkr_oe" \
  --loops "$OUT/ctcf_smoke.bedpe" \
  --output "$OUT"

echo "=== VERIFY OUTPUT MATRICES ==="

"$PYTHON" - "$OUT/ctcf_smoke.validation.json" <<'PY'
import json
import sys

with open(sys.argv[1]) as f:
    result = json.load(f)

eligible = result["statistics"].get("common_eligible_loops", 0)

print("Common eligible loops:", eligible)

for method, values in result["methods"].items():
    matrix = values["symmetric_matrix"]
    normalized = values["normalized_matrix"]

    assert len(matrix) == 5
    assert all(len(row) == 5 for row in matrix)

    assert values["intra"] + values["inter"] == eligible

    for i in range(5):
        for j in range(5):
            assert matrix[i][j] == matrix[j][i]

    normalized_sum = sum(sum(row) for row in normalized)

    if eligible:
        assert abs(normalized_sum - 1.0) < 1e-10

    print(
        f"{method}: "
        f"intra={values['intra']}, "
        f"inter={values['inter']}, "
        f"ratio={values['intra_inter_ratio']}"
    )

print("\nALL OUTPUT CHECKS PASSED")
PY

echo "=== TEST COMPLETED ==="
