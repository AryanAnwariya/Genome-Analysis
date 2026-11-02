#!/bin/bash
#SBATCH --job-name=rnapii_gap_qc
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/rnapii_gap_qc_%j.log

set -euo pipefail

/home/adas/.conda/envs/genome/bin/python - <<'PY'
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(
    0, "/home/adas/sci-hic-validation/scripts"
)

from validate_chia_pet import BIN, load_labels, assign_anchor

base = Path("/scratch/adas/genome-analysis")

labels = load_labels(
    base / "results/gm12878_interkr_oe"
)

path = (
    base / "data/chia_pet/prepared/rnapii_geo_v1/"
    "GSM1872887_RNAPII.pet5.hg19.bedpe"
)

gaps = Counter()
eligible = 0
same_bin = 0

with path.open() as f:
    for line in f:
        c1, s1, e1, c2, s2, e2, pet = line.split()

        s1, e1, s2, e2, pet = map(
            int, (s1, e1, s2, e2, pet)
        )

        if c1 != c2 or pet < 5:
            raise ValueError("Unexpected filtered interaction")

        a = assign_anchor(c1, s1, e1, labels)
        b = assign_anchor(c2, s2, e2, labels)

        if a is None or b is None:
            continue

        eligible += 1

        gap = abs(a[1] - b[1]) // BIN
        gaps[gap] += 1

        if gap == 0:
            same_bin += 1

distinct = eligible - same_bin

assert eligible == 52848
assert same_bin == 32889
assert distinct == 19959

groups = [
    ("Adjacent bins (gap=1)", 1, 2),
    ("Gap=2", 2, 3),
    ("Gap=3-4", 3, 5),
    ("Gap=5-9", 5, 10),
    ("Gap>=10", 10, float("inf")),
]

print("=== RNAPII GENOMIC BIN-SEPARATION AUDIT ===")
print("Eligible loops:", eligible)
print("Same-bin loops:", same_bin)
print("Distinct-bin loops:", distinct)

total = 0

for name, low, high in groups:
    n = sum(
        count for gap, count in gaps.items()
        if low <= gap < high
    )

    total += n

    print(
        f"{name}: {n:,} "
        f"({100*n/distinct:.2f}% of distinct-bin loops)"
    )

assert total == distinct

print("\nTOP 10 EXACT BIN SEPARATIONS")

for gap, count in gaps.most_common(11):
    if gap == 0:
        continue
    print(f"Gap {gap}: {count:,} loops")

print("\nBIN-SEPARATION AUDIT COMPLETED")
PY
