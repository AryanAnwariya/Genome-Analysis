#!/bin/bash
#SBATCH --job-name=rnapii_dist_qc
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/rnapii_dist_qc_%j.log

set -euo pipefail

/home/adas/.conda/envs/genome/bin/python - <<'PY'
import sys
from pathlib import Path

sys.path.insert(
    0, "/home/adas/sci-hic-validation/scripts"
)

from validate_chia_pet import (
    load_labels, assign_anchor, METHODS
)

BASE = Path("/scratch/adas/genome-analysis")

INPUT = (
    BASE / "data/chia_pet/prepared/rnapii_geo_v1/"
    "GSM1872887_RNAPII.pet5.hg19.bedpe"
)

labels = load_labels(
    BASE / "results/gm12878_interkr_oe"
)

methods = list(METHODS)

ranges = [
    ("<500kb", 0, 500000),
    ("500kb-2Mb", 500000, 2000000),
    ("2Mb-10Mb", 2000000, 10000000),
    (">=10Mb", 10000000, float("inf")),
]

counts = {
    name: {
        method: [0, 0]
        for method in methods
    }
    for name, _, _ in ranges
}

total = 0
eligible = 0
same_bin = 0

with INPUT.open() as f:

    for line in f:
        if not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 7:
            raise ValueError("Invalid BEDPE record")

        c1, s1, e1, c2, s2, e2, pet = fields

        s1, e1 = int(s1), int(e1)
        s2, e2 = int(s2), int(e2)
        pet = int(pet)

        total += 1

        if c1 != c2 or pet < 5:
            raise ValueError("Unexpected input record")

        a = assign_anchor(c1, s1, e1, labels)
        b = assign_anchor(c2, s2, e2, labels)

        if a is None or b is None:
            continue

        eligible += 1

        if a == b:
            same_bin += 1

        midpoint_a = (s1 + e1) // 2
        midpoint_b = (s2 + e2) // 2

        distance = abs(midpoint_b - midpoint_a)

        for name, low, high in ranges:
            if low <= distance < high:
                group = name
                break

        for method in methods:

            same = (
                labels[a][method] ==
                labels[b][method]
            )

            counts[group][method][
                0 if same else 1
            ] += 1

assert total == 54175
assert eligible == 52848

print("\n=== RNAPII DISTANCE QC ===")
print("Input loops:", total)
print("Eligible loops:", eligible)
print("Unassigned loops:", total - eligible)
print("Same-bin loops:", same_bin)
print("Distinct-bin loops:", eligible - same_bin)

for name, _, _ in ranges:

    n = sum(counts[name]["Rao"])

    print(f"\n=== {name}: {n} loops ===")

    for method in methods:

        intra, inter = counts[group if False else name][method]

        ratio = intra / inter if inter else None

        print(
            f"{method:18s} "
            f"intra={intra:6d} "
            f"inter={inter:6d} "
            f"ratio={ratio}"
        )

print("\nRNAPII DISTANCE QC COMPLETED")
PY
