#!/bin/bash
#SBATCH --job-name=ctcf_dist_qc
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/ctcf_dist_qc_%j.log

set -euo pipefail

/home/adas/.conda/envs/genome/bin/python - <<'PY'
import sys
from pathlib import Path
from collections import defaultdict

sys.path.insert(
    0, "/home/adas/sci-hic-validation/scripts"
)

from validate_chia_pet import load_labels, assign_anchor

base = Path("/scratch/adas/genome-analysis")

labels = load_labels(
    base / "results/gm12878_interkr_oe"
)

methods = [
    "Rao",
    "SCI_published",
    "HOPE",
    "DeepWalk",
    "DW_HOPE",
    "LINE1000M_L2",
]

ranges = [
    ("<500kb", 0, 500000),
    ("500kb-2Mb", 500000, 2000000),
    ("2Mb-10Mb", 2000000, 10000000),
    (">=10Mb", 10000000, float("inf")),
]

files = [
    "4DNFIA3H6KXY",
    "4DNFIS9CCN6R",
]

for file_id in files:

    path = (
        base / "results/chia_pet" /
        f"{file_id}.pet5.hg19.bedpe"
    )

    counts = defaultdict(
        lambda: defaultdict(
            lambda: [0, 0]
        )
    )

    eligible = 0
    same_bin = 0

    with path.open() as f:

        for line in f:

            c1, s1, e1, c2, s2, e2, pet = (
                line.rstrip("\n").split("\t")
            )

            s1, e1 = int(s1), int(e1)
            s2, e2 = int(s2), int(e2)

            key1 = assign_anchor(
                c1, s1, e1, labels
            )

            key2 = assign_anchor(
                c2, s2, e2, labels
            )

            if key1 is None or key2 is None:
                continue

            eligible += 1

            if key1 == key2:
                same_bin += 1

            midpoint1 = (s1 + e1) // 2
            midpoint2 = (s2 + e2) // 2

            distance = abs(
                midpoint2 - midpoint1
            )

            for name, low, high in ranges:

                if low <= distance < high:
                    group = name
                    break

            for method in methods:

                same = (
                    labels[key1][method]
                    ==
                    labels[key2][method]
                )

                counts[group][method][
                    0 if same else 1
                ] += 1

    print(f"\n=== {file_id} ===")
    print("Eligible loops:", eligible)
    print("Both anchors in same 100kb bin:", same_bin)

    for name, _, _ in ranges:

        total = sum(
            counts[name]["Rao"]
        )

        print(
            f"\nDistance: {name}, "
            f"loops: {total}"
        )

        for method in methods:

            intra, inter = counts[name][method]

            ratio = (
                intra / inter
                if inter else None
            )

            print(
                f"  {method:18s} "
                f"intra={intra:6d} "
                f"inter={inter:6d} "
                f"ratio={ratio}"
            )

print("\nDISTANCE QC COMPLETED")
PY
