#!/bin/bash
#SBATCH --job-name=label_continuity
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/label_continuity_%j.log

set -euo pipefail

/home/adas/.conda/envs/genome/bin/python - <<'PY'
import csv
from collections import Counter

path = (
    "/scratch/adas/genome-analysis/results/"
    "gm12878_interkr_oe/eval/"
    "HOPE_s1_clusters.tsv"
)

rows = []

with open(path) as f:
    for row in csv.DictReader(f, delimiter="\t"):

        if row["rao"] not in {"A1", "A2", "B1", "B2", "B3"}:
            continue

        if row["sci_pub"] not in {"C1", "C2", "C3", "C4", "C5"}:
            continue

        rows.append({
            "chrom": row["chrom"],
            "start": int(row["start"]),
            "Rao": row["rao"],
            "SCI": row["sci_pub"],
        })

rows.sort(key=lambda x: (x["chrom"], x["start"]))

assert len(rows) == 26226

print("=== LABEL CONTINUITY AUDIT ===")
print("Common bins:", len(rows))

for method in ["Rao", "SCI"]:

    counts = Counter(row[method] for row in rows)

    adjacent = 0
    same = 0

    for a, b in zip(rows, rows[1:]):

        if a["chrom"] != b["chrom"]:
            continue

        if b["start"] - a["start"] != 100000:
            continue

        adjacent += 1

        if a[method] == b[method]:
            same += 1

    n = len(rows)

    random_same = sum(
        count * (count - 1)
        for count in counts.values()
    ) / (n * (n - 1))

    print(f"\n=== {method} ===")
    print("Class sizes:", dict(sorted(counts.items())))
    print("Adjacent bin pairs checked:", adjacent)
    print("Adjacent pairs with same label:", same)

    print(
        "Adjacent same-label percentage:",
        round(100 * same / adjacent, 2)
    )

    print(
        "Random-pair same-label percentage:",
        round(100 * random_same, 2)
    )

print("\nLABEL CONTINUITY AUDIT COMPLETED")
PY
