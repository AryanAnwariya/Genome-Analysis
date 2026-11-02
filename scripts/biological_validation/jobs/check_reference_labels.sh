#!/bin/bash
#SBATCH --job-name=chia_ref_qc
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/chia_ref_qc_%j.log

set -euo pipefail

/home/adas/.conda/envs/genome/bin/python - <<'PY'
import csv
import gzip
from bisect import bisect_right
from collections import defaultdict

base = "/scratch/adas/genome-analysis"

rao_file = (
    base + "/reference/"
    "GSE63525_GM12878_subcompartments.bed.gz"
)

sci_file = (
    base + "/sci/predictions/"
    "GM12878_SCI_sub_compartments.bed"
)

cluster_file = (
    base + "/results/gm12878_interkr_oe/"
    "eval/HOPE_s1_clusters.tsv"
)

# Load the ORIGINAL Rao BED intervals.

rao_intervals = defaultdict(list)

with gzip.open(rao_file, "rt") as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        label = fields[3]

        rao_intervals[chrom].append(
            (start, end, label)
        )

rao_starts = {}

for chrom in rao_intervals:
    rao_intervals[chrom].sort()
    rao_starts[chrom] = [
        row[0] for row in rao_intervals[chrom]
    ]

# Load the ORIGINAL published SCI BED intervals.

sci_labels = {}

with open(sci_file) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        key = (
            fields[0],
            int(fields[1]),
            int(fields[2])
        )

        if key in sci_labels:
            raise ValueError(
                f"Duplicate SCI interval: {key}"
            )

        sci_labels[key] = fields[3]


def original_rao_label(chrom, position):
    intervals = rao_intervals.get(chrom, [])
    starts = rao_starts.get(chrom, [])

    index = bisect_right(starts, position) - 1

    if index < 0:
        return None

    start, end, label = intervals[index]

    if start <= position < end:
        return label

    return None


checked = 0
rao_mismatches = []
sci_mismatches = []

with open(cluster_file) as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        if row["rao"] not in {
            "A1", "A2", "B1", "B2", "B3"
        }:
            continue

        if row["sci_pub"] not in {
            "C1", "C2", "C3", "C4", "C5"
        }:
            continue

        chrom = row["chrom"]
        start = int(row["start"])
        end = int(row["end"])

        checked += 1

        # Rao uses variable-length BED intervals.
        # Check the original Rao label at bin center.

        midpoint = (start + end) // 2

        expected_rao = original_rao_label(
            chrom, midpoint
        )

        if expected_rao != row["rao"]:
            rao_mismatches.append((
                chrom, start, end,
                row["rao"], expected_rao
            ))

        # Published SCI uses exact 100-kb bins.

        expected_sci = sci_labels.get(
            (chrom, start, end)
        )

        if expected_sci != row["sci_pub"]:
            sci_mismatches.append((
                chrom, start, end,
                row["sci_pub"], expected_sci
            ))

print("=== ORIGINAL REFERENCE LABEL QC ===")
print("Common bins checked:", checked)
print("Rao mismatches:", len(rao_mismatches))
print("Published SCI mismatches:", len(sci_mismatches))

if rao_mismatches:
    print("\nFirst Rao mismatches:")
    for row in rao_mismatches[:10]:
        print(row)

if sci_mismatches:
    print("\nFirst published SCI mismatches:")
    for row in sci_mismatches[:10]:
        print(row)

if checked != 26226:
    raise SystemExit(
        f"Unexpected common-bin count: {checked}"
    )

if rao_mismatches or sci_mismatches:
    raise SystemExit(
        "REFERENCE LABEL CONSISTENCY CHECK FAILED"
    )

print("\nALL REFERENCE LABEL CHECKS PASSED")
PY
