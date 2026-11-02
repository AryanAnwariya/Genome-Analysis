#!/usr/bin/env python3

import csv
import sys
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(
    0, "/home/adas/sci-hic-validation/scripts"
)

from validate_chia_pet import (
    BIN, METHODS, LABELS, load_labels, assign_anchor
)

ROOT = Path("/scratch/adas/genome-analysis")

INPUT = (
    ROOT / "data/chia_pet/prepared/rnapii_geo_v1/"
    "GSM1872887_RNAPII.pet5.hg19.bedpe"
)

RAW = (
    ROOT / "results/chia_pet/"
    "rnapii_geo_raw_v1/summary.csv"
)

OUT = (
    ROOT / "results/chia_pet/"
    "rnapii_geo_distance_matched_v1"
)

SUMMARY = OUT / "distance_matched_summary.csv"
BY_CHROM = OUT / "distance_matched_by_chromosome.csv"

if SUMMARY.exists() or BY_CHROM.exists():
    raise SystemExit(
        "ERROR: Output already exists. Refusing to overwrite."
    )

OUT.mkdir(parents=True, exist_ok=True)

methods = list(METHODS)

labels = load_labels(
    ROOT / "results/gm12878_interkr_oe"
)

encoders = {
    method: {
        label: i
        for i, label in enumerate(LABELS[method])
    }
    for method in methods
}

# Construct the same common-bin background
# used in our CTCF distance-matched analysis.

arrays = {}

chromosomes = sorted({
    chrom for chrom, start in labels
})

for chrom in chromosomes:

    positions = [
        start // BIN
        for c, start in labels
        if c == chrom
    ]

    arrays[chrom] = np.full(
        (max(positions) + 1, len(methods)),
        -1,
        dtype=np.int8
    )

for (chrom, start), values in labels.items():

    arrays[chrom][start // BIN, :] = [
        encoders[method][values[method]]
        for method in methods
    ]

print("=== RNAPII DISTANCE-MATCHED VALIDATION ===")
print("Common genomic bins:", len(labels))
print("Methods:", ", ".join(methods))

# Count observed interactions by chromosome
# and exact genomic-bin separation.

groups = Counter()
compared_chrom = Counter()

observed = defaultdict(
    lambda: np.zeros(len(methods), dtype=np.int64)
)

input_loops = 0
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

        s1, e1, s2, e2, pet = map(
            int, (s1, e1, s2, e2, pet)
        )

        input_loops += 1

        if c1 != c2 or pet < 5:
            raise ValueError(
                "Unexpected record in filtered RNAPII data"
            )

        a = assign_anchor(c1, s1, e1, labels)
        b = assign_anchor(c2, s2, e2, labels)

        if a is None or b is None:
            continue

        eligible += 1

        if a == b:
            same_bin += 1
            continue

        gap = abs(b[1] - a[1]) // BIN

        if gap <= 0:
            raise AssertionError(
                "Invalid distinct-bin separation"
            )

        groups[(c1, gap)] += 1
        compared_chrom[c1] += 1

        for j, method in enumerate(methods):

            if labels[a][method] == labels[b][method]:
                observed[c1][j] += 1

compared = eligible - same_bin

assert input_loops == 54175
assert eligible == 52848
assert same_bin == 32889
assert compared == 19959
assert sum(groups.values()) == compared

print("\n=== INPUT VERIFICATION ===")
print("Input loops:", input_loops)
print("Eligible loops:", eligible)
print("Same-bin excluded:", same_bin)
print("Distinct-bin loops:", compared)
print("Chromosomes represented:", len(compared_chrom))

# Calculate the expected number of same-label
# interactions from eligible genomic-bin pairs
# matched by chromosome and exact separation.

expected = defaultdict(
    lambda: np.zeros(len(methods), dtype=float)
)

for (chrom, gap), nloops in groups.items():

    arr = arrays[chrom]

    left = arr[:-gap]
    right = arr[gap:]

    valid = (
        (left[:, 0] >= 0) &
        (right[:, 0] >= 0)
    )

    possible = int(np.count_nonzero(valid))

    if possible == 0:
        raise RuntimeError(
            f"No matching background pairs: {chrom}, {gap}"
        )

    for j in range(len(methods)):

        same = (
            valid &
            (left[:, j] == right[:, j])
        )

        probability = (
            np.count_nonzero(same) / possible
        )

        expected[chrom][j] += (
            nloops * probability
        )

# Verify observed counts against the
# previous raw RNAPII validation.

with RAW.open() as f:
    raw_rows = {
        row["method"]: row
        for row in csv.DictReader(f)
    }

if set(raw_rows) != set(methods):
    raise AssertionError(
        "Methods differ from previous RNAPII analysis"
    )

total_observed = np.sum(
    list(observed.values()),
    axis=0
)

total_expected = np.sum(
    list(expected.values()),
    axis=0
)

for j, method in enumerate(methods):

    row = raw_rows[method]

    assert int(row["eligible_loops"]) == eligible

    assert (
        int(row["intra"])
        == int(total_observed[j]) + same_bin
    )

    assert (
        int(row["inter"])
        == compared - int(total_observed[j])
    )

print("Previous raw RNAPII results verified: YES")

# Calculate distance-matched enrichment.

summary_rows = []

print("\n=== DISTANCE-MATCHED RESULTS ===")

print(
    f"{'Method':20s} "
    f"{'Observed%':>10s} "
    f"{'Expected%':>10s} "
    f"{'O/E odds':>10s}"
)

for j, method in enumerate(methods):

    obs = int(total_observed[j])
    exp = float(total_expected[j])

    observed_inter = compared - obs
    expected_inter = compared - exp

    if (
        obs <= 0 or
        observed_inter <= 0 or
        exp <= 0 or
        expected_inter <= 0
    ):
        raise RuntimeError(
            f"Cannot calculate finite enrichment: {method}"
        )

    observed_fraction = obs / compared
    expected_fraction = exp / compared

    observed_odds = obs / observed_inter
    expected_odds = exp / expected_inter

    enrichment = observed_odds / expected_odds

    print(
        f"{method:20s} "
        f"{100 * observed_fraction:10.2f} "
        f"{100 * expected_fraction:10.2f} "
        f"{enrichment:10.4f}"
    )

    summary_rows.append({
        "dataset": "GSM1872887",
        "method": method,
        "eligible_loops": eligible,
        "same_bin_excluded": same_bin,
        "compared_loops": compared,
        "observed_intra": obs,
        "observed_inter": observed_inter,
        "observed_intra_fraction": observed_fraction,
        "expected_intra_count": exp,
        "expected_intra_fraction": expected_fraction,
        "distance_matched_odds_enrichment": enrichment,
    })

# Save per-chromosome contributions for later
# chromosome-aware bootstrap analysis.

chromosome_rows = []

for chrom in sorted(compared_chrom):

    for j, method in enumerate(methods):

        chromosome_rows.append({
            "chrom": chrom,
            "method": method,
            "compared_loops": compared_chrom[chrom],
            "observed_intra": int(observed[chrom][j]),
            "expected_intra": float(expected[chrom][j]),
        })

with SUMMARY.open("w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=list(summary_rows[0])
    )

    writer.writeheader()
    writer.writerows(summary_rows)

with BY_CHROM.open("w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=list(chromosome_rows[0])
    )

    writer.writeheader()
    writer.writerows(chromosome_rows)

print("\nSaved:", SUMMARY)
print("Saved:", BY_CHROM)

print("\nRNAPII DISTANCE-MATCHED ANALYSIS COMPLETED")
