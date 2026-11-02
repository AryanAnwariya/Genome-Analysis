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

BASE = Path("/scratch/adas/genome-analysis")
OUT = BASE / "results/chia_pet/distance_matched_v1"
OUT.mkdir(parents=True, exist_ok=True)

methods = list(METHODS)
encoders = {
    method: {label: i for i, label in enumerate(LABELS[method])}
    for method in methods
}

labels = load_labels(BASE / "results/gm12878_interkr_oe")

# Build chromosome arrays covering all common reference bins.
# -1 indicates that a genomic bin is outside the common set.

chrom_starts = defaultdict(list)

for chrom, start in labels:
    chrom_starts[chrom].append(start // BIN)

arrays = {}

for chrom, positions in chrom_starts.items():
    arrays[chrom] = np.full(
        (max(positions) + 1, len(methods)),
        -1,
        dtype=np.int16
    )

for (chrom, start), values in labels.items():
    index = start // BIN

    for j, method in enumerate(methods):
        arrays[chrom][index, j] = (
            encoders[method][values[method]]
        )

print("=== DISTANCE-MATCHED CTCF VALIDATION ===")
print("Common genomic bins:", len(labels))
print("Methods:", ", ".join(methods))

summary = []

datasets = [
    "4DNFIA3H6KXY",
    "4DNFIS9CCN6R",
]

expected_eligible = {
    "4DNFIA3H6KXY": 472,
    "4DNFIS9CCN6R": 67264,
}

for dataset in datasets:

    path = (
        BASE / "results/chia_pet" /
        f"{dataset}.pet5.hg19.bedpe"
    )

    eligible = 0
    same_bin = 0

    # Observed interactions at different chromosomes
    # and exact 100-kb bin separations.

    groups = Counter()
    observed_same = Counter()

    with path.open() as f:

        for line in f:

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) != 7:
                raise ValueError(
                    f"Invalid BEDPE line: {line[:100]}"
                )

            c1, s1, e1, c2, s2, e2, pet = fields

            s1, e1 = int(s1), int(e1)
            s2, e2 = int(s2), int(e2)
            pet = int(pet)

            if c1 != c2 or pet < 5:
                continue

            a = assign_anchor(c1, s1, e1, labels)
            b = assign_anchor(c2, s2, e2, labels)

            if a is None or b is None:
                continue

            eligible += 1

            if a == b:
                same_bin += 1
                continue

            distance_bins = abs(b[1] - a[1]) // BIN

            if distance_bins <= 0:
                raise AssertionError(
                    "Unexpected zero-distance bin pair"
                )

            groups[(c1, distance_bins)] += 1

            for method in methods:
                if labels[a][method] == labels[b][method]:
                    observed_same[method] += 1

    if eligible != expected_eligible[dataset]:
        raise AssertionError(
            f"{dataset}: eligible loops changed "
            f"from {expected_eligible[dataset]} "
            f"to {eligible}"
        )

    compared = eligible - same_bin

    assert sum(groups.values()) == compared

    # For each observed chromosome and bin separation,
    # calculate same-label probability among ALL
    # eligible genomic-bin pairs at that exact separation.
    #
    # Weight by how many observed CTCF loops have
    # that chromosome and separation.

    expected_same = defaultdict(float)

    for (chrom, distance_bins), nloops in groups.items():

        arr = arrays[chrom]

        left = arr[:-distance_bins]
        right = arr[distance_bins:]

        valid = (
            (left[:, 0] >= 0) &
            (right[:, 0] >= 0)
        )

        possible = int(np.count_nonzero(valid))

        if possible == 0:
            raise AssertionError(
                f"No background genomic pairs: "
                f"{chrom}, separation={distance_bins}"
            )

        for j, method in enumerate(methods):

            same = (
                (left[:, j] == right[:, j]) &
                valid
            )

            probability = (
                np.count_nonzero(same) / possible
            )

            expected_same[method] += (
                nloops * probability
            )

    print(f"\n=== {dataset} ===")
    print("Eligible loops:", eligible)
    print("Excluded same-bin loops:", same_bin)
    print("Compared distinct-bin loops:", compared)

    print(
        "\nMethod                 "
        "Observed%   Expected%   O/E odds"
    )

    for method in methods:

        observed = observed_same[method]
        expected = expected_same[method]

        observed_fraction = observed / compared
        expected_fraction = expected / compared

        observed_inter = compared - observed
        expected_inter = compared - expected

        observed_odds = (
            observed / observed_inter
            if observed_inter > 0 else float("inf")
        )

        expected_odds = (
            expected / expected_inter
            if expected_inter > 0 else float("inf")
        )

        enrichment = (
            observed_odds / expected_odds
            if expected_odds > 0 else float("nan")
        )

        print(
            f"{method:22s} "
            f"{100 * observed_fraction:8.2f} "
            f"{100 * expected_fraction:10.2f} "
            f"{enrichment:10.4f}"
        )

        summary.append({
            "dataset": dataset,
            "method": method,
            "eligible_loops": eligible,
            "same_bin_excluded": same_bin,
            "compared_loops": compared,
            "observed_intra": observed,
            "observed_inter": observed_inter,
            "observed_intra_fraction": observed_fraction,
            "expected_intra_count": expected,
            "expected_intra_fraction": expected_fraction,
            "distance_matched_odds_enrichment": enrichment,
        })

output = OUT / "distance_matched_summary.csv"

with output.open("w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=list(summary[0])
    )

    writer.writeheader()
    writer.writerows(summary)

print("\nSaved:", output)
print("DISTANCE-MATCHED ANALYSIS COMPLETED")
