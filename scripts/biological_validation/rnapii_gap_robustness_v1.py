#!/usr/bin/env python3
"""RNAPII sensitivity analysis: fixed bin-gap thresholds 2 and 3."""

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, "/home/adas/sci-hic-validation/scripts")

from validate_chia_pet import (
    BIN, METHODS, LABELS, assign_anchor, load_labels
)

ROOT = Path("/scratch/adas/genome-analysis")

INPUT = (
    ROOT / "data/chia_pet/prepared/rnapii_geo_v1/"
    "GSM1872887_RNAPII.pet5.hg19.bedpe"
)

OUTPUT = (
    ROOT / "results/chia_pet/"
    "rnapii_gap_robustness_v1.csv"
)

B = 20000
methods = list(METHODS)

labels = load_labels(
    ROOT / "results/gm12878_interkr_oe"
)

chroms = sorted({c for c, _ in labels})

assert len(chroms) == 22

enc = {
    m: {
        v: i for i, v in enumerate(LABELS[m])
    }
    for m in methods
}

arrays = {}

for chrom in chroms:

    positions = [
        s // BIN
        for c, s in labels
        if c == chrom
    ]

    arrays[chrom] = np.full(
        (max(positions) + 1, len(methods)),
        -1,
        dtype=np.int8
    )

for (chrom, start), record in labels.items():

    arrays[chrom][start // BIN] = [
        enc[m][record[m]]
        for m in methods
    ]

# Group observed interactions by chromosome
# and exact 100-kb bin separation.

group_n = Counter()

group_o = defaultdict(
    lambda: np.zeros(
        len(methods),
        dtype=np.int64
    )
)

input_n = 0
eligible = 0
same_bin = 0

with INPUT.open() as stream:

    for line in stream:

        fields = line.split()

        if len(fields) != 7:
            raise ValueError("Invalid BEDPE record")

        c1, s1, e1, c2, s2, e2, pet = fields

        s1, e1, s2, e2, pet = map(
            int, (s1, e1, s2, e2, pet)
        )

        input_n += 1

        if c1 != c2 or pet < 5:
            raise ValueError(
                "Unexpected unfiltered RNAPII record"
            )

        a = assign_anchor(
            c1, s1, e1, labels
        )

        b = assign_anchor(
            c2, s2, e2, labels
        )

        if a is None or b is None:
            continue

        eligible += 1

        gap = abs(a[1] - b[1]) // BIN

        if gap == 0:
            same_bin += 1
            continue

        key = (c1, gap)

        group_n[key] += 1

        group_o[key] += np.array(
            [
                int(
                    labels[a][m] == labels[b][m]
                )
                for m in methods
            ],
            dtype=np.int64
        )

# Verify against the earlier RNAPII analyses.

assert (input_n, eligible, same_bin) == (
    54175, 52848, 32889
)

assert sum(group_n.values()) == 19959

assert sum(
    n for (c, gap), n in group_n.items()
    if gap == 1
) == 15597

assert sum(
    n for (c, gap), n in group_n.items()
    if gap == 2
) == 2527

# Calculate the expected intra-compartment count
# for each chromosome and exact bin separation.

group_e = {}

for (chrom, gap), n in group_n.items():

    arr = arrays[chrom]

    left = arr[:-gap]
    right = arr[gap:]

    valid = (
        (left[:, 0] >= 0) &
        (right[:, 0] >= 0)
    )

    possible = int(valid.sum())

    if possible == 0:
        raise RuntimeError(
            f"No background pairs: {chrom}, gap={gap}"
        )

    probabilities = np.array([
        np.count_nonzero(
            valid & (left[:, j] == right[:, j])
        ) / possible
        for j in range(len(methods))
    ])

    group_e[(chrom, gap)] = n * probabilities


def odds_enrichment(observed, expected, count):

    with np.errstate(
        divide="ignore",
        invalid="ignore"
    ):

        return (
            observed / (count - observed)
        ) / (
            expected / (count - expected)
        )


rng = np.random.default_rng(20261102)

picks = rng.integers(
    0,
    len(chroms),
    size=(B, len(chroms))
)

rao = methods.index("Rao")

rows = []

print("=== RNAPII NONADJACENT-BIN ROBUSTNESS ===")
print("Common genomic bins:", len(labels))
print("Chromosomes:", len(chroms))
print("Bootstrap samples:", B)
print("Fixed thresholds: gap >= 2 and gap >= 3")
print("Bonferroni interval: 99.375% for 8 comparisons")

for threshold, expected_n in [
    (2, 4362),
    (3, 1835)
]:

    n = np.zeros(
        len(chroms),
        dtype=np.int64
    )

    o = np.zeros(
        (len(chroms), len(methods))
    )

    e = np.zeros_like(o)

    for i, chrom in enumerate(chroms):

        for (c, gap), count in group_n.items():

            if c == chrom and gap >= threshold:

                n[i] += count
                o[i] += group_o[(c, gap)]
                e[i] += group_e[(c, gap)]

    total = int(n.sum())

    assert total == expected_n

    full = odds_enrichment(
        o.sum(axis=0),
        e.sum(axis=0),
        total
    )

    # Resample whole chromosomes.

    bn = n[picks].sum(axis=1)[:, None]
    bo = o[picks].sum(axis=1)
    be = e[picks].sum(axis=1)

    bootstrap = odds_enrichment(
        bo, be, bn
    )

    if (
        not np.all(np.isfinite(bootstrap))
        or np.any(bootstrap <= 0)
    ):
        raise RuntimeError(
            f"Invalid bootstrap results for gap >= {threshold}"
        )

    print(f"\n=== GAP >= {threshold} ===")
    print("Interactions:", total)

    for j, method in enumerate(methods):

        relative = float(
            full[j] / full[rao]
        )

        low = None
        high = None

        # Four embedding methods across two thresholds:
        # eight Bonferroni-adjusted comparisons.

        if method not in ("Rao", "SCI_published"):

            sampled = (
                bootstrap[:, j] /
                bootstrap[:, rao]
            )

            low, high = map(
                float,
                np.percentile(
                    sampled,
                    [0.3125, 99.6875]
                )
            )

        print(f"\n{method}")
        print(f"  O/E odds: {full[j]:.4f}")
        print(f"  Relative to Rao: {relative:.4f}")

        if low is not None:

            print(
                f"  Adjusted 99.375% interval: "
                f"[{low:.4f}, {high:.4f}]"
            )

            print(
                f"  Entirely above 1: {low > 1}"
            )

        rows.append({
            "minimum_bin_gap": threshold,
            "loops": total,
            "method": method,
            "distance_matched_oe_odds": float(full[j]),
            "relative_to_rao": relative,
            "adjusted_99_375_low": (
                low if low is not None else ""
            ),
            "adjusted_99_375_high": (
                high if high is not None else ""
            ),
        })

OUTPUT.parent.mkdir(
    parents=True,
    exist_ok=True
)

# Never overwrite earlier results.

with OUTPUT.open("x", newline="") as stream:

    writer = csv.DictWriter(
        stream,
        fieldnames=list(rows[0])
    )

    writer.writeheader()
    writer.writerows(rows)

print("\nSaved:", OUTPUT)
print("RNAPII GAP ROBUSTNESS COMPLETED")
