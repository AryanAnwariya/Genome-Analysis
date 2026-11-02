#!/usr/bin/env python3
"""
Chromosome-aware bootstrap for CTCF distance-matched enrichment.

Reuses:
- The original 26,226 common genomic bins
- Existing classifications
- Exact chromosome/bin-distance matching
- Same-bin exclusion

Checks results against distance_matched_v1 before bootstrapping.
"""

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(
    0, "/home/adas/sci-hic-validation/scripts"
)

from validate_chia_pet import (
    BIN, METHODS, LABELS, load_labels, assign_anchor
)

ROOT = Path("/scratch/adas/genome-analysis")
DATA = ROOT / "results/chia_pet"
OUT = DATA / "chromosome_bootstrap_v1"

OUT.mkdir(parents=True, exist_ok=True)

B = 4000

methods = list(METHODS)

labels = load_labels(
    ROOT / "results/gm12878_interkr_oe"
)

# Build the same common-bin genomic background
# used in the previous distance-matched analysis.

arrays = {}

for chrom in sorted({c for c, _ in labels}):

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

enc = {
    method: {
        label: i
        for i, label in enumerate(LABELS[method])
    }
    for method in methods
}

for (chrom, start), values in labels.items():

    arrays[chrom][start // BIN, :] = [
        enc[method][values[method]]
        for method in methods
    ]

# Read previous results for strict consistency checks.

prior = {}

with (
    DATA /
    "distance_matched_v1/distance_matched_summary.csv"
).open() as f:

    for row in csv.DictReader(f):
        prior[(row["dataset"], row["method"])] = row


datasets = [
    "4DNFIS9CCN6R",
    "4DNFIA3H6KXY",
]

for dindex, dataset in enumerate(datasets):

    groups = Counter()

    observed = defaultdict(
        lambda: np.zeros(len(methods), dtype=np.int64)
    )

    compared = Counter()

    eligible = 0
    same_bin = 0

    path = DATA / f"{dataset}.pet5.hg19.bedpe"

    with path.open() as f:

        for line in f:

            fields = line.split()

            if len(fields) != 7:
                raise ValueError("Invalid BEDPE record")

            c1, s1, e1, c2, s2, e2, support = fields

            s1, e1, s2, e2, support = map(
                int, (s1, e1, s2, e2, support)
            )

            if c1 != c2 or support < 5:
                continue

            a = assign_anchor(c1, s1, e1, labels)
            b = assign_anchor(c2, s2, e2, labels)

            if a is None or b is None:
                continue

            eligible += 1

            if a == b:
                same_bin += 1
                continue

            gap = abs(b[1] - a[1]) // BIN

            assert gap > 0

            groups[(c1, gap)] += 1
            compared[c1] += 1

            for j, method in enumerate(methods):

                observed[c1][j] += (
                    labels[a][method] ==
                    labels[b][method]
                )

    # Compute expected intra-loops separately
    # for each chromosome and exact bin separation.

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

        possible = int(valid.sum())

        if possible == 0:
            raise RuntimeError(
                f"No background pairs: {chrom}, {gap}"
            )

        for j in range(len(methods)):

            probability = np.count_nonzero(
                valid & (left[:, j] == right[:, j])
            ) / possible

            expected[chrom][j] += (
                nloops * probability
            )

    chroms = sorted(compared)

    counts = np.array(
        [compared[c] for c in chroms]
    )

    obs = np.stack(
        [observed[c] for c in chroms]
    )

    exp = np.stack(
        [expected[c] for c in chroms]
    )

    total = int(counts.sum())

    assert total == eligible - same_bin

    # Verify exact agreement with previous analysis.

    for j, method in enumerate(methods):

        previous = prior[(dataset, method)]

        assert eligible == int(
            previous["eligible_loops"]
        )

        assert same_bin == int(
            previous["same_bin_excluded"]
        )

        assert total == int(
            previous["compared_loops"]
        )

        assert int(obs[:, j].sum()) == int(
            previous["observed_intra"]
        )

        assert abs(
            float(exp[:, j].sum()) -
            float(previous["expected_intra_count"])
        ) < 1e-5

    print(f"\n=== {dataset} ===", flush=True)

    print("Eligible loops:", eligible)
    print("Same-bin excluded:", same_bin)
    print("Distinct-bin loops:", total)
    print("Chromosomes:", len(chroms))
    print("Previous results verified: YES")

    # Bootstrap whole chromosomes with replacement.

    rng = np.random.default_rng(
        20261102 + dindex
    )

    picks = rng.integers(
        0,
        len(chroms),
        size=(B, len(chroms))
    )

    boot_n = counts[picks].sum(axis=1)[:, None]
    boot_o = obs[picks].sum(axis=1)
    boot_e = exp[picks].sum(axis=1)

    with np.errstate(
        divide="ignore",
        invalid="ignore"
    ):

        boot_enrichment = (
            (boot_o / (boot_n - boot_o)) /
            (boot_e / (boot_n - boot_e))
        )

        full_o = obs.sum(axis=0)
        full_e = exp.sum(axis=0)

        full_enrichment = (
            (full_o / (total - full_o)) /
            (full_e / (total - full_e))
        )

    rao_index = methods.index("Rao")

    results = []

    print("\n=== CHROMOSOME BOOTSTRAP RESULTS ===")

    for j, method in enumerate(methods):

        samples = boot_enrichment[:, j]

        valid = samples[
            np.isfinite(samples) &
            (samples > 0)
        ]

        if len(valid) < B * 0.9:
            raise RuntimeError(
                f"Too many invalid samples: {method}"
            )

        low, high = np.percentile(
            valid, [2.5, 97.5]
        )

        # Paired method comparison against Rao.

        with np.errstate(
            divide="ignore",
            invalid="ignore"
        ):

            paired = (
                samples /
                boot_enrichment[:, rao_index]
            )

        paired = paired[
            np.isfinite(paired) &
            (paired > 0)
        ]

        if len(paired) < B * 0.9:
            raise RuntimeError(
                f"Too many invalid paired samples: {method}"
            )

        paired_low, paired_high = np.percentile(
            paired, [2.5, 97.5]
        )

        relative = (
            full_enrichment[j] /
            full_enrichment[rao_index]
        )

        results.append({
            "dataset": dataset,
            "method": method,
            "chromosomes": len(chroms),
            "distinct_bin_loops": total,
            "enrichment": float(full_enrichment[j]),
            "bootstrap_95_low": float(low),
            "bootstrap_95_high": float(high),
            "enrichment_relative_to_Rao": float(relative),
            "relative_to_Rao_95_low": float(paired_low),
            "relative_to_Rao_95_high": float(paired_high),
            "valid_bootstraps": len(valid),
        })

        print(
            f"{method:18s} "
            f"O/E={full_enrichment[j]:.4f} "
            f"[{low:.4f}, {high:.4f}] "
            f"| vs Rao={relative:.4f} "
            f"[{paired_low:.4f}, {paired_high:.4f}]"
        )

    output = OUT / f"{dataset}.bootstrap.csv"

    with output.open("w", newline="") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=list(results[0])
        )

        writer.writeheader()
        writer.writerows(results)

    print("\nSaved:", output)

print("\nCHROMOSOME BOOTSTRAP COMPLETED")
