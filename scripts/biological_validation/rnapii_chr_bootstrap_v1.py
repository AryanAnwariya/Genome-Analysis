#!/usr/bin/env python3

import csv
from pathlib import Path

import numpy as np

ROOT = Path("/scratch/adas/genome-analysis/results/chia_pet")

SOURCE = ROOT / "rnapii_geo_distance_matched_v1"

SUMMARY = SOURCE / "distance_matched_summary.csv"
BY_CHROM = SOURCE / "distance_matched_by_chromosome.csv"

OUT = ROOT / "rnapii_geo_chromosome_bootstrap_v1"
OUTPUT = OUT / "rnapii_bootstrap_summary.csv"

B = 4000
SEED = 20261102

METHODS = [
    "Rao",
    "SCI_published",
    "HOPE",
    "DeepWalk",
    "DW_HOPE",
    "LINE1000M_L2",
]

if OUTPUT.exists():
    raise SystemExit(
        "ERROR: Bootstrap results already exist. "
        "Refusing to overwrite."
    )

with SUMMARY.open() as f:
    previous = {
        row["method"]: row
        for row in csv.DictReader(f)
    }

with BY_CHROM.open() as f:
    records = list(csv.DictReader(f))

if set(previous) != set(METHODS):
    raise ValueError("Unexpected methods in previous results")

chroms = sorted(set(row["chrom"] for row in records))

if len(chroms) != 22:
    raise ValueError(
        f"Expected 22 chromosomes, found {len(chroms)}"
    )

lookup = {}

for row in records:

    key = (row["chrom"], row["method"])

    if key in lookup:
        raise ValueError(f"Duplicate record: {key}")

    lookup[key] = row

if len(lookup) != len(chroms) * len(METHODS):
    raise ValueError("Incomplete chromosome-method table")

n = np.zeros(len(chroms), dtype=np.int64)

obs = np.zeros(
    (len(chroms), len(METHODS)),
    dtype=np.float64
)

exp = np.zeros_like(obs)

for i, chrom in enumerate(chroms):

    for j, method in enumerate(METHODS):

        row = lookup[(chrom, method)]

        count = int(row["compared_loops"])
        observed = int(row["observed_intra"])
        expected = float(row["expected_intra"])

        if j == 0:
            n[i] = count
        elif count != n[i]:
            raise ValueError(
                f"Inconsistent loop counts: {chrom}"
            )

        if not (
            0 <= observed <= count
            and 0 <= expected <= count
        ):
            raise ValueError(
                f"Invalid counts: {chrom}, {method}"
            )

        obs[i, j] = observed
        exp[i, j] = expected

total = int(n.sum())

if total != 19959:
    raise ValueError(
        f"Expected 19,959 loops, found {total}"
    )

def odds_enrichment(observed, expected, count):

    with np.errstate(
        divide="ignore",
        invalid="ignore"
    ):

        return (
            (observed / (count - observed))
            /
            (expected / (count - expected))
        )

full_obs = obs.sum(axis=0)
full_exp = exp.sum(axis=0)

full_oe = odds_enrichment(
    full_obs,
    full_exp,
    total
)

# Verify agreement with job 429034.

for j, method in enumerate(METHODS):

    row = previous[method]

    if int(row["compared_loops"]) != total:
        raise AssertionError(
            f"Loop count mismatch: {method}"
        )

    if int(row["observed_intra"]) != int(full_obs[j]):
        raise AssertionError(
            f"Observed count mismatch: {method}"
        )

    if abs(
        float(row["expected_intra_count"]) -
        full_exp[j]
    ) > 1e-5:
        raise AssertionError(
            f"Expected count mismatch: {method}"
        )

    if abs(
        float(row["distance_matched_odds_enrichment"]) -
        full_oe[j]
    ) > 1e-8:
        raise AssertionError(
            f"Enrichment mismatch: {method}"
        )

print("=== RNAPII CHROMOSOME BOOTSTRAP ===")
print("Chromosomes:", len(chroms))
print("Distinct-bin loops:", total)
print("Bootstrap samples:", B)
print("Random seed:", SEED)
print("Previous results verified: YES")

# Resample chromosomes with replacement.
# Each resampled chromosome retains all its
# observed and expected loop contributions.

rng = np.random.default_rng(SEED)

indices = rng.integers(
    0,
    len(chroms),
    size=(B, len(chroms))
)

boot_n = n[indices].sum(axis=1)[:, None]
boot_obs = obs[indices].sum(axis=1)
boot_exp = exp[indices].sum(axis=1)

boot_oe = odds_enrichment(
    boot_obs,
    boot_exp,
    boot_n
)

if not np.all(np.isfinite(boot_oe)):
    raise RuntimeError(
        "Nonfinite bootstrap enrichment encountered"
    )

rao_idx = METHODS.index("Rao")
sci_idx = METHODS.index("SCI_published")

results = []

print()
print("=== BOOTSTRAP RESULTS ===")

for j, method in enumerate(METHODS):

    samples = boot_oe[:, j]

    low, high = np.percentile(
        samples,
        [2.5, 97.5]
    )

    # Paired comparison with Rao.

    relative_rao = (
        full_oe[j] / full_oe[rao_idx]
    )

    paired_rao = (
        boot_oe[:, j] /
        boot_oe[:, rao_idx]
    )

    rao_low, rao_high = np.percentile(
        paired_rao,
        [2.5, 97.5]
    )

    # Paired comparison with published SCI.

    relative_sci = (
        full_oe[j] / full_oe[sci_idx]
    )

    paired_sci = (
        boot_oe[:, j] /
        boot_oe[:, sci_idx]
    )

    sci_low, sci_high = np.percentile(
        paired_sci,
        [2.5, 97.5]
    )

    print()
    print(method)
    print(
        f"  Enrichment: {full_oe[j]:.4f} "
        f"[{low:.4f}, {high:.4f}]"
    )
    print(
        f"  Relative to Rao: {relative_rao:.4f} "
        f"[{rao_low:.4f}, {rao_high:.4f}]"
    )
    print(
        f"  Relative to SCI: {relative_sci:.4f} "
        f"[{sci_low:.4f}, {sci_high:.4f}]"
    )

    results.append({
        "method": method,
        "chromosomes": len(chroms),
        "distinct_bin_loops": total,
        "bootstrap_samples": B,
        "random_seed": SEED,
        "enrichment": float(full_oe[j]),
        "ci95_low": float(low),
        "ci95_high": float(high),
        "relative_to_rao": float(relative_rao),
        "rao_ci95_low": float(rao_low),
        "rao_ci95_high": float(rao_high),
        "relative_to_sci": float(relative_sci),
        "sci_ci95_low": float(sci_low),
        "sci_ci95_high": float(sci_high),
    })

OUT.mkdir(parents=True, exist_ok=True)

with OUTPUT.open("w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=list(results[0])
    )

    writer.writeheader()
    writer.writerows(results)

print()
print("Saved:", OUTPUT)
print("RNAPII CHROMOSOME BOOTSTRAP COMPLETED")
