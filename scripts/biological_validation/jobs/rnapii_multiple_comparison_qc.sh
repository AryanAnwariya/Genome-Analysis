#!/bin/bash
#SBATCH --job-name=rnapii_multi_qc
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/rnapii_multi_qc_%j.log

set -euo pipefail

export OPENBLAS_NUM_THREADS=1

/home/adas/.conda/envs/genome/bin/python - <<'PY'
import csv
from pathlib import Path
import numpy as np

BASE = Path("/scratch/adas/genome-analysis/results/chia_pet")
SOURCE = BASE / "rnapii_geo_distance_matched_v1/distance_matched_by_chromosome.csv"
OUT = BASE / "rnapii_multiple_comparison_v1.csv"

METHODS = ["Rao", "HOPE", "DeepWalk", "DW_HOPE", "LINE1000M_L2"]
B = 20000

if OUT.exists():
    raise SystemExit("Output exists; refusing to overwrite.")

with SOURCE.open() as f:
    rows = list(csv.DictReader(f))

data = {}

for row in rows:
    key = (row["chrom"], row["method"])
    if key in data:
        raise ValueError(f"Duplicate record: {key}")
    data[key] = row

chroms = sorted({r["chrom"] for r in rows})

assert len(chroms) == 22

n = np.zeros(len(chroms))
obs = np.zeros((len(chroms), len(METHODS)))
exp = np.zeros_like(obs)

for i, chrom in enumerate(chroms):
    for j, method in enumerate(METHODS):
        row = data[(chrom, method)]

        count = int(row["compared_loops"])
        if j == 0:
            n[i] = count
        else:
            assert count == n[i]

        obs[i, j] = int(row["observed_intra"])
        exp[i, j] = float(row["expected_intra"])

assert int(n.sum()) == 19959

def enrichment(o, e, total):
    return (o / (total - o)) / (e / (total - e))

full = enrichment(obs.sum(axis=0), exp.sum(axis=0), n.sum())

rng = np.random.default_rng(20261102)
picks = rng.integers(0, len(chroms), size=(B, len(chroms)))

bn = n[picks].sum(axis=1)[:, None]
bo = obs[picks].sum(axis=1)
be = exp[picks].sum(axis=1)

boot = enrichment(bo, be, bn)

assert np.all(np.isfinite(boot))

results = []

print("=== MULTIPLE-COMPARISON BOOTSTRAP ===")
print("Chromosomes:", len(chroms))
print("Bootstrap samples:", B)
print("Comparisons against Rao: 4")
print("Familywise alpha: 0.05")
print("Bonferroni interval level: 98.75%")

for j, method in enumerate(METHODS[1:], 1):
    point = full[j] / full[0]
    ratios = boot[:, j] / boot[:, 0]

    ci95 = np.percentile(ratios, [2.5, 97.5])
    ci_adjusted = np.percentile(ratios, [0.625, 99.375])

    print(f"\n{method}")
    print(f"Relative enrichment: {point:.4f}")
    print(f"Unadjusted 95% interval: {ci95[0]:.4f}–{ci95[1]:.4f}")
    print(f"Bonferroni 98.75% interval: {ci_adjusted[0]:.4f}–{ci_adjusted[1]:.4f}")
    print(f"Adjusted interval above 1: {ci_adjusted[0] > 1}")

    results.append({
        "method": method,
        "relative_to_rao": point,
        "ci95_low": ci95[0],
        "ci95_high": ci95[1],
        "adjusted_ci_low": ci_adjusted[0],
        "adjusted_ci_high": ci_adjusted[1],
        "adjusted_interval_above_1": ci_adjusted[0] > 1,
    })

with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(results[0]))
    writer.writeheader()
    writer.writerows(results)

print("\nSaved:", OUT)
print("MULTIPLE-COMPARISON QC COMPLETED")
PY
