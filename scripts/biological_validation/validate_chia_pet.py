#!/usr/bin/env python3
"""
GM12878 ChIA-PET validation against Rao, published SCI,
and four pair-O/E graph embedding classifications.

Methods:
- PET >= 5 (already filtered)
- Common 26,226-bin reference set
- Anchor overlap strictly > 50%
- Common eligible interactions for every classification
- Unweighted intra/inter-loop ratio
- Symmetric 5x5 interaction matrix normalized by total
- Replicates analyzed separately

IMPORTANT:
Results remain PROVISIONAL until original ChIA-PET file
assembly and column-7 PET-support definition are verified.
"""

import csv
import json
import argparse
from pathlib import Path
from collections import Counter

BIN = 100000

METHODS = {
    "Rao": None,
    "SCI_published": None,
    "HOPE": ("HOPE", 1),
    "DeepWalk": ("DeepWalk", 1),
    "DW_HOPE": ("DW_HOPE", 1),
    "LINE1000M_L2": ("LINE1000M_l2each", 3),
}

LABELS = {
    "Rao": ["A1", "A2", "B1", "B2", "B3"],
    "SCI_published": ["C1", "C2", "C3", "C4", "C5"],
}

for method in METHODS:
    if method not in LABELS:
        LABELS[method] = ["0", "1", "2", "3", "4"]


def load_labels(base):
    """
    Load genomic coordinates and all classifications.

    The cluster files contain the original Rao and SCI
    labels, so no independent node-ID reconstruction
    is necessary.
    """

    all_labels = {}
    common_keys = None

    for method, spec in METHODS.items():

        if spec is None:
            continue

        name, seed = spec

        path = (
            base / "eval" /
            f"{name}_s{seed}_clusters.tsv"
        )

        if not path.exists():
            raise FileNotFoundError(path)

        current = {}

        with path.open() as f:
            reader = csv.DictReader(f, delimiter="\t")

            required = {
                "chrom", "start", "end",
                "cluster", "rao", "sci_pub"
            }

            if not required.issubset(reader.fieldnames or []):
                raise ValueError(
                    f"Missing columns in {path}"
                )

            for row in reader:

                chrom = row["chrom"]
                start = int(row["start"])
                end = int(row["end"])

                if end - start != BIN or start % BIN:
                    raise ValueError(
                        f"Unexpected bin coordinates: {row}"
                    )

                if row["rao"] not in LABELS["Rao"]:
                    continue

                if row["sci_pub"] not in LABELS["SCI_published"]:
                    continue

                cluster = row["cluster"]

                if cluster not in LABELS[method]:
                    raise ValueError(
                        f"Invalid cluster label: {cluster}"
                    )

                key = (chrom, start)

                if key in current:
                    raise ValueError(
                        f"Duplicate genomic bin: {key}"
                    )

                current[key] = row

        keys = set(current)

        if common_keys is None:
            common_keys = keys
        else:
            common_keys &= keys

        all_labels[method] = current

    if len(common_keys) != 26226:
        raise ValueError(
            "Expected exactly 26,226 common bins; "
            f"found {len(common_keys)}"
        )

    # Verify reference labels are identical across files.

    first_method = next(iter(all_labels))

    reference = all_labels[first_method]

    for method, rows in all_labels.items():
        for key in common_keys:

            a = reference[key]
            b = rows[key]

            if (
                a["rao"] != b["rao"] or
                a["sci_pub"] != b["sci_pub"]
            ):
                raise ValueError(
                    f"Reference mismatch: {method}, {key}"
                )

    final = {}

    for key in common_keys:

        row = reference[key]

        values = {
            "Rao": row["rao"],
            "SCI_published": row["sci_pub"],
        }

        for method, rows in all_labels.items():
            values[method] = rows[key]["cluster"]

        final[key] = values

    print(
        f"Loaded {len(final):,} common genomic bins",
        flush=True
    )

    return final


def assign_anchor(chrom, start, end, labels):
    """
    Assign anchor if a single eligible 100-kb bin
    overlaps strictly more than half the anchor.

    Multiple neighboring bins are checked when
    an anchor crosses bin boundaries.
    """

    if start < 0 or end <= start:
        return None

    length = end - start

    first = start // BIN
    last = (end - 1) // BIN

    best_key = None
    best_overlap = 0

    for index in range(first, last + 1):

        bin_start = index * BIN
        bin_end = bin_start + BIN

        overlap = max(
            0,
            min(end, bin_end) - max(start, bin_start)
        )

        key = (chrom, bin_start)

        if key not in labels:
            continue

        if overlap > best_overlap:
            best_overlap = overlap
            best_key = key

    if best_overlap * 2 <= length:
        return None

    return best_key


def evaluate_file(path, labels):
    """
    Compute method-specific matrices on a common
    set of eligible loops.
    """

    matrices = {}

    for method in METHODS:
        matrices[method] = [
            [0 for _ in range(5)]
            for _ in range(5)
        ]

    index_maps = {
        method: {
            label: index
            for index, label in enumerate(LABELS[method])
        }
        for method in METHODS
    }

    stats = Counter()

    with path.open() as f:

        for line_number, line in enumerate(f, 1):

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 7:
                raise ValueError(
                    f"{path}:{line_number}: "
                    "expected 7 columns"
                )

            c1, s1, e1, c2, s2, e2, support = fields

            s1, e1 = int(s1), int(e1)
            s2, e2 = int(s2), int(e2)
            support = int(support)

            stats["input_loops"] += 1

            if support < 5:
                stats["below_pet_threshold"] += 1
                continue

            if c1 != c2:
                stats["trans_loops"] += 1
                continue

            key_a = assign_anchor(c1, s1, e1, labels)
            key_b = assign_anchor(c2, s2, e2, labels)

            if key_a is None or key_b is None:
                stats["missing_assignment"] += 1
                continue

            stats["common_eligible_loops"] += 1

            for method in METHODS:

                label_a = labels[key_a][method]
                label_b = labels[key_b][method]

                i = index_maps[method][label_a]
                j = index_maps[method][label_b]

                matrices[method][i][j] += 1

    results = {}

    for method, matrix in matrices.items():

        total = sum(
            sum(row) for row in matrix
        )

        intra = sum(
            matrix[i][i] for i in range(5)
        )

        inter = total - intra

        # SCI's symmetric matrix M + M.T.

        symmetric = [
            [
                matrix[i][j] + matrix[j][i]
                for j in range(5)
            ]
            for i in range(5)
        ]

        symmetric_total = sum(
            sum(row) for row in symmetric
        )

        normalized = [
            [
                value / symmetric_total
                if symmetric_total else 0.0
                for value in row
            ]
            for row in symmetric
        ]

        results[method] = {
            "labels": LABELS[method],
            "intra": intra,
            "inter": inter,
            "total": total,
            "intra_fraction": (
                intra / total if total else None
            ),
            "intra_inter_ratio": (
                intra / inter if inter else None
            ),
            "matrix": matrix,
            "symmetric_matrix": symmetric,
            "normalized_matrix": normalized,
        }

        if total != stats["common_eligible_loops"]:
            raise AssertionError(
                f"Loop count mismatch for {method}"
            )

    return dict(stats), results


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--clusters-dir",
        required=True
    )

    parser.add_argument(
        "--loops",
        nargs="+",
        required=True
    )

    parser.add_argument(
        "--output",
        required=True
    )

    args = parser.parse_args()

    cluster_dir = Path(args.clusters_dir)
    output_dir = Path(args.output)

    output_dir.mkdir(parents=True, exist_ok=True)

    labels = load_labels(cluster_dir)

    summary_rows = []

    for loop_file in args.loops:

        path = Path(loop_file)

        print(
            f"\nEvaluating {path.name}",
            flush=True
        )

        stats, results = evaluate_file(path, labels)

        name = path.name.replace(".bedpe", "")

        report = {
            "file": str(path),
            "status": "PROVISIONAL",
            "source_assembly": (
                "hg38 assumed; converted to hg19"
            ),
            "pet_column": (
                "Column 7 assumed PET support"
            ),
            "common_bins": len(labels),
            "anchor_overlap_rule": (
                "strictly greater than 50% "
                "of anchor length"
            ),
            "counts_are_pet_weighted": False,
            "statistics": stats,
            "methods": results,
        }

        json_path = output_dir / f"{name}.validation.json"

        with json_path.open("w") as f:
            json.dump(report, f, indent=2)

        print("Statistics:", stats)

        for method, values in results.items():

            print(
                method,
                "intra=", values["intra"],
                "inter=", values["inter"],
                "ratio=", values["intra_inter_ratio"]
            )

            summary_rows.append({
                "file": name,
                "method": method,
                "eligible_loops": stats[
                    "common_eligible_loops"
                ],
                "intra": values["intra"],
                "inter": values["inter"],
                "intra_fraction": values[
                    "intra_fraction"
                ],
                "intra_inter_ratio": values[
                    "intra_inter_ratio"
                ],
            })

    summary_path = output_dir / "summary.csv"

    with summary_path.open("w", newline="") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=[
                "file",
                "method",
                "eligible_loops",
                "intra",
                "inter",
                "intra_fraction",
                "intra_inter_ratio",
            ]
        )

        writer.writeheader()
        writer.writerows(summary_rows)

    print("\nSaved:", summary_path)


if __name__ == "__main__":
    main()
