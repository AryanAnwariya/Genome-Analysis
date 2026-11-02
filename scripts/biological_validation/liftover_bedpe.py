#!/usr/bin/env python3
"""
Convert paired genomic intervals from BEDPE using UCSC liftOver.

- Preserves the relationship between both anchors.
- Preserves the original seventh-column support value.
- Rejects missing, duplicated, and cross-chromosome mappings.
- Checks coordinates and chromosome consistency.
- Records conversion statistics in JSON.
- Does not modify the original input file.

IMPORTANT:
This is provisional hg38-to-hg19 preprocessing.
The source assembly and interpretation of column 7 must be
confirmed from the individual 4DN file metadata before
biological results are reported.
"""

import argparse
import collections
import json
import subprocess
import tempfile
from pathlib import Path


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Convert BEDPE interactions using UCSC liftOver."
    )

    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--audit", required=True)
    parser.add_argument("--liftover", required=True)
    parser.add_argument("--chain", required=True)
    parser.add_argument("--min-match", type=float, default=0.95)

    return parser.parse_args()


def main():
    args = parse_arguments()

    if not 0 < args.min_match <= 1:
        raise ValueError("--min-match must be between 0 and 1")

    input_path = Path(args.input).resolve()
    output_path = Path(args.output).resolve()
    audit_path = Path(args.audit).resolve()

    liftover = Path(args.liftover).resolve()
    chain = Path(args.chain).resolve()

    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    if not liftover.is_file():
        raise FileNotFoundError(liftover)

    if not chain.is_file():
        raise FileNotFoundError(chain)

    if len({input_path, output_path, audit_path}) != 3:
        raise ValueError("Input, output and audit paths must differ")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.parent.mkdir(parents=True, exist_ok=True)

    stats = collections.Counter()
    records = []

    with tempfile.TemporaryDirectory(prefix="chia_liftover_") as tmp:
        tmp = Path(tmp)

        anchors_input = tmp / "anchors_hg38.bed"
        anchors_output = tmp / "anchors_hg19.bed"
        unmapped = tmp / "unmapped.bed"

        # STEP 1: Parse BEDPE and extract both anchors.

        with input_path.open() as src, anchors_input.open("w") as dst:

            for line_number, line in enumerate(src, 1):

                if not line.strip() or line.startswith("#"):
                    continue

                fields = line.rstrip("\n").split("\t")

                if len(fields) != 7:
                    raise ValueError(
                        f"Line {line_number}: expected 7 columns, "
                        f"found {len(fields)}"
                    )

                c1, s1, e1, c2, s2, e2, support = fields

                try:
                    s1, e1 = int(s1), int(e1)
                    s2, e2 = int(s2), int(e2)
                    support = int(support)
                except ValueError as exc:
                    raise ValueError(
                        f"Line {line_number}: invalid integer"
                    ) from exc

                if not (0 <= s1 < e1 and 0 <= s2 < e2):
                    raise ValueError(
                        f"Line {line_number}: invalid coordinates"
                    )

                if support < 0:
                    raise ValueError(
                        f"Line {line_number}: negative support"
                    )

                if c1 != c2:
                    stats["input_trans"] += 1

                index = len(records)

                records.append((
                    c1, s1, e1,
                    c2, s2, e2,
                    support
                ))

                dst.write(
                    f"{c1}\t{s1}\t{e1}\tloop{index}_A\n"
                )

                dst.write(
                    f"{c2}\t{s2}\t{e2}\tloop{index}_B\n"
                )

        stats["input_loops"] = len(records)
        stats["input_anchors"] = 2 * len(records)

        print(f"Input loops: {len(records):,}", flush=True)

        # STEP 2: Run UCSC liftOver.

        command = [
            str(liftover),
            f"-minMatch={args.min_match}",
            str(anchors_input),
            str(chain),
            str(anchors_output),
            str(unmapped),
        ]

        print("Running UCSC liftOver...", flush=True)

        subprocess.run(command, check=True)

        # STEP 3: Read converted anchors.

        mappings = collections.defaultdict(list)

        with anchors_output.open() as src:

            for line in src:

                fields = line.rstrip("\n").split("\t")

                if len(fields) < 4:
                    raise ValueError("Invalid liftOver output")

                chrom, start, end, identifier = fields[:4]

                start = int(start)
                end = int(end)

                mappings[identifier].append(
                    (chrom, start, end)
                )

        stats["mapped_anchor_records"] = sum(
            len(values) for values in mappings.values()
        )

        # STEP 4: Reconstruct interactions.

        converted = []

        for index, record in enumerate(records):

            c1, s1, e1, c2, s2, e2, support = record

            anchor_a = mappings.get(f"loop{index}_A", [])
            anchor_b = mappings.get(f"loop{index}_B", [])

            if not anchor_a or not anchor_b:
                stats["loops_missing_anchor"] += 1
                continue

            if len(anchor_a) != 1 or len(anchor_b) != 1:
                stats["loops_ambiguous_mapping"] += 1
                continue

            new_c1, new_s1, new_e1 = anchor_a[0]
            new_c2, new_s2, new_e2 = anchor_b[0]

            if new_s1 < 0 or new_s2 < 0:
                stats["loops_invalid_coordinates"] += 1
                continue

            if new_e1 <= new_s1 or new_e2 <= new_s2:
                stats["loops_invalid_coordinates"] += 1
                continue

            if new_c1 != new_c2:
                stats["loops_cross_chromosome"] += 1
                continue

            original_length_a = e1 - s1
            original_length_b = e2 - s2

            new_length_a = new_e1 - new_s1
            new_length_b = new_e2 - new_s2

            # Flag large interval-length changes.
            # This is an additional conservative QC rule,
            # separate from liftOver's minMatch threshold.

            if not (
                0.8 <= new_length_a / original_length_a <= 1.25
                and
                0.8 <= new_length_b / original_length_b <= 1.25
            ):
                stats["loops_large_length_change"] += 1
                continue

            converted.append((
                new_c1, new_s1, new_e1,
                new_c2, new_s2, new_e2,
                support
            ))

        stats["retained_loops"] = len(converted)

        # STEP 5: Write outputs only after successful processing.

        temporary_output = output_path.with_name(
            output_path.name + ".tmp"
        )

        with temporary_output.open("w") as dst:

            for record in converted:
                dst.write("\t".join(map(str, record)) + "\n")

        temporary_output.replace(output_path)

        # STEP 6: Produce audit report.

        unmapped_count = 0

        with unmapped.open() as src:
            for line in src:
                if line.startswith("#") or not line.strip():
                    continue
                unmapped_count += 1

        stats["unmapped_anchor_records"] = unmapped_count

        audit = {
            "input": str(input_path),
            "output": str(output_path),
            "chain": str(chain),
            "liftover": str(liftover),
            "min_match": args.min_match,
            "source_assembly": "hg38 (provisional)",
            "target_assembly": "hg19",
            "support_column": "column 7 (definition unverified)",
            "length_ratio_filter": [0.8, 1.25],
            "statistics": dict(stats),
        }

        with audit_path.open("w") as dst:
            json.dump(audit, dst, indent=2)

        print("\n=== CONVERSION SUMMARY ===")

        for key, value in sorted(stats.items()):
            print(f"{key}: {value:,}")

        print(f"\nOutput: {output_path}")
        print(f"Audit: {audit_path}")


if __name__ == "__main__":
    main()
