#!/bin/bash
#SBATCH --job-name=hic_archive_inventory
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/hic_archive_inventory_%j.log

set -euo pipefail

python3 - <<'PY'
import tarfile
import re
from pathlib import Path
from collections import defaultdict

archive = Path(
    "/scratch/adas/genome-analysis/data/"
    "GSE63525_GM12878_combined_interchromosomal_contact_matrices.tar.gz"
)

if not archive.is_file():
    raise FileNotFoundError(archive)

resolution_stats = defaultdict(lambda: [0, 0])
mapq_stats = defaultdict(lambda: [0, 0])
filetype_stats = defaultdict(lambda: [0, 0])
examples = defaultdict(list)

total_files = 0
total_bytes = 0

pattern = re.compile(
    r"(\d+(?:kb|mb))_resolution_interchromosomal"
)

print("=== FULL Hi-C ARCHIVE INVENTORY ===", flush=True)
print("Archive:", archive, flush=True)
print(
    "Compressed size:",
    round(archive.stat().st_size / 1024**3, 2),
    "GiB",
    flush=True
)

with tarfile.open(archive, "r|gz") as tar:

    for member in tar:

        if not member.isfile():
            continue

        total_files += 1
        total_bytes += member.size

        name = member.name

        match = pattern.search(name)

        resolution = (
            match.group(1)
            if match else "OTHER"
        )

        components = name.split("/")

        mapq = next(
            (
                part for part in components
                if part.startswith("MAPQ")
            ),
            "OTHER"
        )

        filetype = (
            name.rsplit(".", 1)[-1]
            if "." in name else "OTHER"
        )

        resolution_stats[resolution][0] += 1
        resolution_stats[resolution][1] += member.size

        mapq_stats[mapq][0] += 1
        mapq_stats[mapq][1] += member.size

        filetype_stats[filetype][0] += 1
        filetype_stats[filetype][1] += member.size

        if (
            filetype == "RAWobserved"
            and len(examples[resolution]) < 2
        ):
            examples[resolution].append(name)

def show_stats(title, stats):

    print(f"\n=== {title} ===")

    print(
        f"{'Category':18s}"
        f"{'Files':>12s}"
        f"{'Uncompressed GiB':>20s}"
    )

    for category, (count, size) in sorted(stats.items()):

        print(
            f"{category:18s}"
            f"{count:12,d}"
            f"{size / 1024**3:20.2f}"
        )

show_stats("RESOLUTIONS", resolution_stats)
show_stats("MAPQ CATEGORIES", mapq_stats)
show_stats("FILE TYPES", filetype_stats)

print("\n=== EXAMPLE RAW MATRICES ===")

for resolution, paths in sorted(examples.items()):

    print(f"\n{resolution}:")

    for path in paths:
        print(" ", path)

print("\n=== ARCHIVE TOTAL ===")
print("Total files:", f"{total_files:,}")
print(
    "Total uncompressed size:",
    round(total_bytes / 1024**3, 2),
    "GiB"
)

print("\nARCHIVE INVENTORY COMPLETED")
PY
