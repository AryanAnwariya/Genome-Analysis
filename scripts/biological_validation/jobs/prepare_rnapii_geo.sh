#!/bin/bash
#SBATCH --job-name=rnapii_prepare
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/rnapii_prepare_%j.log

set -euo pipefail

BASE=/scratch/adas/genome-analysis/data/chia_pet

INPUT="$BASE/geo_GSE72816/GSM1872887_GM12878_RNAPII_PET_clusters.txt.gz"
OUTDIR="$BASE/prepared/rnapii_geo_v1"
OUTPUT="$OUTDIR/GSM1872887_RNAPII.pet5.hg19.bedpe"

mkdir -p "$OUTDIR"

if [ -e "$OUTPUT" ]; then
    echo "ERROR: Output already exists; refusing to overwrite."
    exit 1
fi

echo "=== VERIFY ORIGINAL FILE ==="
gzip -t "$INPUT"

echo "=== FILTER PET >= 5, INTRACHROMOSOMAL ==="

gzip -cd "$INPUT" | awk -F'\t' -v OFS='\t' '
{
    total++

    if (NF != 7 ||
        $2 !~ /^[0-9]+$/ ||
        $3 !~ /^[0-9]+$/ ||
        $5 !~ /^[0-9]+$/ ||
        $6 !~ /^[0-9]+$/ ||
        $7 !~ /^[0-9]+$/ ||
        $2 >= $3 ||
        $5 >= $6) {

        printf "Invalid record at line %d\n", NR > "/dev/stderr"
        exit 2
    }

    if ($1 != $4) {
        trans++
        next
    }

    cis++

    if ($7 + 0 >= 5) {
        print
        retained++
    } else {
        below5++
    }
}
END {
    print "=== INPUT AUDIT ===" > "/dev/stderr"
    printf "Total records: %d\n", total > "/dev/stderr"
    printf "Intrachromosomal: %d\n", cis > "/dev/stderr"
    printf "Interchromosomal: %d\n", trans > "/dev/stderr"
    printf "Cis records with PET < 5: %d\n", below5 > "/dev/stderr"
    printf "Retained cis PET >= 5: %d\n", retained > "/dev/stderr"

    if (total == 0) exit 1
}
' > "$OUTPUT.tmp"

mv "$OUTPUT.tmp" "$OUTPUT"

echo "=== FILTERED FILE ==="
wc -l "$OUTPUT"
ls -lh "$OUTPUT"

echo "=== FIRST FIVE FILTERED LOOPS ==="
head -n 5 "$OUTPUT"

echo "=== RNAPII PREPARATION COMPLETED ==="
