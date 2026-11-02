#!/bin/bash
#SBATCH --job-name=ctcf_pet_audit
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/ctcf_pet_audit_%j.log

set -euo pipefail

BASE=/scratch/adas/genome-analysis/data/chia_pet

for ID in 4DNFIA3H6KXY 4DNFIS9CCN6R; do

    echo
    echo "==================================="
    echo "FILE: $ID"
    echo "==================================="

    gzip -cd "$BASE/${ID}.bedpe.gz" | \
    awk -F'\t' '
    {
        total++

        if (NF != 7 || $7 !~ /^[0-9]+$/) {
            invalid++
            next
        }

        value = $7 + 0

        if (value == 0) zero++
        else if (value == 1) one++
        else if (value == 2) two++
        else if (value == 3) three++
        else if (value == 4) four++
        else if (value < 10) five_nine++
        else if (value < 50) ten_fortynine++
        else fifty_plus++

        if (value > max) max = value

        if (value >= 5) pet5++
    }

    END {
        printf "Total records: %d\n", total
        printf "Invalid records: %d\n", invalid
        printf "Column 7 = 0: %d\n", zero
        printf "Column 7 = 1: %d\n", one
        printf "Column 7 = 2: %d\n", two
        printf "Column 7 = 3: %d\n", three
        printf "Column 7 = 4: %d\n", four
        printf "Column 7 = 5-9: %d\n", five_nine
        printf "Column 7 = 10-49: %d\n", ten_fortynine
        printf "Column 7 >= 50: %d\n", fifty_plus
        printf "Column 7 maximum: %d\n", max
        printf "Column 7 >= 5: %d\n", pet5

        if (total > 0) {
            printf "Percentage >= 5: %.4f%%\n", \
                   100 * pet5 / total
        }
    }'

done

echo
echo "=== PET COUNT AUDIT COMPLETED ==="
