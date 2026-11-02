# Dataset Provenance and Processing

## Hi-C contact matrices

Source:

GEO accession GSE63525.

Dataset:

GM12878 combined interchromosomal contact matrices.

The original archive is:

GSE63525_GM12878_combined_interchromosomal_contact_matrices.tar.gz

The complete archive inventory is stored at:

docs/full_archive_inventory_429072.txt

The archive contains nine genomic resolutions:

1 kb, 5 kb, 10 kb, 25 kb, 50 kb, 100 kb,
250 kb, 500 kb, and 1 Mb.

It also includes MAPQG0 and MAPQGE30 categories.

Total archive size:

- Compressed: 35.01 GiB.
- Uncompressed: 209.03 GiB.
- Total files: 52,119.

The graph-embedding experiments documented here
principally use the selected 100 kb MAPQGE30
interchromosomal data.

The entire archive should not be treated as one graph:
different resolutions and mapping-quality categories
represent alternative views of overlapping genomic data.

## CTCF ChIA-PET

Four-Dimensional Nucleome dataset identifiers:

- 4DNFIA3H6KXY
- 4DNFIS9CCN6R

Processing included PET-support filtering, coordinate
conversion, genomic-bin assignment, and validation.

The preparation used an hg38-to-hg19 liftOver chain.

IMPORTANT:

The source genome assembly and exact metadata for
these 4DN files should be verified against the original
dataset records before publication.

The coordinate conversion should not be treated as
independently validated source metadata.

## RNAPII ChIA-PET

GEO sample:

GSM1872887

Cell line:

GM12878

The RNAPII interaction data were processed with a
minimum PET support of five.

Intrachromosomal interactions were mapped to common
100 kb genomic bins.

This is an independent dataset and is not claimed to
be the identical RNAPII dataset used in SCI's original
published biological comparison.

## Evaluation references

The comparison includes:

- Rao reference compartment assignments.
- Published SCI compartment labels.
- Predictions from our HOPE, DeepWalk, LINE, and
  combined-embedding experiments.

For our biological comparisons, the common reference
set contained 26,226 genomic bins.

## Storage

Original datasets, large contact graphs, and generated
embedding files are intentionally excluded from Git.

The HPC scratch filesystem has a 30-day file-retention
policy. Important raw data and derived artifacts must
therefore be backed up separately.

The repository stores code, documentation, small
experimental results, and clustering assignments.

## Scientific limitations

The distance-matched validation implemented here is an
additional evaluation procedure.

Its numerical values must not be directly compared with
the raw loop-ratio percentages reported in the SCI paper.

The 4DN data provenance requires final verification before
using these results in a publication.
