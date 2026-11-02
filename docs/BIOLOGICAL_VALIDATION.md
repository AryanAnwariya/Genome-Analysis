# Biological Validation of Hi-C Graph Embeddings

## Objective

Evaluate whether genomic sub-compartments obtained using
graph embedding methods reflect biologically meaningful
chromatin organization.

The evaluated approaches include:

- HOPE
- DeepWalk
- LINE
- Combined DeepWalk and HOPE embeddings

The analysis compares our predicted clusters against
Rao's sub-compartment labels and published SCI labels.

## Experimental approach

The Hi-C graph experiments used GM12878 interchromosomal
contact data at 100 kb resolution, with MAPQ >= 30.

Embedding vectors were clustered into five groups using
K-means.

Biological validation was conducted using independent
CTCF and RNA Polymerase II (RNAPII) ChIA-PET interactions.

ChIA-PET identifies genomic interactions associated with
specific proteins.

For each interaction, the two genomic anchors were mapped
to 100 kb genomic bins. We then examined whether the
anchors belonged to the same predicted compartment.

## Input filtering

The validation used:

- Intrachromosomal interactions.
- At least five PETs supporting each retained interaction.
- Genomic anchors mapped using a greater-than-50% overlap
  criterion.
- The common set of 26,226 reference genomic bins.

Interactions with both anchors in the same 100 kb bin
were excluded from the main distance-matched comparisons.

## Why distance matching was necessary

Nearby genomic bins frequently share compartment labels.

Therefore, a high fraction of within-compartment
interactions does not necessarily indicate a biologically
superior prediction.

We compared observed interactions with chromosome- and
distance-matched genomic-region pairs.

The resulting enrichment estimates measure the association
between experimentally observed contacts and predicted
compartment membership after accounting for genomic
separation.

These estimates are not directly interchangeable with
the raw intra/inter ratios reported in the SCI paper.

## CTCF validation

Two CTCF ChIA-PET datasets were evaluated.

The larger dataset yielded:

- 67,264 eligible interactions.
- 51,669 interactions between distinct 100 kb bins.

Distance-matched enrichment:

| Method | Enrichment |
|---|---:|
| Rao | 1.1339 |
| Published SCI labels | 1.1174 |
| HOPE | 1.1387 |
| DeepWalk | 1.1768 |
| DeepWalk + HOPE | 1.1530 |
| LINE | 1.1727 |

All tested methods showed modest enrichment.

However, chromosome-bootstrap comparisons did not
establish a statistically clear advantage over Rao.

The smaller CTCF dataset contained only 343 distinct-bin
interactions and was inconclusive.

## RNAPII validation

An independent GM12878 RNAPII ChIA-PET dataset was used.

After filtering:

- 52,848 interactions were eligible.
- 19,959 connected distinct 100 kb bins.

Distance-matched enrichment:

| Method | Enrichment |
|---|---:|
| Rao | 1.5702 |
| Published SCI labels | 1.2172 |
| HOPE | 1.9970 |
| DeepWalk | 1.7862 |
| DeepWalk + HOPE | 2.1144 |
| LINE | 2.1523 |

Relative to Rao, the point estimates were approximately:

- HOPE: 27% higher.
- DeepWalk: 14% higher.
- DeepWalk + HOPE: 35% higher.
- LINE: 37% higher.

These are relative differences in distance-matched
enrichment, not improvements in prediction accuracy.

## Statistical validation

Chromosome-based bootstrap analysis was conducted to
estimate uncertainty.

Multiple-comparison adjustment was applied when testing
several embedding methods against Rao.

HOPE, DeepWalk + HOPE, and LINE showed evidence of an
advantage in the initial RNAPII comparison.

However, approximately 78% of the distinct-bin RNAPII
interactions involved immediately neighboring bins.

Additional analyses excluded short-range interactions.

When adjacent-bin interactions were excluded, LINE
retained evidence of an advantage over Rao.

Under a stricter restriction excluding the two closest
bin separations, none of the embedding methods showed
a statistically established advantage.

These sensitivity analyses limit the strength of the
biological conclusions.

## Interpretation

The experiments demonstrate that graph embedding methods
can identify compartment assignments associated with
independent chromatin-interaction data.

The strongest positive finding was observed for LINE in
the RNAPII validation.

However:

1. CTCF did not establish an improvement over Rao.
2. RNAPII improvements depended partly on genomic distance.
3. The validation datasets and procedures were not identical
   to those used in the original SCI publication.
4. Overall biological superiority over SCI has not been
   demonstrated.
5. Association with ChIA-PET interactions is only one
   aspect of biological validity.

This is a comparative computational investigation with
mixed but informative biological findings.

It should not be represented as a complete reproduction
of SCI's published biological-validation results.

## Reproducibility

Source code:

`scripts/biological_validation/`

Slurm execution scripts:

`scripts/biological_validation/jobs/`

Numerical results:

`results/biological_validation/`

Embedding cluster assignments and clustering metrics:

`results/gm12878_interkr_oe/eval/`

The Slurm scripts contain HPC-specific paths and resource
settings that may require modification on other systems.

Raw genomic datasets and large embedding files are not
stored in GitHub.
