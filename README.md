# Graph embeddings for Hi-C sub-compartments (SCI reproduction and comparison)

Mid-semester project. We reproduce **SCI** (Ashoor et al., *Nature Communications* 2020,
doi:10.1038/s41467-020-14974-x, code: https://github.com/TheJacksonLaboratory/sci) on real
GM12878 Hi-C data (Rao et al. 2014, GEO GSE63525) and compare LINE with weighted DeepWalk,
HOPE (Katz) and DeepWalk+HOPE. Next part of the project: a data structure and query
processor for genome contact data.

## Main findings
- SCI's way of joining the two LINE orders (train separately, concatenate, then L2) collapses:
  order-2 vectors are ~70x longer than order-1. The authors' own container shows the same collapse.
- HOPE and DeepWalk+HOPE give the best agreement with Rao's sub-compartments. A per-chromosome-pair
  observed/expected correction (pair-O/E) helps every well-formed clustering.
  Best: DeepWalk+HOPE on pair-O/E, ARI 0.633 vs Rao (published SCI labels: 0.503).
- The improvement is modest: ~0.06-0.11 ARI over an improved LINE; on replication timing the best
  methods match Rao's labels but do not beat them.
- We could not reproduce the published SCI labels. The authors' container does not reproduce its
  own saved demo result and gives different results on repeated runs.

## Layout
- `scripts/`  graph building, embeddings (LINE runs, weighted DeepWalk, HOPE), evaluation, bootstrap,
  reproduction checks, Slurm job scripts
- `results/`  metrics files, compiled tables, replication-timing, bootstrap and container results
- `docs/`     `DECISIONS.txt` (dated decision log), `LINE_BUILD_NOTES.txt`, conda env, extraction list
- `report/`   mid-semester report (LaTeX + figures)

## Not in this repo (too large or third-party)
| Item | Size | How to get it |
|---|---|---|
| GSE63525 GM12878 interchromosomal archive | 36 GB (224 GB uncompressed) | GEO GSE63525 |
| Extracted 100 kb autosomal files (MAPQ>=30) | 5.82 GB | `tar -xzf ... -T docs/extract_list_100kb_autosomes.txt` |
| Rao sub-compartment labels | small | GEO GSE63525 (`GSE63525_GM12878_subcompartments.bed.gz`) |
| Published SCI labels | small | SCI repo, `predictions/` |
| Graphs (edge lists) | 8.6-11 GB each | `scripts/build_sci_graph.py` |
| Embeddings | 25 MB each | `scripts/run_line*.sh`, `run_deepwalk.sh`, `run_hope.sh` |
| SCI authors' container | 3.0 GB | `apptainer pull docker://yuz12012/sci_container:latest` |
| ENCODE bigWigs (hg19) | 10-663 MB | UCSC hg19 encodeDCC (see `docs/DECISIONS.txt`) |
| LINE binary | - | build SCI's `LINE/line.cpp` (see `docs/LINE_BUILD_NOTES.txt`) |

## Environment
BITS Hyderabad HPC (Slurm). Python 3.10 env in `docs/genome_env.yml`.
Scripts contain absolute paths for our cluster (`/scratch/adas/genome-analysis`); change them for other systems.

---

## Biological Validation — CTCF and RNAPII

We extended the graph-embedding experiments with
independent biological validation using CTCF and RNAPII
ChIA-PET interaction datasets.

The analysis compares HOPE, DeepWalk, LINE, and combined
DeepWalk + HOPE predictions with Rao and published SCI
sub-compartment labels.

The main RNAPII experiment showed stronger
distance-matched enrichment for several embedding
approaches, particularly LINE. However, the results
were sensitive to genomic separation, and the CTCF
analysis did not demonstrate a clear advantage over Rao.

**These findings do not establish overall superiority
over the published SCI method.**

Detailed documentation:

- [Biological validation](docs/BIOLOGICAL_VALIDATION.md)
- [Data provenance](docs/DATA_PROVENANCE.md)
- [Full Hi-C archive inventory](docs/full_archive_inventory_429072.txt)

Code:

- `scripts/biological_validation/`

Results:

- `results/biological_validation/`
- `results/gm12878_interkr_oe/eval/`

The repository intentionally excludes large raw datasets
and generated graph/embedding files.
