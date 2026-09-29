#!/usr/bin/env python3
"""
Build an SCI-compatible interchromosomal Hi-C graph from Rao et al. 2014
(GSE63525) sparse RAWobserved files, optionally normalised with a Rao norm vector.

Reproduces the rules of sci/hic.py (TheJacksonLaboratory/sci):
  - bins per chromosome = size // res + 1; contacts beyond that are dropped
  - duplicate cells are summed
  - only interchromosomal cells are used
  - an edge is written only if value >= 1, weight = ln(value),
    in both directions, formatted "%d\t%d\t%f"
Node IDs follow the chromosome order of the chrom sizes file (same layout as
SCI's GW_meta_data). The full node table is written to <prefix>_nodes.tsv.
Normalised value = RAW / (v_i * v_j)  (Rao 2014 README).
"""
import argparse
import gzip
import os
from itertools import combinations

import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--rao-dir", required=True,
               help="path to .../100kb_resolution_interchromosomal")
p.add_argument("--chrom-sizes", required=True)
p.add_argument("--res", type=int, default=100000)
p.add_argument("--norm", default="INTERKRnorm",
               help="RAW, or a norm-vector suffix such as INTERKRnorm")
p.add_argument("--mapq", default="MAPQGE30")
p.add_argument("--out-prefix", required=True)
p.add_argument("--write-7col", action="store_true",
               help="also write SCI 7-column input (gzipped) for auditing")
p.add_argument("--test", action="store_true",
               help="process only the chr1-chr2 pair")
p.add_argument("--pair-oe", action="store_true",
               help="divide each chromosome-pair block by its mean over valid bins (observed/expected)")
a = p.parse_args()

res = a.res
tag = "%dkb" % (res // 1000)
THRESHOLDS = [0.5, 1, 2, 5, 10]


def num(c):
    return int(c.replace("chr", ""))


def lookup(vec, idx):
    out = np.full(idx.shape, np.nan)
    ok = idx < len(vec)
    out[ok] = vec[idx[ok]]
    return out


# ---- node table, same layout as SCI GW_meta_data ----
chroms, nbins, offset = [], {}, {}
with open(a.chrom_sizes) as f:
    for line in f:
        c, size = line.split()[:2]
        chroms.append(c)
        nbins[c] = int(size) // res + 1
n_nodes = 0
with open(a.out_prefix + "_nodes.tsv", "w") as f:
    f.write("node\tchrom\tstart\tend\n")
    for c in chroms:
        offset[c] = n_nodes
        for b in range(nbins[c]):
            f.write("%d\t%s\t%d\t%d\n" % (n_nodes, c, b * res, b * res + res))
            n_nodes += 1
print("Chromosomes: %d, total possible nodes: %d" % (len(chroms), n_nodes))

pairs = list(combinations(chroms, 2))
if a.test:
    pairs = [("chr1", "chr2")]

stats = dict(read=0, out_of_range=0, nonfinite_or_zero=0, duplicates=0, cells=0)
ge = {t: 0 for t in THRESHOLDS}
used = np.zeros(n_nodes, dtype=bool)
gfile = open(a.out_prefix + "_HiC_graph.txt", "w")
sfile = gzip.open(a.out_prefix + "_sci_input.txt.gz", "wt") if a.write_7col else None

for k, (c1, c2) in enumerate(pairs, 1):
    A, B = (c1, c2) if num(c1) < num(c2) else (c2, c1)
    d = os.path.join(a.rao_dir, "%s_%s" % (A, B), a.mapq)
    df = pd.read_csv(os.path.join(d, "%s_%s_%s.RAWobserved" % (A, B.replace("chr", ""), tag)),
                     sep=r"\s+", header=None, names=["i", "j", "m"])
    stats["read"] += len(df)
    bi = df["i"].values.astype(np.int64) // res
    bj = df["j"].values.astype(np.int64) // res
    val = df["m"].values.astype(np.float64)

    ok = (bi < nbins[A]) & (bj < nbins[B])
    stats["out_of_range"] += int((~ok).sum())
    bi, bj, val = bi[ok], bj[ok], val[ok]

    if a.norm != "RAW":
        vA = np.loadtxt(os.path.join(d, "%s_%s.%s" % (A, tag, a.norm)))
        vB = np.loadtxt(os.path.join(d, "%s_%s.%s" % (B, tag, a.norm)))
        with np.errstate(divide="ignore", invalid="ignore"):
            val = val / (lookup(vA, bi) * lookup(vB, bj))

    good = np.isfinite(val) & (val > 0)
    stats["nonfinite_or_zero"] += int((~good).sum())
    bi, bj, val = bi[good], bj[good], val[good]
    if a.pair_oe:
        if a.norm != "RAW":
            nA = int(np.sum(np.isfinite(vA) & (vA > 0)))
            nB = int(np.sum(np.isfinite(vB) & (vB > 0)))
        else:
            nA, nB = nbins[A], nbins[B]
        val = val / (val.sum() / (nA * nB))

    # SCI accumulates repeated cells with +=
    s = pd.Series(val).groupby(bi * nbins[B] + bj).sum()
    stats["duplicates"] += len(val) - len(s)
    key = s.index.values
    val = s.values
    bi, bj = key // nbins[B], key % nbins[B]
    stats["cells"] += len(val)
    for t in THRESHOLDS:
        ge[t] += int((val >= t).sum())

    if sfile is not None:
        pd.DataFrame({"c1": A, "s1": bi * res, "e1": bi * res + res,
                      "c2": B, "s2": bj * res, "e2": bj * res + res,
                      "v": val}).to_csv(sfile, sep="\t", header=False,
                                        index=False, float_format="%.6f")

    keep = val >= 1
    n1 = offset[A] + bi[keep]
    n2 = offset[B] + bj[keep]
    w = np.log(val[keep])
    used[n1] = True
    used[n2] = True
    m = len(w)
    src = np.empty(2 * m, dtype=np.int64)
    dst = np.empty(2 * m, dtype=np.int64)
    src[0::2], dst[0::2] = n1, n2
    src[1::2], dst[1::2] = n2, n1
    pd.DataFrame({"a": src, "b": dst, "w": np.repeat(w, 2)}).to_csv(
        gfile, sep="\t", header=False, index=False, float_format="%.6f")
    print("[%3d/%d] %s_%s  entries=%d  kept_edges=%d"
          % (k, len(pairs), A, B, len(df), m), flush=True)

gfile.close()
if sfile is not None:
    sfile.close()

lines = [
    "norm: %s" % a.norm,
    "pairs processed: %d" % len(pairs),
    "RAWobserved entries read: %d" % stats["read"],
    "dropped (beyond chrom length): %d" % stats["out_of_range"],
    "dropped (NaN/inf/zero after norm): %d" % stats["nonfinite_or_zero"],
    "duplicate cells summed: %d" % stats["duplicates"],
    "cells with value > 0: %d" % stats["cells"],
] + ["cells with value >= %s: %d" % (t, ge[t]) for t in THRESHOLDS] + [
    "undirected edges (value >= 1): %d" % ge[1],
    "directed edges written: %d" % (2 * ge[1]),
    "nodes with >= 1 edge: %d of %d" % (used.sum(), n_nodes),
]
with open(a.out_prefix + "_graph_stats.txt", "w") as f:
    f.write("\n".join(lines) + "\n")
print("\n".join(lines))
print("Done.")
