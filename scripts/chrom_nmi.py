#!/usr/bin/env python3
"""NMI(cluster, chromosome) on Rao A1-B3 bins, for every *_s<seed>_clusters.tsv in a directory."""
import sys, glob, os, re
import pandas as pd
from sklearn.metrics import normalized_mutual_info_score as nmi

d = sys.argv[1]
labels = [
    ("LINE25M_scimerge",  "LINE - paper settings (25M, SCI merge)"),
    ("LINE_o1_1000M",     "LINE - order 1 only (1000M)"),
    ("LINE1000M_l2each",  "LINE - 1000M, per-order L2"),
    ("DeepWalk",          "Weighted DeepWalk"),
    ("HOPE",              "HOPE (Katz)"),
    ("DW_HOPE",           "DeepWalk + HOPE (L2 each)"),
]
EVAL = ["A1", "A2", "B1", "B2", "B3"]
print("Directory:", d)
print("NMI(label, chromosome) on Rao A1-B3 bins; lowest-inertia seed per method\n")
h = "%-38s %6s %10s %14s"
print(h % ("Method", "seed", "NMI(chr)", "min-max seeds"))
print("-" * 72)
ref_df = None
for key, lab in labels:
    files = sorted(glob.glob(os.path.join(d, key + "_s*_clusters.tsv")))
    if not files:
        continue
    res = []
    for p in files:
        df = pd.read_csv(p, sep="\t", keep_default_na=False)
        ev = df[df["rao"].isin(EVAL)]
        ref_df = ev
        mt = open(p.replace("_clusters.tsv", "_metrics.txt")).read()
        inertia = float(re.search(r"kmeans inertia: ([\d.]+)", mt).group(1))
        seed = int(re.search(r"_s(\d+)_clusters", p).group(1))
        res.append((inertia, seed, nmi(ev["chrom"], ev["cluster"]), len(ev)))
    best = min(res)
    vals = [r[2] for r in res]
    print(h % (lab, best[1], "%.4f" % best[2], "%.3f-%.3f" % (min(vals), max(vals))))
print("-" * 72)
print(h % ("Rao labels (reference)", "", "%.4f" % nmi(ref_df["chrom"], ref_df["rao"]), ""))
sci = ref_df[ref_df["sci_pub"] != "none"]
print(h % ("Published SCI labels (reference)", "", "%.4f" % nmi(sci["chrom"], sci["sci_pub"]), "n=%d" % len(sci)))
print("\nbins used: %d" % len(ref_df))
