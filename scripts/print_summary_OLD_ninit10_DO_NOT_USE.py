#!/usr/bin/env python3
"""Print a screenshot-ready comparison table of all GM12878 INTERKR results."""
import re
R = "results/gm12878_interkr"
rows = [
    ("LINE - paper settings (25M, SCI merge)", "eval_LINE"),
    ("LINE - order 1 only (25M)",              "eval_LINE_order1only"),
    ("LINE - order 2 only (25M)",              "eval_LINE_order2only"),
    ("LINE - 1000M, per-order L2 merge",       "eval_LINE1000M_l2each"),
    ("Weighted DeepWalk",                      "eval_DeepWalk"),
    ("HOPE (Katz)",                            "eval_HOPE"),
    ("DeepWalk + HOPE (L2 each, concat)",      "eval_DW_HOPE_l2each"),
]
print("GM12878 Hi-C (Rao 2014), 100 kb, MAPQ>=30, INTERKR-normalised interchromosomal graph")
print("227.5M edges, 26,743 nodes | K-means k=5, seed 12345, n_init=10 | evaluated on 26,226 Rao A1-B3 bins\n")
h = "%-40s %10s %8s %10s %10s %10s %10s"
print(h % ("Method", "Silhouette", "DB", "ARI(Rao)", "NMI(Rao)", "ARI(SCI)", "NMI(SCI)"))
print("-" * 104)
ref = None
for label, f in rows:
    t = open("%s/%s_metrics.txt" % (R, f)).read()
    sil = re.search(r"silhouette \(euclidean\): ([\d.-]+)", t).group(1)
    db = re.search(r"davies-bouldin: ([\d.-]+)", t).group(1)
    ar, nr = re.search(r"vs Rao A1-B3 \(\d+ nodes\): ARI=([\d.-]+) NMI=([\d.-]+)", t).groups()
    asci, nsci = re.search(r"vs SCI published \(\d+ nodes\): ARI=([\d.-]+) NMI=([\d.-]+)", t).groups()
    ref = re.search(r"REFERENCE published SCI vs Rao \(\d+ bins\): ARI=([\d.-]+) NMI=([\d.-]+)", t).groups()
    print(h % (label, "%.3f" % float(sil), "%.2f" % float(db), ar, nr, asci, nsci))
print("-" * 104)
print(h % ("Published SCI labels (reference)", "", "", ref[0], ref[1], "1", "1"))
