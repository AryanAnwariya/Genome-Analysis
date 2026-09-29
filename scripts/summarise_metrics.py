#!/usr/bin/env python3
"""Summarise *_s<seed>_metrics.txt files in a directory.
Per method: lowest-kmeans-inertia seed is reported; ARI range over seeds is shown."""
import re, sys, glob, os

d = sys.argv[1]
labels = [
    ("LINE25M_scimerge",    "LINE - paper settings (25M, SCI merge)"),
    ("LINE_o1_1000M",       "LINE - order 1 only (1000M)"),
    ("LINE1000M_l2each",    "LINE - 1000M, per-order L2"),
    ("DeepWalk",            "Weighted DeepWalk"),
    ("HOPE",                "HOPE (Katz)"),
    ("DW_HOPE",             "DeepWalk + HOPE (L2 each)"),
]

def f(pat, t):
    m = re.search(pat, t)
    return m.groups() if m else None

data = {}
ref = None
cfg = None
for path in sorted(glob.glob(os.path.join(d, "*_s[0-9]*_metrics.txt"))):
    m = re.match(r"(.+)_s(\d+)_metrics\.txt$", os.path.basename(path))
    if not m:
        continue
    meth, seed = m.group(1), int(m.group(2))
    t = open(path).read()
    inertia = float(f(r"kmeans inertia: ([\d.]+)", t)[0])
    sil = float(f(r"silhouette \(euclidean\): ([\d.-]+)", t)[0])
    db = float(f(r"davies-bouldin: ([\d.-]+)", t)[0])
    sizes = f(r"cluster sizes: \[(.*?)\]", t)[0]
    n, ar, nr = f(r"vs Rao A1-B3 \((\d+) nodes\): ARI=([\d.-]+) NMI=([\d.-]+)", t)
    _, asci, nsci = f(r"vs SCI published \((\d+) nodes\): ARI=([\d.-]+) NMI=([\d.-]+)", t)
    ref = f(r"REFERENCE published SCI vs Rao \(\d+ bins\): ARI=([\d.-]+) NMI=([\d.-]+)", t)
    cfg = f(r"nodes clustered: (\d+), dims: (\d+), k=(\d+), seed=\d+, n_init=(\d+)", t)
    data.setdefault(meth, []).append(dict(seed=seed, inertia=inertia, sil=sil, db=db,
        sizes=sizes, n=n, ar=float(ar), nr=float(nr), asci=float(asci), nsci=float(nsci)))

print("Directory:", d)
print("nodes clustered %s | k=%s | n_init=%s | evaluated on %s Rao A1-B3 bins" % (cfg[0], cfg[2], cfg[3], n))
print("Reported row = lowest-inertia seed per method\n")
h = "%-38s %5s %7s %6s %8s %8s %8s %8s %13s"
print(h % ("Method", "seeds", "Silh", "DB", "ARI(Rao)", "NMI(Rao)", "ARI(SCI)", "NMI(SCI)", "ARI min-max"))
print("-" * 112)
for key, lab in labels:
    if key not in data:
        continue
    rs = data[key]
    b = min(rs, key=lambda r: r["inertia"])
    ars = [r["ar"] for r in rs]
    print(h % (lab, len(rs), "%.3f" % b["sil"], "%.2f" % b["db"], "%.4f" % b["ar"], "%.4f" % b["nr"],
               "%.4f" % b["asci"], "%.4f" % b["nsci"], "%.3f-%.3f" % (min(ars), max(ars))))
print("-" * 112)
print(h % ("Published SCI vs Rao (reference)", "", "", "", ref[0], ref[1], "", "", ""))
print("\nCluster sizes (lowest-inertia seed):")
for key, lab in labels:
    if key in data:
        b = min(data[key], key=lambda r: r["inertia"])
        print("  %-38s seed %d  [%s]" % (lab, b["seed"], b["sizes"]))
