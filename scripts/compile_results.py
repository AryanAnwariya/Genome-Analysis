#!/usr/bin/env python3
"""Compile n_init=100 results for both graphs into one CSV + Markdown file.
Lowest-inertia seed per method; chromosome NMI computed from the clusters file."""
import re, glob, os
import pandas as pd
from sklearn.metrics import normalized_mutual_info_score as nmi

graphs = [("INTERKR", "results/gm12878_interkr/ninit100"),
          ("Pair-O/E", "results/gm12878_interkr_oe/eval")]
methods = [("LINE25M_scimerge", "LINE paper settings (25M, SCI merge)"),
           ("LINE_o1_1000M", "LINE order 1 only (1000M)"),
           ("LINE1000M_l2each", "LINE 1000M, per-order L2"),
           ("DeepWalk", "Weighted DeepWalk"),
           ("HOPE", "HOPE (Katz)"),
           ("DW_HOPE", "DeepWalk + HOPE (L2 each)")]
EVAL = ["A1", "A2", "B1", "B2", "B3"]

def g(pat, t):
    return re.search(pat, t).groups()

rows = []
ref = None
for gname, d in graphs:
    for key, lab in methods:
        files = sorted(glob.glob(os.path.join(d, key + "_s[0-9]_metrics.txt")))
        if not files:
            continue
        recs = []
        for p in files:
            t = open(p).read()
            seed = int(re.search(r"_s(\d+)_metrics", p).group(1))
            inertia = float(g(r"kmeans inertia: ([\d.]+)", t)[0])
            sil = float(g(r"silhouette \(euclidean\): ([\d.-]+)", t)[0])
            db = float(g(r"davies-bouldin: ([\d.-]+)", t)[0])
            sizes = [int(x) for x in g(r"cluster sizes: \[(.*?)\]", t)[0].split(",")]
            ar, nr = map(float, g(r"vs Rao A1-B3 \(\d+ nodes\): ARI=([\d.-]+) NMI=([\d.-]+)", t))
            asci, nsci = map(float, g(r"vs SCI published \(\d+ nodes\): ARI=([\d.-]+) NMI=([\d.-]+)", t))
            ref = tuple(map(float, g(r"REFERENCE published SCI vs Rao \(\d+ bins\): ARI=([\d.-]+) NMI=([\d.-]+)", t)))
            cdf = pd.read_csv(p.replace("_metrics.txt", "_clusters.tsv"), sep="\t", keep_default_na=False)
            ev = cdf[cdf["rao"].isin(EVAL)]
            recs.append(dict(graph=gname, method=lab, seed=seed, inertia=inertia, silhouette=sil, DB=db,
                             ARI_Rao=ar, NMI_Rao=nr, ARI_SCI=asci, NMI_SCI=nsci,
                             NMI_chrom=nmi(ev["chrom"], ev["cluster"]),
                             largest_cluster_frac=max(sizes) / sum(sizes)))
        best = min(recs, key=lambda r: r["inertia"])
        best["n_seeds"] = len(recs)
        best["ARI_min"] = min(r["ARI_Rao"] for r in recs)
        best["ARI_max"] = max(r["ARI_Rao"] for r in recs)
        rows.append(best)

df = pd.DataFrame(rows)
df.to_csv("results/compiled/step49_compiled.csv", index=False, float_format="%.4f")

cols = ["method", "seed", "silhouette", "DB", "ARI_Rao", "NMI_Rao", "ARI_SCI", "NMI_SCI", "NMI_chrom", "largest_cluster_frac"]
hdr = ["Method", "Seed", "Silh", "DB", "ARI(Rao)", "NMI(Rao)", "ARI(SCI)", "NMI(SCI)", "NMI(chr)", "Largest cluster"]
out = ["# GM12878 results (n_init=100, lowest-inertia seed of 3, 26,226 Rao A1-B3 bins)", ""]
for gname, _ in graphs:
    out += ["## " + gname + " graph", "", "| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    for _, r in df[df.graph == gname].iterrows():
        vals = [r["method"], str(int(r["seed"])), "%.3f" % r["silhouette"], "%.2f" % r["DB"],
                "%.3f" % r["ARI_Rao"], "%.3f" % r["NMI_Rao"], "%.3f" % r["ARI_SCI"], "%.3f" % r["NMI_SCI"],
                "%.3f" % r["NMI_chrom"], "%.0f%%" % (100 * r["largest_cluster_frac"])]
        out.append("| " + " | ".join(vals) + " |")
    out.append("")
out += ["Reference: published SCI vs Rao ARI %.3f, NMI %.3f. NMI(chr) reference: Rao 0.100, published SCI 0.035." % ref, "",
        "Caveats:",
        "1. Agreement with Rao is not the same as biologically better; biological validation is pending.",
        "2. DW+HOPE and pair-O/E were informed by inspecting earlier results; each was pre-stated and tested once.",
        "3. LINE sample count and merge were tuned after seeing scores; paper-settings LINE is the main LINE result.",
        "4. Own choices: INTERKR normalisation, HOPE beta and dimension, DeepWalk parameters.",
        "5. On INTERKR, HOPE, DW+HOPE and LINE 1000M have chromosome NMI 1.5-2x Rao's 0.100.",
        "6. HOPE and DW+HOPE are roughly tied; no single best method is claimed.",
        "7. LINE paper settings on INTERKR contains a 16-node outlier cluster; LINE order-1-only on pair-O/E has one cluster holding 67% of bins."]
open("results/compiled/step49_compiled.md", "w").write("\n".join(out) + "\n")
print("\n".join(out))
