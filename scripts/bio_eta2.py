#!/usr/bin/env python3
"""Eta-squared of bin-level signal explained by cluster labels (5 clusters) on Rao A1-B3 bins."""
import os, glob, re
import numpy as np, pandas as pd, pyBigWig

BIO = "data/bio"
TRACKS = [("RepliTiming", "wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig"),
          ("H3K27me3", "wgEncodeBroadHistoneGm12878H3k27me3StdSigV2.bigWig"),
          ("H3K9me3", "wgEncodeBroadHistoneGm12878H3k9me3StdSig.bigWig"),
          ("H3K4me3", "wgEncodeBroadHistoneGm12878H3k4me3StdSig.bigWig")]
names = [n for n, _ in TRACKS]
graphs = [("INTERKR", "results/gm12878_interkr/ninit100"),
          ("Pair-O/E", "results/gm12878_interkr_oe/eval")]
methods = [("LINE25M_scimerge", "LINE paper settings"),
           ("LINE_o1_1000M", "LINE order 1 only 1000M"),
           ("LINE1000M_l2each", "LINE 1000M per-order L2"),
           ("DeepWalk", "Weighted DeepWalk"),
           ("HOPE", "HOPE (Katz)"),
           ("DW_HOPE", "DeepWalk + HOPE")]
EVAL = ["A1", "A2", "B1", "B2", "B3"]

base = pd.read_csv("results/gm12878_interkr/ninit100/DW_HOPE_s1_clusters.tsv", sep="\t", keep_default_na=False)
sig_path = "results/bio/bin_signals.tsv"
if os.path.exists(sig_path):
    sig = pd.read_csv(sig_path, sep="\t")
else:
    sig = base[["node", "chrom", "start", "end"]].copy()
    for name, fn in TRACKS:
        bw = pyBigWig.open(os.path.join(BIO, fn))
        L = bw.chroms()
        vals = []
        for c, s, e in zip(sig.chrom, sig.start, sig.end):
            e2 = min(int(e), L.get(c, 0))
            if c not in L or int(s) >= e2:
                vals.append(np.nan)
                continue
            v = bw.stats(c, int(s), e2, type="mean")[0]
            vals.append(np.nan if v is None else v)
        sig[name] = vals
        bw.close()
        print("read track", name, flush=True)
    sig.to_csv(sig_path, sep="\t", index=False)

def eta2(y, g):
    ok = ~np.isnan(y)
    y, g = y[ok], g[ok]
    _, inv = np.unique(g, return_inverse=True)
    n = np.bincount(inv)
    m = np.bincount(inv, weights=y) / n
    mu = y.mean()
    return (n * (m - mu) ** 2).sum() / ((y - mu) ** 2).sum()

def eta2_wc(y, g, chrom):
    ok = ~np.isnan(y)
    y, g, chrom = y[ok], g[ok], chrom[ok]
    r = y - pd.Series(y).groupby(chrom).transform("mean").values
    return eta2(r, g)

def load(path):
    df = pd.read_csv(path, sep="\t", keep_default_na=False)
    df = df[df["rao"].isin(EVAL)].merge(sig[["node"] + names], on="node", how="left")
    return df

rows = []
def add(graph, label, df, col):
    r = {"graph": graph, "clustering": label}
    for nm in names:
        y = df[nm].values.astype(float)
        r[nm] = eta2(y, df[col].values)
        r[nm + "_wc"] = eta2_wc(y, df[col].values, df["chrom"].values)
    rows.append(r)

ref = load("results/gm12878_interkr/ninit100/DW_HOPE_s1_clusters.tsv")
print("evaluation bins: %d; non-missing signal per track: %s" %
      (len(ref), ", ".join("%s=%d" % (n, ref[n].notna().sum()) for n in names)))
add("reference", "Rao labels", ref, "rao")
add("reference", "Published SCI labels", ref, "sci_pub")
add("reference", "Chromosome as label (22 groups)", ref, "chrom")
rng = np.random.default_rng(0)
sh = []
for _ in range(100):
    tmp = ref.copy()
    tmp["shuf"] = rng.permutation(ref["rao"].values)
    sh.append([eta2(tmp[nm].values.astype(float), tmp["shuf"].values) for nm in names])
print("shuffled-label sanity (mean raw eta2 over 100): " + ", ".join("%s=%.5f" % (n, v) for n, v in zip(names, np.mean(sh, axis=0))))

for gname, d in graphs:
    for key, lab in methods:
        files = sorted(glob.glob(os.path.join(d, key + "_s[0-9]_metrics.txt")))
        if not files:
            continue
        best = min(files, key=lambda p: float(re.search(r"kmeans inertia: ([\d.]+)", open(p).read()).group(1)))
        add(gname, lab, load(best.replace("_metrics.txt", "_clusters.tsv")), "cluster")

res = pd.DataFrame(rows)
res.to_csv("results/bio/eta2_table.csv", index=False, float_format="%.4f")
fmt = lambda v: "%.3f" % v
print("\nRAW eta-squared (share of signal variance explained by the labels)")
print(res[["graph", "clustering"] + names].to_string(index=False, float_format=fmt))
print("\nWITHIN-CHROMOSOME eta-squared (signal centred per chromosome first)")
print(res[["graph", "clustering"] + [n + "_wc" for n in names]].rename(columns=lambda c: c.replace("_wc", "")).to_string(index=False, float_format=fmt))
