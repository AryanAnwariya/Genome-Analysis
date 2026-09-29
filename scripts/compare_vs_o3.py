#!/usr/bin/env python3
"""Paired chromosome-bootstrap differences with LINE order-3 (25M) as baseline. Also chromosome NMI and sizes for order 3."""
import re, glob, os
import numpy as np, pandas as pd
from sklearn.metrics import adjusted_rand_score as sk_ari, normalized_mutual_info_score as sk_nmi

B = 1000
rng = np.random.default_rng(2028)
EVAL = ["A1", "A2", "B1", "B2", "B3"]
graphs = [("INTERKR", "results/gm12878_interkr/ninit100", "results/gm12878_interkr/order3"),
          ("Pair-O/E", "results/gm12878_interkr_oe/eval", "results/gm12878_interkr_oe/order3")]
keys = {"LINE25M_scimerge": "LINE paper", "LINE1000M_l2each": "LINE 1000M L2",
        "DeepWalk": "DeepWalk", "HOPE": "HOPE", "DW_HOPE": "DW+HOPE"}
sig = pd.read_csv("results/bio/bin_signals.tsv", sep="\t")[["node", "RepliTiming"]]

def best_file(d, key):
    fs = sorted(glob.glob(os.path.join(d, key + "_s[0-9]_metrics.txt")))
    assert fs, (d, key)
    b = min(fs, key=lambda p: float(re.search(r"kmeans inertia: ([\d.]+)", open(p).read()).group(1)))
    return b.replace("_metrics.txt", "_clusters.tsv")

def load(path):
    df = pd.read_csv(path, sep="\t", keep_default_na=False)
    df = df[df["rao"].isin(EVAL)].merge(sig, on="node", how="left")
    return df.sort_values("node").reset_index(drop=True)

def ari_c(c):
    c = c.astype(float); n = c.sum()
    sc = (c * (c - 1) / 2).sum(); a = c.sum(1); b = c.sum(0)
    sa = (a * (a - 1) / 2).sum(); sb = (b * (b - 1) / 2).sum()
    ex = sa * sb / (n * (n - 1) / 2)
    return (sc - ex) / ((sa + sb) / 2 - ex)

def nmi_c(c):
    c = c.astype(float); n = c.sum(); p = c / n
    pa = p.sum(1, keepdims=True); pb = p.sum(0, keepdims=True)
    nz = p > 0
    mi = (p[nz] * np.log(p[nz] / (pa @ pb)[nz])).sum()
    ha = -(pa[pa > 0] * np.log(pa[pa > 0])).sum(); hb = -(pb[pb > 0] * np.log(pb[pb > 0])).sum()
    return mi / ((ha + hb) / 2)

def eta2(y, g, k=5):
    n = np.bincount(g, minlength=k); m = np.bincount(g, weights=y, minlength=k) / np.maximum(n, 1)
    mu = y.mean()
    return (n * (m - mu) ** 2).sum() / ((y - mu) ** 2).sum()

ref = None
L = {}
raw = {}
for gname, d, d3 in graphs:
    for k, lab in keys.items():
        df = load(best_file(d, k))
        ref = df if ref is None else ref
        assert (df["node"].values == ref["node"].values).all()
        L[(gname, lab)] = pd.factorize(df["cluster"])[0]
    df3 = load(best_file(d3, "LINE_o3_25M"))
    assert (df3["node"].values == ref["node"].values).all()
    L[(gname, "LINE o3 25M")] = pd.factorize(df3["cluster"])[0]
    raw[gname] = df3
rao = pd.factorize(ref["rao"])[0]
y = ref["RepliTiming"].values.astype(float)
chrom = ref["chrom"].values
cl = sorted(set(chrom))
by = {c: np.where(chrom == c)[0] for c in cl}

def met(idx, lab):
    c = np.bincount(rao[idx] * 5 + lab[idx], minlength=25).reshape(5, 5)
    return np.array([ari_c(c), nmi_c(c), eta2(y[idx], lab[idx])])

allidx = np.arange(len(y))
for kk, lab in L.items():
    m = met(allidx, lab)
    assert abs(m[0] - sk_ari(rao, lab)) < 1e-9 and abs(m[1] - sk_nmi(rao, lab)) < 1e-9, kk
print("fast ARI/NMI match scikit-learn on all 10 clusterings")
for g in ["INTERKR", "Pair-O/E"]:
    f = met(allidx, L[(g, "LINE o3 25M")])
    print("%s LINE order 3: ARI %.4f NMI %.4f timing eta2 %.3f NMI(chr) %.4f sizes %s" % (
        g, f[0], f[1], f[2], sk_nmi(raw[g]["chrom"], raw[g]["cluster"]), sorted(np.bincount(L[(g, 'LINE o3 25M')]).tolist(), reverse=True)))

full = {k: met(allidx, v) for k, v in L.items()}
boots = {k: np.zeros((B, 3)) for k in L}
for b in range(B):
    idx = np.concatenate([by[c] for c in rng.choice(cl, size=len(cl), replace=True)])
    for k, v in L.items():
        boots[k][b] = met(idx, v)

pairs = []
for g in ["INTERKR", "Pair-O/E"]:
    for m in ["DeepWalk", "HOPE", "DW+HOPE", "LINE 1000M L2"]:
        pairs.append(((g, m), (g, "LINE o3 25M")))
    pairs.append(((g, "LINE o3 25M"), (g, "LINE paper")))
print("\nPAIRED DIFFERENCES = first minus second; 95 pct CI; * = interval excludes 0")
print("%-9s %-14s minus %-14s | %-24s | %-24s | %-24s" % ("graph", "method", "baseline", "dARI", "dNMI", "d timing eta2"))
for a, b_ in pairs:
    d = full[a] - full[b_]
    bd = boots[a] - boots[b_]
    lo = np.percentile(bd, 2.5, axis=0); hi = np.percentile(bd, 97.5, axis=0)
    cells = ["%+.3f [%+.3f,%+.3f]%s" % (d[i], lo[i], hi[i], "*" if (lo[i] > 0 or hi[i] < 0) else " ") for i in range(3)]
    print("%-9s %-14s minus %-14s | %s | %s | %s" % (a[0], a[1], b_[1], cells[0], cells[1], cells[2]))
