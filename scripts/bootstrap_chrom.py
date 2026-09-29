#!/usr/bin/env python3
"""Chromosome-level bootstrap (resample chromosomes with replacement) for ARI vs Rao and replication-timing eta2."""
import re, glob, os
import numpy as np, pandas as pd
from sklearn.metrics import adjusted_rand_score as ari

B = 1000
rng = np.random.default_rng(2026)
EVAL = ["A1", "A2", "B1", "B2", "B3"]
graphs = [("INTERKR", "results/gm12878_interkr/ninit100"), ("Pair-O/E", "results/gm12878_interkr_oe/eval")]
methods = [("LINE25M_scimerge", "LINE paper"), ("LINE_o1_1000M", "LINE order1 1000M"),
           ("LINE1000M_l2each", "LINE 1000M L2"), ("DeepWalk", "DeepWalk"),
           ("HOPE", "HOPE"), ("DW_HOPE", "DW+HOPE")]
sig = pd.read_csv("results/bio/bin_signals.tsv", sep="\t")[["node", "RepliTiming"]]

def load(path):
    df = pd.read_csv(path, sep="\t", keep_default_na=False)
    df = df[df["rao"].isin(EVAL)].merge(sig, on="node", how="left")
    return df.sort_values("node").reset_index(drop=True)

def best_file(d, key):
    fs = sorted(glob.glob(os.path.join(d, key + "_s[0-9]_metrics.txt")))
    if not fs:
        return None
    b = min(fs, key=lambda p: float(re.search(r"kmeans inertia: ([\d.]+)", open(p).read()).group(1)))
    return b.replace("_metrics.txt", "_clusters.tsv")

def eta2(y, g):
    n = np.bincount(g)
    m = np.bincount(g, weights=y) / np.maximum(n, 1)
    mu = y.mean()
    return (n * (m - mu) ** 2).sum() / ((y - mu) ** 2).sum()

ref = None
L = {}
for gname, d in graphs:
    for key, lab in methods:
        f = best_file(d, key)
        if f is None:
            continue
        df = load(f)
        if ref is None:
            ref = df
        assert (df["node"].values == ref["node"].values).all()
        L[gname + " | " + lab] = pd.factorize(df["cluster"])[0]
rao = pd.factorize(ref["rao"])[0]
sci = pd.factorize(ref["sci_pub"])[0]
y = ref["RepliTiming"].values.astype(float)
chrom = ref["chrom"].values
cl = sorted(set(chrom))
by = {c: np.where(chrom == c)[0] for c in cl}
names = ["ref | Rao", "ref | Published SCI"] + list(L.keys())

def metrics(idx):
    yy, r = y[idx], rao[idx]
    out = {"ref | Rao": (np.nan, eta2(yy, r)), "ref | Published SCI": (ari(r, sci[idx]), eta2(yy, sci[idx]))}
    for k, v in L.items():
        out[k] = (ari(r, v[idx]), eta2(yy, v[idx]))
    return out

full = metrics(np.arange(len(y)))
A = {k: np.zeros(B) for k in names}
E = {k: np.zeros(B) for k in names}
for b in range(B):
    pick = rng.choice(cl, size=len(cl), replace=True)
    m = metrics(np.concatenate([by[c] for c in pick]))
    for k in names:
        A[k][b], E[k][b] = m[k]
print("chromosomes resampled: %d groups, B=%d, bins=%d" % (len(cl), B, len(y)))

def ci(a):
    return np.nanpercentile(a, 2.5), np.nanpercentile(a, 97.5)

print("\n%-30s %8s %17s   %8s %17s" % ("Clustering", "ARI", "95% CI", "eta2 time", "95% CI"))
rows = []
for k in names:
    a, e = full[k]
    la, ha = ci(A[k]) if not np.isnan(a) else (np.nan, np.nan)
    le, he = ci(E[k])
    print("%-30s %8.3f  [%6.3f,%6.3f]   %8.3f  [%6.3f,%6.3f]" % (k, a, la, ha, e, le, he))
    rows.append([k, a, la, ha, e, le, he])
pd.DataFrame(rows, columns=["clustering", "ARI", "ARI_lo", "ARI_hi", "eta2", "eta2_lo", "eta2_hi"]).to_csv("results/bio/bootstrap_table.csv", index=False, float_format="%.4f")

def diff(label, da, de, pa, pe):
    ca = ci(da) if da is not None else None
    ce = ci(de)
    s = "%-46s" % label
    if da is not None:
        s += " dARI %+.3f [%+.3f,%+.3f]%s" % (pa, ca[0], ca[1], "*" if ca[0] > 0 or ca[1] < 0 else " ")
    s += "  d_eta2 %+.3f [%+.3f,%+.3f]%s" % (pe, ce[0], ce[1], "*" if ce[0] > 0 or ce[1] < 0 else " ")
    print(s)

print("\nPAIRED DIFFERENCES (* = 95%% interval excludes 0). Reference: ARI vs published SCI; eta2 vs Rao.")
for k in L:
    diff(k + "  minus reference", A[k] - A["ref | Published SCI"], E[k] - E["ref | Rao"],
         full[k][0] - full["ref | Published SCI"][0], full[k][1] - full["ref | Rao"][1])
print("\nPair-O/E minus INTERKR")
for key, lab in methods:
    a, b2 = "INTERKR | " + lab, "Pair-O/E | " + lab
    if a in L and b2 in L:
        diff(lab, A[b2] - A[a], E[b2] - E[a], full[b2][0] - full[a][0], full[b2][1] - full[a][1])
print("\nDW+HOPE minus HOPE")
for g, _ in graphs:
    a, b2 = g + " | HOPE", g + " | DW+HOPE"
    diff(g, A[b2] - A[a], E[b2] - E[a], full[b2][0] - full[a][0], full[b2][1] - full[a][1])
