#!/usr/bin/env python3
"""Cluster only the bins the published SCI labels cover; score vs published SCI and Rao."""
import os
import numpy as np, pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score as ari, normalized_mutual_info_score as nmi

base = pd.read_csv("results/gm12878_interkr/ninit100/DW_HOPE_s1_clusters.tsv", sep="\t", keep_default_na=False)
ref = base[base["sci_pub"] != "none"].copy()
ids = ref["node"].astype(str).values
sci = ref["sci_pub"].values
rao = ref["rao"].values
print("bins with a published SCI label:", len(ref), "| Rao labels among them:", sorted(set(rao)))

def read_emb(p):
    df = pd.read_csv(p, sep=r"\s+", skiprows=1, header=None, index_col=0)
    df.index = df.index.astype(str)
    return df.loc[ids].values.astype(float)

def l2(x):
    return x / np.maximum(np.linalg.norm(x, axis=1, keepdims=True), 1e-12)

I = "results/gm12878_interkr/gm12878_INTERKR_order_%s_samples_25M.embedding"
variants = [
    ("INTERKR order 3, 25M", [I % "3"], "concat_l2"),
    ("INTERKR orders 1+2 concatenated (SCI merge)", [I % "1", I % "2"], "concat_l2"),
    ("INTERKR orders 1+2, per-order L2", [I % "1", I % "2"], "l2_each"),
    ("Pair-O/E order 3, 25M", ["results/gm12878_interkr_oe/gm12878_INTERKR_OE_order_3_samples_25M.embedding"], "concat_l2"),
]
print("\n%-46s %5s %8s %8s %8s %8s  %-22s %s" % ("Variant", "seed", "ARI_SCI", "NMI_SCI", "ARI_Rao", "NMI_Rao", "ARI_SCI range (seeds)", "sizes"))
for name, files, mode in variants:
    if not all(os.path.exists(f) for f in files):
        print("%-46s SKIPPED (missing file)" % name)
        continue
    mats = [read_emb(f) for f in files]
    X = np.hstack([l2(m) for m in mats]) if mode == "l2_each" else l2(np.hstack(mats))
    res = []
    for s in (1, 2, 3):
        km = KMeans(n_clusters=5, n_init=100, random_state=s).fit(X)
        res.append((km.inertia_, s, km.labels_))
    best = min(res, key=lambda r: r[0])
    a_all = [ari(sci, r[2]) for r in res]
    print("%-46s %5d %8.4f %8.4f %8.4f %8.4f  %-22s %s" % (
        name, best[1], ari(sci, best[2]), nmi(sci, best[2]), ari(rao, best[2]), nmi(rao, best[2]),
        "%.3f-%.3f" % (min(a_all), max(a_all)), sorted(np.bincount(best[2]).tolist(), reverse=True)), flush=True)
