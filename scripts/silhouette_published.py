#!/usr/bin/env python3
import os
import numpy as np, pandas as pd
from sklearn.metrics import silhouette_score, davies_bouldin_score

base = pd.read_csv("results/gm12878_interkr/ninit100/DW_HOPE_s1_clusters.tsv", sep="\t", keep_default_na=False)
ref = base[base["sci_pub"] != "none"]
ids = ref["node"].astype(str).values
sci = ref["sci_pub"].values
rao = ref["rao"].values
print("bins:", len(ref))

def read_emb(p):
    df = pd.read_csv(p, sep=r"\s+", skiprows=1, header=None, index_col=0)
    df.index = df.index.astype(str)
    return df.loc[ids].values.astype(float)

def l2(x):
    return x / np.maximum(np.linalg.norm(x, axis=1, keepdims=True), 1e-12)

I = "results/gm12878_interkr/"
O = "results/gm12878_interkr_oe/"
L = I + "gm12878_INTERKR_order_%s_samples_25M.embedding"
variants = [
    ("INTERKR LINE order 3", [L % "3"], "cat"),
    ("INTERKR LINE 1+2 concatenated", [L % "1", L % "2"], "cat"),
    ("INTERKR LINE 1+2 per-order L2", [L % "1", L % "2"], "each"),
    ("INTERKR LINE order 1 only", [L % "1"], "cat"),
    ("INTERKR weighted DeepWalk", [I + "gm12878_INTERKR_weighted_deepwalk.embedding"], "cat"),
    ("INTERKR HOPE", [I + "gm12878_INTERKR_hope_katz.embedding"], "cat"),
    ("Pair-O/E LINE order 3", [O + "gm12878_INTERKR_OE_order_3_samples_25M.embedding"], "cat"),
    ("Pair-O/E weighted DeepWalk", [O + "gm12878_INTERKR_OE_weighted_deepwalk.embedding"], "cat"),
    ("Pair-O/E HOPE", [O + "gm12878_INTERKR_OE_hope_katz.embedding"], "cat"),
]
print("%-32s | %8s %8s %6s | %8s %8s %6s" % ("Embedding", "silh SCI", "silh Rao", "ratio", "DB SCI", "DB Rao", "ratio"))
for name, files, mode in variants:
    if not all(os.path.exists(f) for f in files):
        print("%-32s SKIPPED (missing: %s)" % (name, [f for f in files if not os.path.exists(f)]))
        continue
    mats = [read_emb(f) for f in files]
    X = np.hstack([l2(m) for m in mats]) if mode == "each" else l2(np.hstack(mats))
    ss = silhouette_score(X, sci, sample_size=10000, random_state=0)
    sr = silhouette_score(X, rao, sample_size=10000, random_state=0)
    ds = davies_bouldin_score(X, sci)
    dr = davies_bouldin_score(X, rao)
    print("%-32s | %8.4f %8.4f %6.2f | %8.3f %8.3f %6.2f" % (name, ss, sr, ss / sr if sr else float("nan"), ds, dr, ds / dr), flush=True)
