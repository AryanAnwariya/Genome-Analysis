#!/usr/bin/env python3
"""
Cluster graph embeddings with the SCI recipe and evaluate against
Rao et al. 2014 subcompartments and the SCI authors' published GM12878 labels.

SCI recipe (sci/Compartments.py): concatenate LINE order-1 and order-2
embeddings, L2-normalise rows, KMeans(k, random_state=12345).
n_init=10 is set explicitly (scikit-learn 0.19 default used by SCI).
Multiple embedding files are joined by node ID, not by row order.
"""
import argparse
import gzip

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import (adjusted_rand_score, davies_bouldin_score,
                             normalized_mutual_info_score, silhouette_score)
from sklearn.preprocessing import normalize

p = argparse.ArgumentParser()
p.add_argument("--emb", nargs="+", required=True,
               help="embedding file(s), word2vec text format; joined by node ID")
p.add_argument("--nodes", required=True)
p.add_argument("--rao", required=True)
p.add_argument("--sci", required=True)
p.add_argument("--k", type=int, default=5)
p.add_argument("--seed", type=int, default=12345)
p.add_argument("--n-init", type=int, default=10)
p.add_argument("--no-l2", action="store_true")
p.add_argument("--l2-each", action="store_true", help="L2-normalise each embedding file before concatenation (official LINE pipeline)")
p.add_argument("--out-prefix", required=True)
a = p.parse_args()

MAIN = {"A1", "A2", "B1", "B2", "B3"}


def read_emb(path):
    with open(path) as f:
        n, d = map(int, f.readline().split())
        ids, vecs = [], []
        for line in f:
            parts = line.split()
            ids.append(int(parts[0]))
            vecs.append(np.array(parts[1:], dtype=float))
    df = pd.DataFrame(np.vstack(vecs), index=ids)
    print("%s: header %d x %d, read %d x %d" % (path, n, d, df.shape[0], df.shape[1]))
    return df


def read_bed(path, res):
    opener = gzip.open if path.endswith(".gz") else open
    lab = {}
    with opener(path, "rt") as f:
        for line in f:
            c, s, e, l = line.rstrip("\n").split("\t")[:4]
            for b in range(int(s), int(e), res):
                lab[(c, b)] = l
    return lab


def ari_nmi(x, y):
    return adjusted_rand_score(x, y), normalized_mutual_info_score(x, y)


# ---- embeddings, joined by node ID ----
embs = [read_emb(e) for e in a.emb]
common = embs[0].index
for e in embs[1:]:
    common = common.intersection(e.index)
common = common.sort_values()
X = np.hstack([normalize(e.loc[common].values) if a.l2_each else e.loc[common].values for e in embs])
if not a.no_l2:
    X = normalize(X, norm="l2")
print("Embedding matrix: %d nodes x %d dims" % X.shape)

# ---- clustering ----
km = KMeans(n_clusters=a.k, random_state=a.seed, n_init=a.n_init)
clus = km.fit_predict(X)

# ---- labels ----
nodes = pd.read_csv(a.nodes, sep="\t").set_index("node")
res_bp = int(nodes["end"].iloc[0] - nodes["start"].iloc[0])
rao = read_bed(a.rao, res_bp)
sci = read_bed(a.sci, res_bp)

out = nodes.loc[common].copy()
out["cluster"] = clus
keys = list(zip(out["chrom"], out["start"]))
out["rao"] = [rao.get(k, "none") for k in keys]
out["sci_pub"] = [sci.get(k, "none") for k in keys]
out.to_csv(a.out_prefix + "_clusters.tsv", sep="\t", index_label="node")

# ---- metrics ----
lines = [
    "embedding files: " + ", ".join(a.emb),
    "nodes clustered: %d, dims: %d, k=%d, seed=%d, n_init=%d, l2=%s"
    % (X.shape[0], X.shape[1], a.k, a.seed, a.n_init, not a.no_l2),
    "cluster sizes: %s" % np.bincount(clus).tolist(),
    "kmeans inertia: %.4f" % km.inertia_,
    "silhouette (euclidean): %.6f" % silhouette_score(X, clus),
    "davies-bouldin: %.6f" % davies_bouldin_score(X, clus),
]

m = out["rao"].isin(MAIN)
ari, nmi = ari_nmi(out.loc[m, "rao"], out.loc[m, "cluster"])
lines.append("vs Rao A1-B3 (%d nodes): ARI=%.4f NMI=%.4f" % (m.sum(), ari, nmi))

m2 = out["sci_pub"] != "none"
ari, nmi = ari_nmi(out.loc[m2, "sci_pub"], out.loc[m2, "cluster"])
lines.append("vs SCI published (%d nodes): ARI=%.4f NMI=%.4f" % (m2.sum(), ari, nmi))

shared = [k for k in sci if rao.get(k) in MAIN]
ari, nmi = ari_nmi([rao[k] for k in shared], [sci[k] for k in shared])
lines.append("REFERENCE published SCI vs Rao (%d bins): ARI=%.4f NMI=%.4f"
             % (len(shared), ari, nmi))

with open(a.out_prefix + "_metrics.txt", "w") as f:
    f.write("\n".join(lines) + "\n\n")
    f.write("Clusters vs Rao (A1-B3 nodes):\n")
    f.write(pd.crosstab(out.loc[m, "cluster"], out.loc[m, "rao"]).to_string() + "\n\n")
    f.write("Clusters vs SCI published:\n")
    f.write(pd.crosstab(out.loc[m2, "cluster"], out.loc[m2, "sci_pub"]).to_string() + "\n")
print("\n".join(lines))
print("Done.")
