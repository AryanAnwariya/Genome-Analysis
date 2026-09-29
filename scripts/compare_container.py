#!/usr/bin/env python3
import numpy as np, pandas as pd
from sklearn.metrics import adjusted_rand_score as ari, normalized_mutual_info_score as nmi

ref = pd.read_csv("results/gm12878_interkr/ninit100/DW_HOPE_s1_clusters.tsv", sep="\t", keep_default_na=False)
ref["key"] = ref["chrom"] + ":" + ref["start"].astype(str)
ref = ref.set_index("key")
pub = ref["sci_pub"].where(ref["sci_pub"] != "none")
rao = ref["rao"].where(ref["rao"].isin(["A1", "A2", "B1", "B2", "B3"]))

def read_lab(p):
    d = {}
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) < 4 or not f[1].isdigit():
            continue
        d[f[0] + ":" + f[1]] = f[3]
    return pd.Series(d)

A = read_lab("results/container_demo/demo_compartments.txt")
B = read_lab("results/container_demo/saved_demo_run_compartments.txt")
bad = {"NA", "nan", "NaN", ""}
A = A[~A.isin(bad)]
B = B[~B.isin(bad)]
print("labelled bins: fresh run %d | saved run %d | published %d | Rao A1-B3 %d" % (len(A), len(B), pub.notna().sum(), rao.notna().sum()))
print("cluster sizes fresh:", sorted(A.value_counts().tolist(), reverse=True))
print("cluster sizes saved:", sorted(B.value_counts().tolist(), reverse=True))

def cmp(name, x, y):
    k = x.index.intersection(y.dropna().index)
    if len(k) == 0:
        print(name, "no overlap")
        return
    print("%-34s bins %6d  ARI %.4f  NMI %.4f" % (name, len(k), ari(x[k], y[k]), nmi(x[k], y[k])))

cmp("fresh run vs saved run", A, B)
cmp("fresh run vs published SCI", A, pub)
cmp("saved run vs published SCI", B, pub)
cmp("fresh run vs Rao A1-B3", A, rao)
cmp("saved run vs Rao A1-B3", B, rao)
cmp("published SCI vs Rao (reference)", pub.dropna(), rao)
k = A.index.intersection(pub.dropna().index)
print("\nfresh run (rows) vs published (columns):")
print(pd.crosstab(A[k], pub[k]))
