#!/usr/bin/env python3
import os
import pandas as pd
from sklearn.metrics import adjusted_rand_score as ari

def read_lab(p):
    d = {}
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) < 4 or not f[1].isdigit():
            continue
        d[f[0] + ":" + f[1]] = f[3]
    s = pd.Series(d)
    return s[~s.isin({"NA", "nan", "NaN", ""})]

saved = read_lab("results/container_demo/saved_demo_run_compartments.txt")
fresh = read_lab("results/container_demo/demo_compartments.txt")
runs = [("first run: -o both -s 25", fresh)]
for n, lab in [("o1_s25", "-o 1 -s 25"), ("both_s1", "-o both -s 1"), ("both_s25_rep", "-o both -s 25 (repeat)")]:
    p = "results/container_demo_%s/demo_compartments.txt" % n
    if os.path.exists(p):
        runs.append((lab, read_lab(p)))
    else:
        print("missing:", p)
print("%-26s %8s %14s %14s   %s" % ("run", "bins", "ARI vs saved", "ARI vs first", "cluster sizes"))
for lab, s in runs:
    k = s.index.intersection(saved.index)
    k2 = s.index.intersection(fresh.index)
    print("%-26s %8d %14.4f %14.4f   %s" % (lab, len(s), ari(s[k], saved[k]), ari(s[k2], fresh[k2]), sorted(s.value_counts().tolist(), reverse=True)))
