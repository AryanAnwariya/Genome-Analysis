#!/usr/bin/env python3
"""
Weighted DeepWalk on the SCI Hi-C graph (the same graph file LINE uses).

Next node is chosen with probability proportional to edge weight
(SCI weight = ln(normalised contact)). Walks are embedded with Word2Vec
skip-gram (gensim). Defaults match the demo run: walk length 40,
5 walks per node, dim 100, window 5, sg=1, epochs 5, workers 1, seed 12345.
The graph is held as CSR numpy arrays so that ~455M edges fit in memory.
"""
import argparse
import os
import time

import numpy as np
import pandas as pd
from gensim.models import Word2Vec

p = argparse.ArgumentParser()
p.add_argument("--graph", required=True)
p.add_argument("--out", required=True, help="embedding output, word2vec text format")
p.add_argument("--csr-cache", default=None, help=".npz to save/load the CSR graph")
p.add_argument("--walk-length", type=int, default=40)
p.add_argument("--walks-per-node", type=int, default=5)
p.add_argument("--dim", type=int, default=100)
p.add_argument("--window", type=int, default=5)
p.add_argument("--epochs", type=int, default=5)
p.add_argument("--workers", type=int, default=1)
p.add_argument("--seed", type=int, default=12345)
p.add_argument("--chunksize", type=int, default=20000000)
a = p.parse_args()
t0 = time.time()


def log(msg):
    print("[%8.1fs] %s" % (time.time() - t0, msg), flush=True)


def build_csr(path):
    s_parts, d_parts, w_parts, total = [], [], [], 0
    for ch in pd.read_csv(path, sep="\t", header=None, names=["s", "d", "w"],
                          dtype={"s": np.int32, "d": np.int32, "w": np.float32},
                          chunksize=a.chunksize):
        s_parts.append(ch["s"].values)
        d_parts.append(ch["d"].values)
        w_parts.append(ch["w"].values)
        total += len(ch)
        log("read %d edges" % total)
    src = np.concatenate(s_parts)
    del s_parts
    dst = np.concatenate(d_parts)
    del d_parts
    w = np.concatenate(w_parts)
    del w_parts
    n = int(max(src.max(), dst.max())) + 1
    log("sorting %d edges by source" % len(src))
    order = np.argsort(src, kind="stable")
    dst = dst[order]
    w = w[order]
    counts = np.bincount(src, minlength=n)
    del src, order
    indptr = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(counts, out=indptr[1:])
    return indptr, dst, w


if a.csr_cache and os.path.exists(a.csr_cache):
    z = np.load(a.csr_cache)
    indptr, dst, w = z["indptr"], z["dst"], z["w"]
    log("loaded CSR cache %s" % a.csr_cache)
else:
    indptr, dst, w = build_csr(a.graph)
    if a.csr_cache:
        np.savez(a.csr_cache, indptr=indptr, dst=dst, w=w)
        log("saved CSR cache %s" % a.csr_cache)

n = len(indptr) - 1
cum = np.cumsum(w, dtype=np.float64)
cum0 = np.concatenate(([0.0], cum))
before = cum0[indptr[:-1]]
total_w = cum0[indptr[1:]] - before
starts = np.flatnonzero(total_w > 0)
log("nodes: %d, with positive-weight edges: %d, directed edges: %d, zero-weight edges: %d"
    % (n, len(starts), len(dst), int((w == 0).sum())))

# ---- weighted random walks, all walkers advanced together ----
rng = np.random.default_rng(a.seed)
first = np.concatenate([rng.permutation(starts) for _ in range(a.walks_per_node)])
W = len(first)
walks = np.full((W, a.walk_length), -1, dtype=np.int32)
walks[:, 0] = first
cur = first.astype(np.int64)
alive = np.ones(W, dtype=bool)
for step in range(1, a.walk_length):
    alive &= total_w[cur] > 0
    idx_alive = np.flatnonzero(alive)
    c = cur[idx_alive]
    r = before[c] + rng.random(len(c)) * total_w[c]
    e = np.searchsorted(cum, r, side="right")
    e = np.minimum(e, indptr[c + 1] - 1)
    nxt = dst[e].astype(np.int64)
    walks[idx_alive, step] = nxt
    cur[idx_alive] = nxt
truncated = int((walks[:, -1] < 0).sum())
log("walks: %d x length %d, truncated early: %d" % (W, a.walk_length, truncated))

sentences = [[str(x) for x in row if x >= 0] for row in walks]
del walks, cum, cum0, dst, w
log("training Word2Vec")
model = Word2Vec(sentences, vector_size=a.dim, window=a.window, min_count=1,
                 sg=1, workers=a.workers, epochs=a.epochs, seed=a.seed)
model.wv.save_word2vec_format(a.out)
log("saved %d embeddings to %s" % (len(model.wv), a.out))
print("Done.")
