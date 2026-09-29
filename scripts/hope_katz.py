#!/usr/bin/env python3
"""
HOPE (Ou et al. 2016) with Katz similarity, following the GEM implementation:
  S = (I - beta*A)^-1 (beta*A);  rank-(d/2) SVD  S ~ U s V^T;
  embedding = [U*sqrt(s), V*sqrt(s)]  (d dims in total).
Katz requires beta < 1/spectral_radius(A). GEM's default beta=0.01 diverges on
this dense weighted Hi-C graph, so beta = beta_frac / spectral_radius (default 0.5).
Input: CSR cache written by weighted_deepwalk_csr.py (same graph as LINE).
"""
import argparse
import time

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh, svds

p = argparse.ArgumentParser()
p.add_argument("--csr", required=True)
p.add_argument("--out", required=True)
p.add_argument("--dim", type=int, default=100)
p.add_argument("--beta-frac", type=float, default=0.5)
a = p.parse_args()
t0 = time.time()


def log(msg):
    print("[%8.1fs] %s" % (time.time() - t0, msg), flush=True)


z = np.load(a.csr)
n = len(z["indptr"]) - 1
A = sp.csr_matrix((z["w"].astype(np.float64), z["dst"], z["indptr"]), shape=(n, n))
keep = np.flatnonzero(np.asarray(A.sum(axis=1)).ravel() > 0)
A = A[keep][:, keep]
m = len(keep)
log("nodes kept: %d of %d, nonzeros: %d" % (m, n, A.nnz))

rho = eigsh(A, k=1, which="LA", return_eigenvectors=False)[0]
beta = a.beta_frac / rho
log("spectral radius: %.4f, beta: %.6e" % (rho, beta))

Ad = A.toarray()
del A
M = np.eye(m) - beta * Ad
Ad *= beta
S = np.linalg.solve(M, Ad)
del M, Ad
log("Katz matrix computed (%d x %d)" % S.shape)

k = a.dim // 2
U, s, Vt = svds(S, k=k)
order = np.argsort(s)[::-1]
U, s, Vt = U[:, order], s[order], Vt[order, :]
X = np.hstack([U * np.sqrt(s), Vt.T * np.sqrt(s)])
log("top singular values: %s" % np.round(s[:5], 4).tolist())

with open(a.out, "w") as f:
    f.write("%d %d\n" % X.shape)
    for node, row in zip(keep, X):
        f.write("%d %s\n" % (node, " ".join("%.6f" % v for v in row)))
log("saved %d x %d embedding to %s" % (X.shape[0], X.shape[1], a.out))
print("Done.")
