#!/usr/bin/env python3
"""junctions.py -- coupling rank of a junction block J: the rank of the bilinear form
(u, v) -> H(u J v) over all prefixes u and suffixes v, computed as rank(Q M_J P) with P a basis of the
reachable space (26 vectors for span 2) and Q a basis of the co-reachable functionals (14).
Rank 1 = ordinal sum, rank 3 = one-arc join; which ranks occur among junctions of length <= 4?"""
import itertools, numpy as np
from ribbon import transfer_system, letters, word_str
w = 2
order, states, mats, final = transfer_system(w)
d = len(order); L = letters(w)
M = {}
for a in L:
    A = np.zeros((d, d), dtype=np.float64)
    for (j, i), m in mats[a].items(): A[j, i] = m
    M[a] = A
phi = np.zeros(d); phi[final] = 1
def span_basis(vectors):
    B = []
    for v in vectors:
        if B:
            if np.linalg.matrix_rank(np.array(B + [v]), tol=1e-9) > len(B): B.append(v)
        elif np.any(v): B.append(v)
    return np.array(B)
cur = [np.eye(d)[0]]; allv = list(cur)
for n in range(1, 8):
    cur = [M[a] @ v for v in cur for a in L]
    allv.extend(cur)
    if len(allv) > 6000: break
P = span_basis(allv).T
cur = [phi]; allf = [phi]
for n in range(1, 8):
    cur = [f @ M[a] for f in cur for a in L]
    allf.extend(cur)
    if len(allf) > 6000: break
Q = span_basis(allf)
print("reachable dim", P.shape[1], "co-reachable dim", Q.shape[0], "full coupling rank", np.linalg.matrix_rank(Q @ P))
res = {}
for k in range(1, 5):
    for J in itertools.product(L, repeat=k):
        MJ = np.eye(d)
        for a in J: MJ = M[a] @ MJ
        rk = np.linalg.matrix_rank(Q @ MJ @ P, tol=1e-6)
        res.setdefault(rk, []).append(word_str(list(J), w))
for rk in sorted(res):
    print(f"coupling rank {rk}: {len(res[rk])} junctions, e.g. {res[rk][:10]}")
