#!/usr/bin/env python3
"""Audit A: N(28) without materialising S_27 (memory-lean extension of a2_census.py).
T_28 = 3 S_27 + 2^28 with S_27 = S_26 u T_27, so T_28 = (3 S_26 + 2^28) u (3 T_27 + 2^28) and
hits_28 = T_28 n S_27 is the union of four sorted-array intersections, done in chunks."""
import math, time
import numpy as np

t0 = time.time()
S = np.array([-1], dtype=np.int64)
P = np.array([], dtype=np.int64)
for e in range(1, 27):
    T = 3 * S + (1 << e)
    idx = np.searchsorted(S, T); idx[idx == len(S)] = len(S) - 1
    hit = T[S[idx] == T]; del idx
    P = np.union1d(P, hit)
    S = np.union1d(S, T); del T
# now S = S_26, P = P_26
T27 = 3 * S + (1 << 27)
idx = np.searchsorted(S, T27); idx[idx == len(S)] = len(S) - 1
hit27 = T27[S[idx] == T27]; del idx
P27 = np.union1d(P, hit27)
print(f"N(26) = {len(P)}, N(27) = {len(P27)}  |S_26| = {len(S)}, |T_27| = {len(T27)}  [{time.time()-t0:.1f}s]", flush=True)


def inter(A_src_mult, A_src, B, chunk=4_000_000):
    """sorted intersection of (3*A_src + 2^28) with sorted unique B, chunked"""
    out = []
    for i in range(0, len(A_src), chunk):
        X = 3 * A_src[i:i + chunk] + (1 << 28)
        j = np.searchsorted(B, X); j[j == len(B)] = len(B) - 1
        out.append(X[B[j] == X])
    return np.concatenate(out) if out else np.array([], dtype=np.int64)


hits = []
for src in (S, T27):
    for B in (S, T27):
        hits.append(inter(3, src, B))
hit28 = np.unique(np.concatenate(hits))
P28 = np.union1d(P27, hit28)
for E, N in ((26, len(P)), (27, len(P27)), (28, len(P28))):
    print(f"E={E} N(E)={N}  N sqrt(E)/2^E = {N*math.sqrt(E)/2**E:.6f}")
b = -math.log(len(P28) / (2 * len(P27))) / math.log(28 / 27)
print(f"beta_local(28) = {b:.3f}   [{time.time()-t0:.1f}s]")
