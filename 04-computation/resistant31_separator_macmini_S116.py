#!/usr/bin/env python3
"""resistant31_separator_macmini_S116.py -- mac-mini-2026-07-16-S116.
THE RESISTANT 31 (boxeph-S31 / owner): at n=7, 256 within-level cospectral class pairs;
score-of-scores (rowsums of A^2) separates 225; 31 pairs resist. FIND THE SEPARATOR.

Census method: enumerate all 2^21 labeled tournaments; group by exact integer charpoly
(class invariant); within each group, peel iso classes by orbit generation (5040 perms).
Then test the invariant tower on the resistant pairs:
  I3  = sorted rowsums of A^3        (score-of-scores-of-scores)
  D3  = sorted diag(A^3)             (per-vertex 3-cycle participation, /1)
  I4  = sorted rowsums of A^4
  TAU = sorted arborescence vector   (Matrix-Tree: det of reduced out-Laplacian per root;
                                      klein cont.9's tau-vector)
  H   = Hamiltonian path count       (the house invariant)
  MIX = sorted rowsums of A A^T      (common-victims profile)
"""
import sys, math
import numpy as np
from itertools import combinations, permutations
sys.stdout.reconfigure(line_buffering=True)

n = 7
pairs = list(combinations(range(n), 2))
m = len(pairs)
N = 1 << m
print(f"enumerating {N} labeled tournaments on {n} vertices...")

# adjacency bit codes -> per-labeled charpoly via Newton identities (batched in chunks)
perms = list(permutations(range(n)))
pidx = {p: i for i, p in enumerate(pairs)}
# permutation action on codes: perm sends pair (u,v)->(pu,pv); bit i (u<v, A[u][v]=bit) maps w/ orientation
perm_maps = []
for pm in perms:
    mp = []
    for (u, v) in pairs:
        a, b = pm[u], pm[v]
        if a < b: mp.append((pidx[(a, b)], 0))
        else: mp.append((pidx[(b, a)], 1))
    perm_maps.append(mp)

def apply_perm(code, mp):
    out = 0
    for i, (j, flip) in enumerate(mp):
        bit = (code >> i) & 1
        if flip: bit ^= 1
        if bit: out |= (1 << j)
    return out

def code_to_A(code):
    A = np.zeros((n, n), dtype=np.int64)
    for i, (u, v) in enumerate(pairs):
        if (code >> i) & 1: A[u, v] = 1
        else: A[v, u] = 1
    return A

def charpoly_key(A):
    # Newton's identities from power traces (exact ints)
    P = [np.trace(np.linalg.matrix_power(A, k)) for k in range(1, n + 1)]
    e = [1]
    for k in range(1, n + 1):
        s = 0
        for i in range(1, k + 1):
            s += (-1) ** (i - 1) * e[k - i] * P[i - 1]
        e.append(s // k)
    return tuple(int(v) for v in e[1:])

# pass 1: group codes by charpoly (chunked, plain python fast enough? 2M x (7 matmuls))
# speed: use numpy batched power traces
CH = 1 << 15
groups = {}
codes_arr = np.arange(N, dtype=np.int64)
bits = ((codes_arr[:, None] >> np.arange(m)[None, :]) & 1).astype(np.int8)
A_all = np.zeros((N, n, n), dtype=np.int8)
for i, (u, v) in enumerate(pairs):
    A_all[:, u, v] = bits[:, i]
    A_all[:, v, u] = 1 - bits[:, i]
del bits
print("computing charpolys (batched)...")
keys = np.zeros((N, n), dtype=np.int64)
for lo in range(0, N, CH):
    hi = min(lo + CH, N)
    A = A_all[lo:hi].astype(np.int64)
    Pk = np.zeros((hi - lo, n), dtype=np.int64)
    Ak = A.copy()
    Pk[:, 0] = np.trace(Ak, axis1=1, axis2=2)
    for k in range(1, n):
        Ak = np.matmul(Ak, A)
        Pk[:, k] = np.trace(Ak, axis1=1, axis2=2)
    e = np.zeros((hi - lo, n + 1), dtype=np.int64); e[:, 0] = 1
    for k in range(1, n + 1):
        s = np.zeros(hi - lo, dtype=np.int64)
        for i in range(1, k + 1):
            s += ((-1) ** (i - 1)) * e[:, k - i] * Pk[:, i - 1]
        e[:, k] = s // k
    keys[lo:hi] = e[:, 1:]
print("grouping...")
from collections import defaultdict
bykey = defaultdict(list)
for c in range(N):
    bykey[tuple(keys[c])].append(c)
print(f"charpoly groups: {len(bykey)}")

# pass 2: peel classes within multi-class groups
print("peeling iso classes within groups...")
classes = []          # (rep_code, charpoly_key)
for key, members in bykey.items():
    mem = set(members)
    while mem:
        rep = min(mem)
        orbit = {apply_perm(rep, mp) for mp in perm_maps}
        mem -= orbit
        classes.append((rep, key))
print(f"total iso classes: {len(classes)} (expect 456)")

# cospectral class pairs (same charpoly)
bykey_cls = defaultdict(list)
for rep, key in classes: bykey_cls[key].append(rep)
def xval(A):
    s = A.sum(1); d = 2 * s - (n - 1)
    return int((d * d).sum())
cos_pairs = []
for key, reps in bykey_cls.items():
    if len(reps) > 1:
        for a, b in combinations(reps, 2):
            Aa, Ab = code_to_A(a), code_to_A(b)
            cos_pairs.append((a, b, xval(Aa) == xval(Ab)))
within = [(a, b) for a, b, sx in cos_pairs if sx]
print(f"cospectral class pairs: {len(cos_pairs)}; within-x-level: {len(within)} (boxeph: 256)")

def inv_I(A, k):
    return tuple(sorted(np.linalg.matrix_power(A.astype(np.int64), k).sum(1)))
def inv_D3(A):
    return tuple(sorted(np.diag(np.linalg.matrix_power(A.astype(np.int64), 3))))
def inv_TAU(A):
    L = np.diag(A.sum(1)) - A          # out-degree Laplacian
    taus = []
    for r in range(n):
        M = np.delete(np.delete(L, r, 0), r, 1).astype(float)
        taus.append(int(round(np.linalg.det(M))))
    return tuple(sorted(taus))
def inv_H(A):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1 << n):
        for v in range(n):
            c = dp[S][v]
            if not c or not (S >> v) & 1: continue
            for u in range(n):
                if not (S >> u) & 1 and A[v][u]:
                    dp[S | 1 << u][u] += c
    return sum(dp[(1 << n) - 1][v] for v in range(n))
def inv_MIX(A):
    return tuple(sorted((A @ A.T).sum(1)))

# resistant under score-of-scores
resist = []
for a, b in within:
    Aa, Ab = code_to_A(a), code_to_A(b)
    if inv_I(Aa, 2) == inv_I(Ab, 2):
        resist.append((a, b))
print(f"score-of-scores separates {len(within) - len(resist)}; RESISTANT: {len(resist)} (boxeph: 31)")

towers = [("I3 rowsums A^3", lambda A: inv_I(A, 3)),
          ("D3 diag A^3 (3-cycle participation)", inv_D3),
          ("I4 rowsums A^4", lambda A: inv_I(A, 4)),
          ("TAU arborescence vector", inv_TAU),
          ("MIX rowsums A A^T", inv_MIX),
          ("H hamiltonian paths", inv_H)]
print("\nSEPARATOR HUNT on the resistant pairs:")
for name, f in towers:
    sep = sum(1 for a, b in resist if f(code_to_A(a)) != f(code_to_A(b)))
    print(f"   {name}: separates {sep}/{len(resist)}")
# best pairings
print("\ncombination check (first invariant that separates each pair):")
remaining = list(resist)
for name, f in towers:
    still = [(a, b) for a, b in remaining if f(code_to_A(a)) == f(code_to_A(b))]
    print(f"   after {name}: {len(still)} remain")
    remaining = still
    if not remaining: break
if remaining:
    print(f"   UNSEPARATED pairs (codes): {remaining[:5]}")
print("\nDONE")
