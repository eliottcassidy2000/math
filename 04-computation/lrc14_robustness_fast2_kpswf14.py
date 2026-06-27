#!/usr/bin/env python3
"""
kind-pasteur-2026-06-22 (kpswf14) FAST robustness bank v2 (optimized H).

Tests the 4 magnitude-aware separators over a curated loose bank of ~2000 rows.
Question: how often does a loose row's iso-fingerprint collide with a TIGHT (AP or GW) value?
A robust separator -> few collisions.

Speed: H counted by an optimized bitmask DP using bit-packed adjacency (one int per row) and
iterating subsets in increasing popcount-friendly order. Only loose rows are scanned (we know
from the proven census that AP, GW are the only tight rows here).
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce
import random
random.seed(99)

def H_paths_fast(adj):
    """adj[v] = bitmask of out-neighbors. Count directed Hamiltonian paths, n<=13."""
    n = len(adj)
    size = 1 << n
    dp = [0] * (size * n)  # dp[mask*n + v]
    for v in range(n):
        dp[(1 << v) * n + v] = 1
    for mask in range(1, size):
        base = mask * n
        for v in range(n):
            c = dp[base + v]
            if c == 0:
                continue
            outs = adj[v] & ~mask
            while outs:
                w = (outs & -outs).bit_length() - 1
                outs &= outs - 1
                dp[(mask | (1 << w)) * n + w] += c
    fb = (size - 1) * n
    return sum(dp[fb:fb + n])

def to_adj(A):
    n = len(A)
    return [sum(1 << j for j in range(n) if A[i][j]) for i in range(n)]

def c3(A):
    n = len(A); cnt = 0
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]: cnt += 1
        if A[a][c] and A[c][b] and A[b][a]: cnt += 1
    return cnt

def floor_odd(S):
    k = len(S); A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            ai = (S[i] // S[j]) % 2 == 1; aj = (S[j] // S[i]) % 2 == 1
            if ai and not aj: A[i][j] = 1
            elif aj and not ai: A[j][i] = 1
            else: A[i][j] = 1 if S[i] < S[j] else 0; A[j][i] = 1 - A[i][j]
    return A

def reciprocal_stack(S):
    k = len(S); A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            votes = 0; tot = 0; di = S[i] - S[j]
            for s in S:
                rel = (F(di, s)) % 1
                if rel == 0 or rel == F(1, 2): continue
                tot += 1
                if rel < F(1, 2): votes += 1
            if 2*votes > tot: A[i][j] = 1
            elif 2*votes < tot: A[j][i] = 1
            else: A[i][j] = 1 if S[i] < S[j] else 0; A[j][i] = 1 - A[i][j]
    return A

def cf_parity(S):
    def cl(a, b):
        n = 0
        while b: a, b = b, a % b; n += 1
        return n
    k = len(S); A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            L = cl(S[i], S[j]); small = (L % 2 == 1)
            if small: A[i][j] = 1 if S[i] < S[j] else 0
            else: A[i][j] = 1 if S[i] > S[j] else 0
            A[j][i] = 1 - A[i][j]
    return A

def divis_lt_half(S):
    k = len(S); raw = [[2*(S[i] % S[j]) < S[j] for j in range(k)] for i in range(k)]
    A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            a, b = raw[i][j], raw[j][i]
            if a and not b: A[i][j] = 1
            elif b and not a: A[j][i] = 1
            else: A[i][j] = 1 if S[i] < S[j] else 0; A[j][i] = 1 - A[i][j]
    return A

CANDS = {
    "floor-odd     ": floor_odd,
    "reciprocal-stk": reciprocal_stack,
    "CF-parity     ": cf_parity,
    "divis<s/2     ": divis_lt_half,
}

def primitive(S): return reduce(gcd, S) == 1
def sig(fn, S):
    A = fn(S); return (c3(A), H_paths_fast(to_adj(A)))

AP = list(range(1, 14)); GW = list(range(1, 12)) + [13, 24]
tightsigH = {}; tightsigC3H = {}
print("TIGHT signatures:")
for cn, fn in CANDS.items():
    sa = sig(fn, AP); sg = sig(fn, GW)
    tightsigC3H[cn] = {sa, sg}
    tightsigH[cn] = {sa[1], sg[1]}
    print(f"  {cn}: AP=(c3={sa[0]},H={sa[1]})  GW=(c3={sg[0]},H={sg[1]})")

# curated loose bank ~2000
rows = set()
for v in range(14, 161):  # GW petal
    S = tuple(sorted(set(list(range(1, 12)) + [13, v])))
    if len(S) == 13 and primitive(S): rows.add(S)
for v in range(13, 161):  # top swap
    S = tuple(sorted(set(list(range(1, 13)) + [v])))
    if len(S) == 13 and primitive(S): rows.add(S)
for rem in range(1, 14):  # single replacement
    for ins in range(14, 45):
        S = tuple(sorted(set(range(1, 14)) - {rem} | {ins}))
        if len(S) == 13 and primitive(S): rows.add(S)
for _ in range(1500):  # random
    S = tuple(sorted(random.sample(range(1, 45), 13)))
    if primitive(S): rows.add(S)
rows.discard(tuple(AP)); rows.discard(tuple(GW))
rows = list(rows)
print(f"\nloose bank: {len(rows)} rows")

collH = {cn: [] for cn in CANDS}; collC3H = {cn: [] for cn in CANDS}
for idx, S in enumerate(rows):
    S = list(S)
    for cn, fn in CANDS.items():
        s = sig(fn, S)
        if s[1] in tightsigH[cn]:
            collH[cn].append(tuple(S))
            if s in tightsigC3H[cn]: collC3H[cn].append(tuple(S))

print("\nROBUSTNESS (loose collisions onto a TIGHT value):")
print(f"{'candidate':16s} | {'H-coll':>7s} | {'(c3,H)-coll':>12s} | rate(c3,H)")
ranked = []
for cn in CANDS:
    nH = len(collH[cn]); nc = len(collC3H[cn])
    ranked.append((nc, nH, cn))
    print(f"{cn:16s} | {nH:>7d} | {nc:>12d} | {nc/len(rows):.4f}")
ranked.sort()
print(f"\nCLEANEST separator by (c3,H)-collision rate: {ranked[0][2].strip()} "
      f"({ranked[0][0]} collisions / {len(rows)})")

print("\nsample (c3,H)-collisions:")
for cn in CANDS:
    print(f"  {cn}: {collC3H[cn][:5]}")
