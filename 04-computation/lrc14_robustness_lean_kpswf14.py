#!/usr/bin/env python3
"""
kind-pasteur-2026-06-22 (kpswf14) LEAN robustness check (fast).

600 curated loose rows; 3 lead magnitude-aware separators. Optimized H (bit-packed DP).
Collision = loose row whose (c3,H) or H equals AP's or GW's value under that candidate.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce
import random
random.seed(7)

def H_paths_fast(adj):
    n = len(adj); size = 1 << n
    dp = [0] * (size * n)
    for v in range(n): dp[(1 << v) * n + v] = 1
    for mask in range(1, size):
        base = mask * n
        for v in range(n):
            c = dp[base + v]
            if c == 0: continue
            outs = adj[v] & ~mask
            while outs:
                w = (outs & -outs).bit_length() - 1
                outs &= outs - 1
                dp[(mask | (1 << w)) * n + w] += c
    fb = (size - 1) * n
    return sum(dp[fb:fb + n])

def adj_c3(A):
    n = len(A)
    adj = [sum(1 << j for j in range(n) if A[i][j]) for i in range(n)]
    c = 0
    for a, b, cc in combinations(range(n), 3):
        if A[a][b] and A[b][cc] and A[cc][a]: c += 1
        if A[a][cc] and A[cc][b] and A[b][a]: c += 1
    return adj, c

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

CANDS = {"floor-odd": floor_odd, "reciprocal-stk": reciprocal_stack, "CF-parity": cf_parity}

def primitive(S): return reduce(gcd, S) == 1
def sg(fn, S):
    adj, c = adj_c3(fn(S)); return (c, H_paths_fast(adj))

AP = list(range(1, 14)); GW = list(range(1, 12)) + [13, 24]
tH = {}; tC3H = {}
for cn, fn in CANDS.items():
    a = sg(fn, AP); g = sg(fn, GW)
    tH[cn] = {a[1], g[1]}; tC3H[cn] = {a, g}
    print(f"{cn:14s}: AP=(c3={a[0]},H={a[1]})  GW=(c3={g[0]},H={g[1]})")

rows = set()
for v in range(14, 121):
    S = tuple(sorted(set(list(range(1, 12)) + [13, v])))
    if len(S) == 13 and primitive(S): rows.add(S)
for v in range(13, 121):
    S = tuple(sorted(set(list(range(1, 13)) + [v])))
    if len(S) == 13 and primitive(S): rows.add(S)
for rem in range(1, 14):
    for ins in range(14, 30):
        S = tuple(sorted(set(range(1, 14)) - {rem} | {ins}))
        if len(S) == 13 and primitive(S): rows.add(S)
for _ in range(400):
    S = tuple(sorted(random.sample(range(1, 40), 13)))
    if primitive(S): rows.add(S)
rows.discard(tuple(AP)); rows.discard(tuple(GW))
rows = list(rows)
print(f"\nloose bank: {len(rows)} rows")

cH = {cn: 0 for cn in CANDS}; cC = {cn: [] for cn in CANDS}
for S in rows:
    S = list(S)
    for cn, fn in CANDS.items():
        s = sg(fn, S)
        if s[1] in tH[cn]:
            cH[cn] += 1
            if s in tC3H[cn]: cC[cn].append(tuple(S))

print("\nROBUSTNESS:")
print(f"{'cand':14s} | {'H-coll':>7s} | {'(c3,H)-coll':>12s}")
for cn in CANDS:
    print(f"{cn:14s} | {cH[cn]:>7d} | {len(cC[cn]):>12d}")
for cn in CANDS:
    if cC[cn]: print(f"  {cn} full-sig collisions: {cC[cn][:5]}")
