#!/usr/bin/env python3
"""
kind-pasteur-2026-06-22 (kpswf14) FAST robustness bank for the magnitude separators.

We already KNOW (proven bounded+single-swap census, HYP-2920/kps-S39) that the ONLY tight
13-rows in these structured families are AP and GW. So instead of recomputing exact M for
hundreds of large-speed rows (slow), we take the loose bank as GIVEN (everything != AP,GW)
and ask the robustness question directly:

  Over a LARGE bank of loose rows, does any loose row's iso-fingerprint (c3,c5,H) under
  candidate C collide with H(AP) or H(GW)?

A robust separator has FEW/ZERO collisions. We rank the 4 magnitude-aware candidates by
collision count. (We DO double check a handful of random rows with exact M to confirm none
are accidentally tight.)

n=13, exact Fraction for reciprocal-stack; H via bitmask DP.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd
from functools import reduce
import random
random.seed(1414)

def H_paths(A):
    n = len(A); size = 1 << n
    dp = [[0]*n for _ in range(size)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(size):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c == 0: continue
            Av = A[v]
            for w in range(n):
                if mask & (1 << w): continue
                if Av[w]: dp[mask | (1 << w)][w] += c
    return sum(dp[size - 1])

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

AP = list(range(1, 14)); GW = list(range(1, 12)) + [13, 24]

# tight signatures per candidate: (c3, H) AND H alone
tightsig = {}
for cn, fn in CANDS.items():
    A_ap = fn(AP); A_gw = fn(GW)
    tightsig[cn] = {
        'H': {H_paths(A_ap), H_paths(A_gw)},
        'c3H': {(c3(A_ap), H_paths(A_ap)), (c3(A_gw), H_paths(A_gw))},
    }
    print(f"{cn}: AP (c3={c3(A_ap)}, H={H_paths(A_ap)})  GW (c3={c3(A_gw)}, H={H_paths(A_gw)})")

# LARGE loose bank (these are loose by the proven census / by being non-AP/GW structured rows)
rows = set()
# GW-petal family: {1..11,13,v}, v in 14..300  (tight ONLY v=24)
for v in range(14, 301):
    S = tuple(sorted(set(list(range(1, 12)) + [13, v])))
    if len(S) == 13 and primitive(S) and S != tuple(GW): rows.add(S)
# top-swap {1..12, v}, v 13..300
for v in range(13, 301):
    S = tuple(sorted(set(list(range(1, 13)) + [v])))
    if len(S) == 13 and primitive(S): rows.add(S)
# AP single replacements rem->ins
for rem in range(1, 14):
    for ins in range(14, 80):
        S = tuple(sorted(set(range(1, 14)) - {rem} | {ins}))
        if len(S) == 13 and primitive(S) and S != tuple(GW): rows.add(S)
# two-swaps {1..11,u,v}
for u in range(12, 50):
    for v in range(u+1, 50):
        S = tuple(sorted(set(list(range(1, 12)) + [u, v])))
        if len(S) == 13 and primitive(S) and S != tuple(GW): rows.add(S)
# random primitive 13-subsets of {1..50}
for _ in range(6000):
    S = tuple(sorted(random.sample(range(1, 51), 13)))
    if primitive(S) and S not in (tuple(AP), tuple(GW)): rows.add(S)

rows.discard(tuple(AP)); rows.discard(tuple(GW))
print(f"\nloose bank size: {len(rows)} rows\n")

collide_H = {cn: [] for cn in CANDS}
collide_c3H = {cn: [] for cn in CANDS}
for S in rows:
    S = list(S)
    for cn, fn in CANDS.items():
        A = fn(S); h = H_paths(A)
        if h in tightsig[cn]['H']:
            collide_H[cn].append(tuple(S))
            cc = c3(A)
            if (cc, h) in tightsig[cn]['c3H']:
                collide_c3H[cn].append(tuple(S))

print("ROBUSTNESS (collisions = loose rows whose invariant equals a TIGHT value):")
print(f"{'candidate':16s} | {'H-collisions':>13s} | {'(c3,H)-collisions':>18s} | collide-rate")
for cn in CANDS:
    nH = len(collide_H[cn]); nc = len(collide_c3H[cn])
    print(f"{cn:16s} | {nH:>13d} | {nc:>18d} | {nH/len(rows):.4f}")

print("\nSample (c3,H)-collisions per candidate (loose rows hitting full tight signature):")
for cn in CANDS:
    print(f"  {cn}: {collide_c3H[cn][:6]}")

# sanity: confirm a few collision rows are genuinely LOOSE (exact M != 1/14)
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def gmin(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def M_exact(S):
    b = F(0)
    for t in cand(S):
        v = gmin(S, t)
        if v > b: b = v
    return b
print("\nSanity (exact M on up-to-5 (c3,H)-collision rows, must be != 1/14):")
shown = 0
for cn in CANDS:
    for S in collide_c3H[cn]:
        if shown >= 5: break
        M = M_exact(list(S)); print(f"  {cn} row {S}: M={M} tight={M==F(1,14)}")
        shown += 1
    if shown >= 5: break
