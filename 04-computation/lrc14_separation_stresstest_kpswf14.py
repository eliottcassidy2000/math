#!/usr/bin/env python3
"""
kind-pasteur-2026-06-22 (kpswf14) STRESS TEST of the two candidates that separated the
small bank: 'divis floor-odd' and 'reciprocal-stack'.

We throw a LARGE bank of loose 13-rows and ask: does any loose row land on a tight H-value
(H(AP) or H(GW))?  If separation is real & robust, NO loose row should collide.  If it was
sparse-bank luck, collisions will appear quickly.

Bank: every primitive 13-subset of {1..M} reachable by:
  - {1..12, v}  for v in 13..400          (single top-swap family, all loose except none here)
  - {1..11, 13, v} for v in 14..400       (GW-petal family, tight only at v=24)
  - {1..11, u, v} two-swaps with u,v <= 60
  - random primitive 13-subsets of {1..40}
Tight rows are detected by exact M and excluded from 'loose'; we verify AP,GW are the ONLY
tight ones found (consistent with the proven bounded census).
Exact arithmetic. H via bitmask DP.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce
import random
random.seed(14)

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
TIGHT = F(1, 14)

def divis_floor_odd(S):
    k = len(S); A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            ai = (S[i] // S[j]) % 2 == 1
            aj = (S[j] // S[i]) % 2 == 1
            if ai and not aj: A[i][j] = 1
            elif aj and not ai: A[j][i] = 1
            else:
                if S[i] < S[j]: A[i][j] = 1
                else: A[j][i] = 1
    return A

def reciprocal_stack(S):
    k = len(S); A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            votes = 0; tot = 0
            di = S[i] - S[j]
            for s in S:
                rel = (F(di, s)) % 1
                if rel == 0 or rel == F(1, 2): continue
                tot += 1
                if rel < F(1, 2): votes += 1
            if 2*votes > tot: A[i][j] = 1
            elif 2*votes < tot: A[j][i] = 1
            else:
                if S[i] < S[j]: A[i][j] = 1
                else: A[j][i] = 1
    return A

def primitive(S): return reduce(gcd, S) == 1

AP = list(range(1, 14)); GW = list(range(1, 12)) + [13, 24]
H_AP_fo = H_paths(divis_floor_odd(AP)); H_GW_fo = H_paths(divis_floor_odd(GW))
H_AP_rs = H_paths(reciprocal_stack(AP)); H_GW_rs = H_paths(reciprocal_stack(GW))
print(f"floor-odd      : H(AP)={H_AP_fo}  H(GW)={H_GW_fo}")
print(f"reciprocal-stk : H(AP)={H_AP_rs}  H(GW)={H_GW_rs}")
tightH_fo = {H_AP_fo, H_GW_fo}; tightH_rs = {H_AP_rs, H_GW_rs}

# build large loose bank
rows = set()
for v in range(13, 401):
    S = tuple(sorted(set(list(range(1, 13)) + [v])))
    if len(S) == 13 and primitive(S): rows.add(S)
for v in range(14, 401):
    S = tuple(sorted(set(list(range(1, 12)) + [13, v])))
    if len(S) == 13 and primitive(S): rows.add(S)
# two-swaps {1..11,u,v}
for u in range(12, 61):
    for v in range(u+1, 61):
        S = tuple(sorted(set(list(range(1, 12)) + [u, v])))
        if len(S) == 13 and primitive(S): rows.add(S)
# random primitive 13-subsets of {1..40}
tries = 0
while tries < 4000:
    S = tuple(sorted(random.sample(range(1, 41), 13)))
    tries += 1
    if primitive(S): rows.add(S)

rows.discard(tuple(AP)); rows.discard(tuple(GW))
print(f"\nloose bank size: {len(rows)} rows")

# scan
coll_fo = []; coll_rs = []
extra_tight = []
checked = 0
for S in rows:
    S = list(S)
    # cheap tight check first (only flag if exactly 1/14)
    M = M_exact(S)
    if M == TIGHT:
        extra_tight.append(S)
        continue  # not loose
    checked += 1
    hfo = H_paths(divis_floor_odd(S))
    if hfo in tightH_fo: coll_fo.append((tuple(S), hfo))
    hrs = H_paths(reciprocal_stack(S))
    if hrs in tightH_rs: coll_rs.append((tuple(S), hrs))

print(f"loose rows checked: {checked}")
print(f"extra tight rows found (should be NONE besides AP/GW): {extra_tight}")
print()
print(f"floor-odd      collisions (loose row hitting H(AP) or H(GW)): {len(coll_fo)}")
for S, h in coll_fo[:25]: print(f"    {S}  H={h}  ({'=AP' if h==H_AP_fo else '=GW'})")
print()
print(f"reciprocal-stk collisions: {len(coll_rs)}")
for S, h in coll_rs[:25]: print(f"    {S}  H={h}  ({'=AP' if h==H_AP_rs else '=GW'})")

print()
print("VERDICT:")
print(f"  floor-odd      separation survives large bank: {len(coll_fo)==0}")
print(f"  reciprocal-stk separation survives large bank: {len(coll_rs)==0}")
