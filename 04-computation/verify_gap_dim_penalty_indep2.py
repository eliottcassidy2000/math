#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Part 2: exhaustive search + single-defect family + spread-GAP penalty."""
from fractions import Fraction as F
import itertools
from math import gcd
from functools import reduce

def N_at(E, x):
    hit = set()
    for e in E:
        v = e * x; v = v - (v.numerator // v.denominator)
        hit.add((v.numerator * 7) // v.denominator)
    return sum(1 for j in range(1, 7) if j not in hit)

def dist_p(E):
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for a in range(0, 7 * e + 1): bps.add(F(a, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [F(0)] * 7
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i + 1]
        if hi == lo: continue
        p[N_at(E, (lo + hi) / 2)] += (hi - lo)
    return p

def g_poly(k):
    g = []
    for t in range(7):
        if k == 8: g.append(F((t-1)*(t-2)*(t-4)*(t-5), 40))
        else: g.append(F(-(t-2)*(t-3)*(t-6), 36))
    return g

def L_y(E, k):
    p = dist_p(E); g = g_poly(k)
    return sum(p[t] * g[t] for t in range(7))

def primitive(E):
    return reduce(gcd, sorted(set(E))) == 1

def is_AP(E):
    E = sorted(set(E))
    if len(E) <= 2: return True
    d = E[1] - E[0]
    return all(E[i+1]-E[i] == d for i in range(len(E)-1))

# consec_9 p_0 check
p = dist_p(list(range(9)))
print(f"consec_9 p_0 = {p[0]} = {float(p[0]):.6f}  (claim 2447/5880)")
print()

# ---- Exhaustive over ALL primitive non-AP 9-sets with max<=18, 0 in E ----
k = 9
best = None; cnt = 0
for combo in itertools.combinations(range(1, 19), k - 1):
    E = (0,) + combo
    if not primitive(E): continue
    if is_AP(E): continue
    cnt += 1
    L = L_y(E, k)
    if best is None or L > best[0]:
        best = (L, E)
print(f"k=9 exhaustive non-AP, max<=18, 0 in E: scanned {cnt} primitive non-AP sets")
print(f"  global max L_y = {best[0]} = {float(best[0]):.6f}  at E={list(best[1])}")
print(f"  (claim: unique max 38681/79380 at [0,1,2,3,4,5,6,7,9])")
print()

# ---- single-defect family supremum: E = {0..L-1} U {L-1+s}, s>=2 ----
print("single-defect family E={0..L-1} U {L-1+s}, k=8,9,10; find sup over s,L:")
for k in [8, 9, 10]:
    gk = (lambda t: F((t-1)*(t-2)*(t-4)*(t-5),40)) if k==8 else (lambda t: F(-(t-2)*(t-3)*(t-6),36))
    g = [gk(t) for t in range(7)]
    def Ly(E):
        p = dist_p(E); return sum(p[t]*g[t] for t in range(7))
    bestf = None
    for L in range(1, k):
        base = list(range(L))
        # remaining k-L points after a gap of size s, then consecutive
        for s in range(2, 40):
            start = L - 1 + s
            E = base + list(range(start, start + (k - L)))
            if len(set(E)) != k: continue
            if not primitive(E) or is_AP(E): continue
            val = Ly(E)
            if bestf is None or val > bestf[0]:
                bestf = (val, E, s, L)
    print(f"  k={k}: sup L_y = {bestf[0]} = {float(bestf[0]):.6f} at E={bestf[1]} (s={bestf[2]},L={bestf[3]})")
print()

# ---- a spread 2-dim GAP for contrast (k=8) ----
E = [0,3,5,6,8,9,11,14]
print(f"spread 2-dim GAP k=8 {E}: L_y={float(L_y(E,8)):.6f} (claim 0.22453)")
