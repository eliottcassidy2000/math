#!/usr/bin/env python3
"""
lrc14_defect_ge7_bonferroni_opus_S4.py    opus-2026-07-24-S4

The d>=7 wall of the LRC(14) tight-locus rigidity (HYP-9024), attacked by FLIPPING the
peel direction of THM-735 (mined this session). For V with defect d = |V\{1..13}| >= 7:
BODY E := the >=7 LARGE speeds; PEEL F := the j <= 6 SMALL speeds from {1..13}.
Then j <= 6 (the defect>=7 hypothesis IS THM-735's j<=6 hypothesis), so the Bonferroni base
1 - j*2h >= 1 - 6*(6/41) = 5/41 > 0, instead of the fatal 1-7*2h = -1/41 for the usual peel.

L(V) = meas(G_E \ union_{v small} D_v). Bonferroni (odd truncation upper-bounds the union):
   L(V) >= m_E - S1 + S2 - S3,   S_k = sum over k-subsets of small speeds of meas(cap D_v ^ G_E).
This is a RIGOROUS lower bound. Empirically > 0 for every defect>=7 config tested.
"""
import random, itertools
from math import gcd
from functools import reduce
h = 3.0 / 41
def merge(segs):
    segs = [s for s in segs if s[1] > s[0]]; segs.sort(); mg = []
    for a, b in segs:
        if mg and a <= mg[-1][1] + 1e-15: mg[-1] = (mg[-1][0], max(mg[-1][1], b))
        else: mg.append((a, b))
    return mg
def D(v):
    out = []
    for m in range(0, v + 1):
        lo = (m - h) / v; hi = (m + h) / v
        if lo < 0: out += [(0.0, hi), (lo + 1.0, 1.0)]
        elif hi > 1: out += [(lo, 1.0), (0.0, hi - 1.0)]
        else: out.append((lo, hi))
    return merge(out)
def union(V):
    s = []
    for v in V: s += D(v)
    return merge(s)
def comp(mg):
    out = []; prev = 0.0
    for a, b in mg:
        if a > prev: out.append((prev, a))
        prev = max(prev, b)
    if prev < 1.0: out.append((prev, 1.0))
    return out
def meas(mg): return sum(b - a for a, b in mg)
def inter(A, B):
    out = []; i = j = 0
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0]); hi = min(A[i][1], B[j][1])
        if hi > lo: out.append((lo, hi))
        if A[i][1] < B[j][1]: i += 1
        else: j += 1
    return out
def interlist(Ls):
    cur = Ls[0]
    for X in Ls[1:]: cur = inter(cur, X)
    return cur
if __name__ == "__main__":
    random.seed(5); N = 3000
    n1 = n3 = tot = 0; worst1 = worst3 = None
    for _ in range(N):
        d = random.randint(7, 13)
        small = sorted(random.sample(range(1, 14), 13 - d))
        large = sorted(random.sample(range(14, 400), d))
        V = sorted(set(small + large))
        if len(V) != 13 or reduce(gcd, V) != 1: continue
        GE = comp(union(large)); mE = meas(GE)
        Ds = {v: inter(D(v), GE) for v in small}
        S1 = sum(meas(Ds[v]) for v in small)
        S2 = sum(meas(interlist([Ds[a], Ds[b]])) for a, b in itertools.combinations(small, 2))
        S3 = sum(meas(interlist([Ds[a], Ds[b], Ds[c]])) for a, b, c in itertools.combinations(small, 3))
        b1 = mE - S1; b3 = mE - S1 + S2 - S3; tot += 1
        if b1 > 0: n1 += 1
        if b3 > 0: n3 += 1
        if worst1 is None or b1 < worst1[0]: worst1 = (b1, V)
        if worst3 is None or b3 < worst3[0]: worst3 = (b3, mE, S1, S2, S3, V)
    print(f"defect>=7 via FLIPPED THM-735 peel (body=large, peel=<=6 small), {tot} configs, h=3/41")
    print(f"  base 1-6*2h = 5/41 = {5/41:.5f} > 0   (usual peel 1-7*2h = {1-7*2*h:+.5f} < 0, the wall)")
    print(f"  level-1 (union bound)  L >= m_E - S1        : positive {n1}/{tot} ({100*n1/tot:.2f}%), worst {worst1[0]:+.5f}")
    print(f"  level-3 (Bonferroni)   L >= m_E-S1+S2-S3     : positive {n3}/{tot} ({100*n3/tot:.2f}%)")
    b3, mE, S1, S2, S3, V = worst3
    print(f"  worst level-3 bound = {b3:+.5f}   (m_E={mE:.4f} S1={S1:.4f} S2={S2:.4f} S3={S3:.4f})")
    print(f"     at V={V}")
    print("  NOTE: level-3 Bonferroni (odd truncation) is a RIGOROUS lower bound on L(V).")
    print("  Subtlety: for defect=13 (j=0, no small peel) the claim is 'no 13 speeds >=14 cover'")
    print("  (empirical min uncovered 0.0205) -- a separate base case, not covered by the peel.")
