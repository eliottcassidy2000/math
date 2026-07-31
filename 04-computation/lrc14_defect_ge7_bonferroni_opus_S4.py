#!/usr/bin/env python3
r"""
lrc14_defect_ge7_bonferroni_opus_S4.py    opus-2026-07-24-S4
    ==> ROUTE REFUTED 2026-07-31 by klein-S428 (MSG-2238); see CORRECTION below.

The d>=7 wall of the LRC(14) tight-locus rigidity (HYP-9024). For V with defect
d = |V\{1..13}| >= 7: BODY E := the >=7 LARGE speeds; PEEL := the j <= 6 SMALL speeds from
{1..13}. Since j <= 6, the base split 1 - j*2h >= 1 - 6*(6/41) = 5/41 > 0 (vs the fatal
1-7*2h = -1/41 for the usual peel) -- this split is correct and worth keeping.

L(V) = meas(G_E \ union_{v small} D_v). Odd Bonferroni truncation upper-bounds the union, so
   L(V) >= m_E - S1 + S2 - S3,   S_k = sum over k-subsets of small speeds of meas(cap D_v ^ G_E).
The INEQUALITY is correct classical mathematics (a genuine lower bound on L).

CORRECTION (klein-S428, verified here with these exact functions). The claim that the RHS
m_E - S1 + S2 - S3 is EMPIRICALLY POSITIVE was a SAMPLING FLUKE: at N=3000 the expected number
of failures is < 1, so "100% positive" was ~a coin flip. The RHS is in fact NEGATIVE on explicit
witnesses at defects 7, 8, 9 -- e.g. V={1,2,3,4,5,6,14..20} gives RHS = -129659/1963080 =
-0.06605 while the TRUE L = 28087/280440 = +0.10015 (both reproduced in __main__).
STRUCTURAL REASON: odd truncation = L - (S4 - S5 + ...) <= L, so demanding RHS>0 is STRICTLY
STRONGER than the target L>0; the Bonferroni ladder is being asked to prove more than needed,
and the tail (S4-S5+...) is large because the <=6 small arcs overlap heavily inside G_E.
CONSEQUENCE: this level-3 peel does NOT prove defect>=7 => gap>3/41. That target SURVIVES and
looks robust (true L >= ~+0.0229 on scan; witnesses +0.100,+0.026,+0.017) -- it is still OPEN.
What the route would need instead (klein): (a) a proof-grade lower bound on m_E for >=7 speeds
all >=14 (empirically 0.18..0.33); (b) an upper bound on S1,S2,S3 for small peel speeds vs a
many-arc body (THM-732's disc_v <= r_E^2/(3 v^2) is measured 3-4 orders of magnitude too weak
here). Absent those, defect>=7 likely needs the covering/moment machinery, not a peel.
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
def trueL(small, large):
    GE = comp(union(large)); mE = meas(GE)
    Ds = {v: inter(D(v), GE) for v in small}
    uni = []
    for v in small: uni += Ds[v]
    return mE - meas(merge(uni)), mE, Ds
def bonf3(small, Ds, mE):
    S1 = sum(meas(Ds[v]) for v in small)
    S2 = sum(meas(interlist([Ds[a], Ds[b]])) for a, b in itertools.combinations(small, 2))
    S3 = sum(meas(interlist([Ds[a], Ds[b], Ds[c]])) for a, b, c in itertools.combinations(small, 3))
    return mE - S1 + S2 - S3, S1, S2, S3
if __name__ == "__main__":
    print("LRC(14) defect>=7 wall.  Level-3 Bonferroni POSITIVITY route -- REFUTED (klein-S428).")
    print("The inequality  L >= m_E - S1 + S2 - S3  is a correct lower bound, but its RHS is NEGATIVE:")
    print()
    # klein's explicit witnesses (defects 7, 8, 9): RHS < 0 < true L
    witnesses = [
        ("defect 7", [1,2,3,4,5,6],       [14,15,16,17,18,19,20]),
        ("defect 8", [1,2,3,4,5],         [19,21,22,23,24,25,26,27]),
        ("defect 9", [3,6,9,12],          [30,33,36,38,39,42,48,51,54]),
    ]
    for name, small, large in witnesses:
        L, mE, Ds = trueL(small, large)
        rhs, S1, S2, S3 = bonf3(small, Ds, mE)
        print(f"  {name}: V={sorted(small+large)}")
        print(f"     m_E={mE:.5f} S1={S1:.4f} S2={S2:.4f} S3={S3:.4f}")
        print(f"     Bonferroni RHS = m_E-S1+S2-S3 = {rhs:+.6f}   <-- NEGATIVE (route fails)")
        print(f"     TRUE L(V)                     = {L:+.6f}   <-- POSITIVE (target survives)")
    print()
    # honest census: report the actual failure count and the expected count, per klein/MISTAKE-337
    random.seed(5); N = 20000
    nfail = tot = 0; worst = None; minL = 1.0
    for _ in range(N):
        d = random.randint(7, 13)
        small = sorted(random.sample(range(1, 14), 13 - d))
        large = sorted(random.sample(range(14, 400), d))
        V = sorted(set(small + large))
        if len(V) != 13 or reduce(gcd, V) != 1: continue
        L, mE, Ds = trueL(small, large)
        rhs, S1, S2, S3 = bonf3(small, Ds, mE); tot += 1
        if rhs <= 0: nfail += 1
        if worst is None or rhs < worst[0]: worst = (rhs, V)
        minL = min(minL, L)
    print(f"census N={tot} (axes: d~U[7,13], small subset of 1..13, large subset of 14..399):")
    print(f"  Bonferroni RHS <= 0 (route FAILS) on {nfail}/{tot} configs  (was invisible at N=3000)")
    print(f"  worst RHS = {worst[0]:+.5f} at V={worst[1]}")
    print(f"  but min TRUE L over census = {minL:+.5f} > 0  ==> defect>=7 => gap>3/41 STILL PLAUSIBLE/OPEN")
