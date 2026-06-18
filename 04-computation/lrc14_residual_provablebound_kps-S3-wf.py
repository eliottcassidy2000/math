#!/usr/bin/env python3
"""
LRC(14) S3 -- a PROVABLE global-witness bound + bounded-Vmax finite reduction.

Two rigorous routes tested:

ROUTE 1 (measure global witness, the clean theorem):
  Let P = small part (subset of {1..13}, covering most q's), L = cluster (all >13),
  Vmin=min L, Vmax=max L, c=|L|, s=Vmax-Vmin.
  The small part safe measure: mu_P = meas{tau: ||u tau||>=1/14 for all u in P}.
  The cluster danger measure over the WHOLE circle: each u in L blocks measure
  u*(1/(7u)) = 1/7. Union over L blocks <= c/7. BUT overlaps: the cluster danger
  has measure exactly meas(union of L-teeth). If mu_P > meas(L-danger restricted to
  P-safe), a global witness exists. Provable form:
      mu_P - (cluster danger measure) > 0  ==> global witness EXISTS (just positivity).
  Refined: we need a cluster-free POINT in P-safe of positive measure. Sufficient:
      mu_P > c/7   (crude)  OR  mu_P > meas(L-teeth)  (exact).
  Test how often mu_P > meas(L-teeth ∩ Psafe). This gives M(S)>=1/14 by MEASURE
  (a point exists) -- but NOT width > 1/(7Vmax) (no robustness margin). For LRC we
  only need M(S)>=1/14, i.e. ONE safe point. So positive-measure global witness SUFFICES!

  *** This is the key realization: for M(S)>=1/14 we need a single tau with min>=1/14,
      i.e. the GLOBAL safe set G_S has positive measure (it is closed, so positive
      measure => nonempty, in fact contains an arc). We do NOT need width>1/(7Vmax).***

  So the clean lemma: meas(G_S) >= mu_P - meas(L-danger) and if this is >0 we are done.
  Test: is mu_P > meas(L-danger) on all S3? (meas(L-danger) <= c/7.)

ROUTE 2 (bounded Vmax finite reduction):
  Find smallest threshold T(P) such that Vmax >= T(P) => global witness guaranteed by a
  provable inequality; then S3 with Vmax < T(P) is a FINITE set to check exactly.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random
from collections import Counter

def circle_meas_union(ivs):
    """measure of union of intervals (a,b) on [0,1), a,b in [0,1), possibly wrapping
    handled by caller (we pass non-wrapping here)."""
    ivs = sorted(ivs); m = []
    for a, b in ivs:
        if m and a <= m[-1][1]: m[-1] = (m[-1][0], max(m[-1][1], b))
        else: m.append((a, b))
    return sum(b - a for a, b in m), m

def teeth(u, h=F(1, 14)):
    """all teeth of runner u on [0,1) as (a,b) with possible wrap split."""
    out = []
    for j in range(u):
        c = F(j, u); a = (c - h/u) % 1; b = (c + h/u) % 1
        if a < b: out.append((a, b))
        else: out.append((a, F(1))); out.append((F(0), b))
    return out

def danger_union(A):
    ivs = []
    for u in A: ivs.extend(teeth(u))
    meas, m = circle_meas_union(ivs)
    return meas, m

def safe_meas_and_arcs(A):
    """return (measure of safe set, list of safe arcs (a,b))."""
    dmeas, dm = danger_union(A)
    safe = []; prev = F(0)
    for a, b in dm:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    smeas = sum(b - a for a, b in safe)
    return smeas, safe

def global_safe_meas(S):
    smeas, safe = safe_meas_and_arcs(S)
    return smeas, safe

def gcd_all(S): return reduce(gcd, S, 0)
def is_covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def classify(S):
    S = sorted(set(S)); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13 * Vmin: return 'S2'
    return 'S3'

def gen_S3(seed=0, target=2000, vmax_cap=3000):
    rng = random.Random(seed); out = []; tries = 0
    while len(out) < target and tries < target * 500:
        tries += 1
        c = rng.choice([2, 2, 3, 3, 4, 5, 6, 7, 8]); nsmall = 13 - c
        if nsmall < 1: continue
        V = rng.choice([rng.randint(40, 80), rng.randint(80, 200),
                        rng.randint(200, 600), rng.randint(600, vmax_cap)])
        spread = rng.choice([14, 20, 28, 35, 42, 56, 70, 90])
        cluster = set()
        while len(cluster) < c: cluster.add(V + rng.randint(0, spread))
        small = set(); pool = list(range(1, 14)); rng.shuffle(pool)
        for x in pool:
            if len(small) >= nsmall: break
            small.add(x)
        if len(small) < nsmall: continue
        S = sorted(small | cluster)
        if len(S) != 13 or gcd_all(S) != 1: continue
        if not is_covering(S) or classify(S) != 'S3': continue
        out.append(S)
    return out

if __name__ == "__main__":
    S3 = gen_S3(seed=11, target=800)
    print(f"n S3 = {len(S3)}")
    print("\nROUTE 1: positive-measure GLOBAL safe set G_S (one safe point => M>=1/14).")
    nz = 0; mins = (F(10), None)
    muP_vs_Ldanger = 0
    crude = 0
    margins = []
    for S in S3:
        small = [u for u in S if u <= 13]; cluster = [u for u in S if u > 13]
        c = len(cluster)
        gm, garcs = global_safe_meas(S)
        margins.append(float(gm))
        if gm > 0: nz += 1
        if gm < mins[0]: mins = (gm, S)
        # crude test: mu_P > c/7 ?
        muP, _ = safe_meas_and_arcs(small)
        Ld, _ = danger_union(cluster)
        if muP > Ld: muP_vs_Ldanger += 1
        if muP > F(c, 7): crude += 1
    n = len(S3)
    def st(a): return f"min={min(a):.4f} mean={sum(a)/len(a):.4f} max={max(a):.4f}"
    print(f"  global safe measure meas(G_S): {st(margins)}")
    print(f"  meas(G_S) > 0 (=> M>=1/14) on {nz}/{n}")
    print(f"  min meas(G_S) = {float(mins[0]):.5f} at S={mins[1]}")
    print(f"  crude bound  mu_P > meas(L-danger): holds {muP_vs_Ldanger}/{n}")
    print(f"  cruder bound mu_P > c/7           : holds {crude}/{n}")
    # also: distribution of c among S3 with the crude bound failing
    print("\n  (mu_P = small-part safe measure; meas(L-danger) <= c/7)")
