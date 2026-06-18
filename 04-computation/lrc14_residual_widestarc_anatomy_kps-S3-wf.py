#!/usr/bin/env python3
"""
LRC(14) S3 residual -- ANATOMY of the widest arc of S\\{max}.

The arc-width lemma needs W(S\\{v}) > 1/(7v) for SOME v. Empirically v=max works.
WHY is the widest safe arc of S\\{max} wide? This script dissects the widest arc:
  - which runner's tooth bounds it on the LEFT, which on the RIGHT
  - the arc's location tau, width, and width*7*max (the criterion margin)
  - whether the bounding teeth come from the cluster, the small part, or mixed

GOAL: find an EXACT, provable lower bound on W(S\\{max}) in terms of the cluster
structure, to convert the empirical C(S) into a theorem on a sub-regime.

KEY IDEA being tested (multi-band-fit): the widest arc lives between two ADJACENT
teeth of the merged danger set. If we can show that SOMEWHERE there is a pair of
adjacent teeth (left edge L, right edge R) with R-L > 1/(7*max), and that the gap
(L,R) is safe for ALL runners (it is, by construction of the merged safe set), we
are done. So we need: the largest gap between consecutive merged-tooth-edges
exceeds 1/(7*max). This is a STATEMENT ABOUT THE MERGED TOOTH SET of S\\{max}.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random
from collections import Counter

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def teeth_edges(A, h=F(1, 14)):
    """Return sorted list of (edge, kind, runner) for all danger-tooth boundaries.
    Each runner u has teeth centered at j/u (j=0..u-1) of half-width h/u.
    We record left edges (j/u - h/u) and right edges (j/u + h/u), mod 1."""
    pts = []
    for u in A:
        for j in range(0, u):
            c = F(j, u)
            a = (c - h/u) % 1; b = (c + h/u) % 1
            pts.append((a, 'L', u)); pts.append((b, 'R', u))
    return pts

def safe_components(A, h=F(1, 14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h/u) % 1; b = (c + h/u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    safe = []; prev = F(0)
    for a, b in merged:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def widest_arc_full(A):
    """Return (width, lo, hi, left_runner, right_runner) for widest safe arc.
    left_runner = runner whose tooth right-edge equals lo (bounds arc on left).
    right_runner = runner whose tooth left-edge equals hi."""
    sc = safe_components(A)
    if not sc: return None
    # build map from edge value to runner(s)
    redge = {}  # right-edge value -> runners (these bound a safe arc on its LEFT)
    ledge = {}  # left-edge value  -> runners (bound a safe arc on its RIGHT)
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - F(1,14)/u) % 1; b = (c + F(1,14)/u) % 1
            redge.setdefault(b, set()).add(u)
            ledge.setdefault(a, set()).add(u)
    cand = []
    for a, b in sc:
        cand.append((b - a, a, b))
    # wrap arc
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1:
        cand.append((sc[0][1] + (1 - sc[-1][0]), sc[-1][0], sc[0][1] + 1))
    width, lo, hi = max(cand)
    lo_m = lo % 1
    hi_m = hi % 1
    lr = redge.get(lo_m, set())   # who bounds the arc on the left
    rr = ledge.get(hi_m, set())   # who bounds it on the right
    return (width, lo, hi, lr, rr)

def gcd_all(S): return reduce(gcd, S, 0)
def is_covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def classify(S):
    S = sorted(set(S)); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13 * Vmin: return 'S2'
    return 'S3'

def gen_S3(seed=0, target=2000):
    rng = random.Random(seed); out = []; tries = 0
    while len(out) < target and tries < target * 400:
        tries += 1
        c = rng.choice([2, 2, 3, 3, 4, 5, 6]); nsmall = 13 - c
        if nsmall < 1: continue
        V = rng.choice([rng.randint(40, 80), rng.randint(80, 200),
                        rng.randint(200, 600), rng.randint(600, 1500)])
        spread = rng.choice([14, 20, 28, 35, 42, 56, 70])
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
    S3 = gen_S3(seed=2, target=500)
    print(f"n S3 = {len(S3)}")
    # Anatomy of widest arc of S\\{max}
    bound_kinds = Counter()   # (left from cluster?, right from cluster?)
    margins = []
    left_is_small = 0; right_is_small = 0; both_small = 0; both_cluster = 0
    minmargin = (F(10), None)
    for S in S3:
        Vmax = S[-1]; rest = [u for u in S if u != Vmax]
        cluster = set(u for u in S if u > 13)
        res = widest_arc_full(rest)
        if res is None: continue
        width, lo, hi, lr, rr = res
        margin = width * 7 * Vmax
        margins.append(float(margin))
        if margin < minmargin[0]: minmargin = (margin, S, lr, rr, width)
        lc = any(u in cluster for u in lr); rc = any(u in cluster for u in rr)
        ls = any(u <= 13 for u in lr); rs = any(u <= 13 for u in rr)
        bound_kinds[(lc, rc)] += 1
        if ls and not lc: left_is_small += 1
        if rs and not rc: right_is_small += 1
    print(f"widest-arc margin (W*7*max) over S3: min={min(margins):.3f} "
          f"mean={sum(margins)/len(margins):.3f} max={max(margins):.3f}")
    print(f"min-margin set: margin={float(minmargin[0]):.4f}")
    print(f"   S={minmargin[1]}")
    print(f"   left-bound runners={sorted(minmargin[2])}, right-bound runners={sorted(minmargin[3])}, width={minmargin[4]}")
    print("\nWho bounds the widest arc of S\\max  (left-from-cluster, right-from-cluster):")
    for k, v in sorted(bound_kinds.items()):
        print(f"   left_cluster={k[0]}, right_cluster={k[1]}: {v}")
    print(f"left bounded purely by small (no cluster): {left_is_small}")
    print(f"right bounded purely by small (no cluster): {right_is_small}")
