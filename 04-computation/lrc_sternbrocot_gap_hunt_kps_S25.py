#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S25: HUNT the simplest gap members in the DIVISIBILITY-RICH
space at the scale the lever demands -- the underlying Stern-Brocot structure.

Two facts combine:
 (1) gap_candidate_has_multiple (S24, GREEN): a gap member has a multiple of every
     k in {2..12} (in particular of 5,7,8,9,11,12).
 (2) mac-mini's lever (HYP-4432): M=c/q => q | (v_i +- v_j) or 2v_i, so max >= q/2.

The gap (1/13, 2/25) is a FAREY interval (det 1*25-2*13 = -1); its values are the
STERN-BROCOT descendants of the pair, simplest = the mediant 3/38 (q=38=13+25).
So a gap member needs q >= 38 => max >= 19.  My S21 census (max<=18) COULD NOT have
found any -- it was structurally blind to gap members.  This targets them: the
divisibility-rich families with max >= 19, hunting M in (1/13, 2/25).
"""
from fractions import Fraction
import numpy as np, random
from math import gcd
from functools import reduce

LO, HI = Fraction(1, 13), Fraction(2, 25)

def M_exact(v):
    S = int(sum(abs(x) for x in v)); Q = min(4*S, 2*max(abs(x) for x in v)+2)
    va = np.array(v, dtype=np.int64); bn, bd = 0, 1
    for q in range(2, Q+1):
        a = np.arange(1, q, dtype=np.int64); r = np.outer(va, a) % q
        d = np.minimum(r, q-r); bq = int(d.min(axis=0).max())
        if bq*bd > bn*q: bn, bd = bq, q
    return Fraction(bn, bd)

# ---- (1) the Stern-Brocot structure of the gap denominators ----
print("=== the gap (1/13, 2/25) as a Farey interval: Stern-Brocot descendants ===", flush=True)
print(f"  det check: 1*25 - 2*13 = {1*25 - 2*13} (Farey neighbors)", flush=True)
def factor(n):
    f = []; d = 2; m = n
    while d*d <= m:
        while m % d == 0: f.append(d); m //= d
        d += 1
    if m > 1: f.append(m)
    return f
# gap fractions c/q in lowest terms with 1/13 < c/q < 2/25, by denominator
gapfracs = []
for q in range(14, 130):
    for c in range(1, q):
        fr = Fraction(c, q)
        if LO < fr < HI and fr.denominator == q:
            gapfracs.append((q, c, fr))
gapfracs.sort()
print("  simplest gap fractions (denominator, value, factorization of q):", flush=True)
for q, c, fr in gapfracs[:12]:
    print(f"    {c}/{q} = {float(fr):.5f}   q = {q} = {'*'.join(map(str,factor(q)))}", flush=True)
print(f"  => mediant 3/38 (q=38=13+25) is the simplest; gap members need q>=38, so max>=19.", flush=True)

# ---- (2) targeted hunt: divisibility-rich families, max in [19, MAXH] ----
print(flush=True)
print("=== HUNT: divisibility-rich families with max >= 19, looking for M in (1/13,2/25) ===", flush=True)
random.seed(20260706)
def mult_of(k, lo, hi):
    lst = [x for x in range(k, hi+1, k) if x >= 1]
    return lst
def build_divrich(MAXH):
    """build a 12-family guaranteed to have a multiple of each k in {5,7,8,9,11,12}, distinct, <= MAXH."""
    need = [5,7,8,9,11,12]
    chosen = set()
    for k in need:
        opts = [x for x in range(k, MAXH+1, k) if x not in chosen]
        if not opts: return None
        chosen.add(random.choice(opts))
    # fill to 12 with random distinct speeds <= MAXH
    pool = [x for x in range(1, MAXH+1) if x not in chosen]
    random.shuffle(pool)
    for x in pool:
        if len(chosen) >= 12: break
        chosen.add(x)
    if len(chosen) != 12: return None
    v = sorted(chosen)
    if reduce(gcd, v) != 1: return None
    return v

hits = []; near = []; n = 0
for MAXH in [24, 28, 32, 38, 45]:
    trials = 40000
    ingap = 0; nn = 0
    for _ in range(trials):
        v = build_divrich(MAXH)
        if v is None or max(v) < 19: continue
        nn += 1; n += 1
        M = M_exact(v)
        if LO < M < HI:
            ingap += 1; hits.append((v, M))
            if len(hits) <= 5: print(f"    *** GAP MEMBER: v={v} M={M}={float(M):.6f} ***", flush=True)
        elif HI <= M <= Fraction(1,12):
            near.append((v, float(M)))
    print(f"  MAXH={MAXH}: {nn} divisibility-rich families (max>=19); gap members = {ingap}", flush=True)

print(flush=True)
print(f"TOTAL: {n} divisibility-rich families probed (max in [19,45]); gap members = {len(hits)}", flush=True)
if not hits:
    print("  NO gap member in the divisibility-rich space at max in [19,45].", flush=True)
    print("  Combined with S21 (max<=18 forced loose/tight), the gap is empty over max<=45", flush=True)
    print("  in the ONLY space where a member could live (divisibility-rich).", flush=True)
# show the tightest near-2/25 divisibility-rich families (the real boundary)
near.sort(key=lambda x: x[1])
print(f"  tightest divisibility-rich families just above 2/25:", flush=True)
for v, m in near[:8]:
    print(f"    v={v}  M={m:.6f}", flush=True)
