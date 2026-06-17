#!/usr/bin/env python3
"""
LRC(14) PROVE — SPREAD THRESHOLD (family c), made rigorous.

CLASSICAL SUFFICIENT CONDITION (well known, here specialized to n=14):
  If 13 speeds satisfy a strong lacunarity / largeness condition, M(S) >= 1/14
  follows by an EXPLICIT lonely point, no search needed.

Two clean, PROVABLE sub-families:

(I) "BIG-AND-COPRIME-TO-LATTICE": pick tau = a/q for a denominator q chosen so
    that q | none of the (v_i +/- v_j) issues. Simplest: tau rational a/N.

(II) GREEDY NESTED INTERVALS (the standard lacunary proof). Process speeds in
    increasing order; maintain an interval I_k of tau-values where v_1..v_k are
    all lonely (dist >= 1/14). The set {tau: ||v tau|| >= 1/14} is a union of
    arcs of total measure 1 - 2/14 = 6/7 = 12/14. If at each step the current
    interval is long enough that v_{k+1} (a higher frequency) is guaranteed to
    have a lonely sub-arc inside it, we can keep a nonempty nested interval to
    the end => a lonely point exists => M >= 1/14.

    Sufficient step condition: an interval of length L contains a full period of
    v_{k+1} (length 1/v_{k+1}) whenever L >= 1/v_{k+1}; a full period always
    contains a lonely sub-arc of length (12/14)/v_{k+1}. So if we shrink to that
    sub-arc, the new length is (6/7)/v_{k+1}. For step k+1 to proceed we need
    (6/7)/v_k >= 1/v_{k+1}, i.e.  v_{k+1} >= (7/6) v_k.

    => THEOREM (specialized): If v_1 < v_2 < ... < v_13 with v_{k+1} >= (7/6) v_k
       for all k, AND the initial interval (length 6/7 around the v_1-lonely arc)
       is available, then M(S) >= 1/14. Ratio 7/6 ~ 1.1667.

We VERIFY this empirically: every S meeting the ratio condition has M >= 1/14,
and we find the true threshold ratio (smallest rho such that ALL S with ratio
>= rho are safe), to quantify how much of the family is provably safe.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def g(S, t): return min(nrm(v * t) for v in S)
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
def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at
def gcd_all(S): return reduce(gcd, S, 0)

THRESH = F(1, 14)

def min_ratio(S):
    S = sorted(S)
    return min(F(S[i+1], S[i]) for i in range(len(S)-1))

print("=== Verify ratio condition v_{k+1} >= (7/6) v_k => M >= 1/14 ===")
random.seed(11)
violations = 0
tested = 0
for _ in range(8000):
    S = [random.randint(1, 5)]
    ok = True
    for _ in range(12):
        # ensure ratio >= 7/6
        lo = (F(7, 6) * S[-1]).__ceil__()
        nxt = lo + random.randint(0, 8)
        S.append(nxt)
    S = sorted(set(S))
    if len(S) != 13: continue
    if gcd_all(S) != 1: continue
    if min_ratio(S) < F(7, 6): continue
    tested += 1
    m, _ = M(S)
    if m < THRESH:
        violations += 1
        print("  VIOLATION:", S, m, "ratio=", min_ratio(S))
print(f"  tested {tested} sets with min-ratio>=7/6; violations of M>=1/14: {violations}")

print()
print("=== Empirical TRUE threshold: scan ratio condition, find worst M per band ===")
random.seed(13)
bands = {}  # band index -> worst M
import math
for _ in range(40000):
    S = sorted(random.sample(range(1, 80), 13))
    if gcd_all(S) != 1: continue
    r = min_ratio(S)
    if r < 1: continue
    m, _ = M(S)
    # band by floor of 100*ratio
    bi = int(r * 20)  # 0.05 granularity
    cur = bands.get(bi)
    if cur is None or m < cur[0]:
        bands[bi] = (m, S, r)
print("  ratio_band(>=)   worstM     float    safe(>=1/14)?  example")
for bi in sorted(bands):
    m, S, r = bands[bi]
    lo = F(bi, 20)
    print(f"  >={float(lo):.2f}        {str(m):>9}  {float(m):.5f}   "
          f"{'YES' if m>=THRESH else 'NO '}        ratio={float(r):.3f}")
