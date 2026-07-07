#!/usr/bin/env python3
"""
kps-2026-07-06-S52 -- THE COARSE / SCALE REDUCTION, AIMED DIRECTLY AT THE SUP (Route 1).

Post opus-S130: Route 2 (J-K -> rank-2 -> (C) 1-D gap) is RETIRED -- J-K bounds
ACCUMULATION POINTS, not the SUPREMUM (MISTAKE-117), so closing (C) would not close
LRC(14).  The correctly-aimed route is ROUTE 1: bound M(v) >= 1/14 DIRECTLY.

The coarse reduction is a GENERAL fact about the sup M, so it survives the Route-2
collapse and applies DIRECTLY to LRC(14) via the SETTLED LRC(<=13).  No J-K.

LRC(14): 13 nonzero speeds v_1..v_13, want M(v) = sup_t min_i ||v_i t|| >= 1/14.

COARSE BOUND.  If v_i = a_i + L*k_i with all k_i >= 1 and |a_i| <= A, and K = {distinct
k_i}, then at the witness t = t_K / L (t_K = K's optimal witness):
    ||v_i t|| = ||a_i t_K/L + k_i t_K|| >= ||k_i t_K|| - |a_i|/L >= M(K) - A/L.
So  M(v) >= M(K) - A/L.   [RIGOROUS -- scaled witness, needs k_i >= 1]

DICHOTOMY on r = #distinct k_i (= #clusters at scale L):
  * r <= 12  (two speeds share a cluster = a close pair at scale L):
        K has <= 12 speeds => M(K) >= 1/13 by LRC(<=13, SETTLED)
        => M(v) >= 1/13 - A/L > 1/14  for L > 182*A.        LONELY -- no circularity.
  * r = 13   (all clusters distinct): K is a 13-speed family of SMALLER height
        (height(K) ~ height(v)/L).  DESCEND.  Terminates at a single-scale family.

RESIDUE (the honest hard core): SINGLE-SCALE / COMPRESSED families -- no scale gap,
so the coarse reduction is vacuous.  These are exactly where decorrelation (Route 1's
Fourier/witness-density floor) is the genuine open analytic core.

This script verifies the coarse bound and the r<=12 => M(v)>1/14 conclusion on the
CORRECT object: 13-speed families, threshold 1/14.
"""
from fractions import Fraction
import numpy as np
from math import gcd
from functools import reduce
import random

ONE_14 = Fraction(1, 14)
ONE_13 = Fraction(1, 13)

def Mw(v):
    """Exact M(v) = max over q of (max over c of min_i dist(v_i c, qZ)) / q."""
    v = [x for x in v if x]
    if not v:
        return Fraction(0)
    S = sum(abs(x) for x in v)
    Q = min(4 * S, 2 * max(abs(x) for x in v) + 2)
    va = np.array(v, dtype=np.int64)
    bn, bd = 0, 1
    for q in range(2, Q + 1):
        a = np.arange(1, q)
        r = np.outer(va, a) % q
        d = np.minimum(r, q - r)
        col = d.min(axis=0)
        qb = int(col.max())
        if qb * bd > bn * q:
            bn, bd = qb, q
    return Fraction(bn, bd)

def reduce_gcd(v):
    g = reduce(gcd, v)
    return [x // g for x in v] if g > 1 else list(v)

print("=" * 78)
print("COARSE / SCALE REDUCTION for LRC(14), aimed DIRECTLY at the sup (Route 1)")
print("Target: M(v) >= 1/14 = %.5f for 13 nonzero speeds.  LRC(<=13): M >= 1/13 = %.5f" %
      (float(ONE_14), float(ONE_13)))
print("=" * 78)

random.seed(7)

print("\n(1) r <= 12 KILLER: <= 12 clusters at scale L => M(v) >= 1/13 - A/L > 1/14.")
print("    (grounds in SETTLED LRC(<=13); the r<=12 branch needs NO new analysis)")
print("    v = {a_i + L*k_i}, 13 speeds, k_i in a set of size r<=12, |a_i| <= A.")
hdr = "    %-4s %-3s %-8s %-10s %-10s %-9s %s"
print(hdr % ("L", "r", "M(K)", "1/13-A/L", "M(v)", ">1/14?", "k-multiset (clusters)"))
for L in (60, 120, 240):
    for trial in range(3):
        # 13 speeds, cluster labels k_i drawn from {1..r} with r<=12 (=> a repeat)
        r_target = random.randint(2, 12)
        ks = [random.randint(1, r_target) for _ in range(13)]
        # ensure at least one repeat (r<=12 guaranteed since 13 labels in <=12 values)
        A = 12
        a = [random.randint(-A, A) for _ in range(13)]
        v = sorted(set(a[i] + L * ks[i] for i in range(13)))
        if len(v) < 13:  # collision after dedup -- skip (still valid family, but keep 13 distinct)
            continue
        v = reduce_gcd(v)
        K = sorted(set(ks))
        MK = Mw(list(K))
        Mv = Mw(v)
        floor = ONE_13 - Fraction(A, L)
        ok = Mv > ONE_14
        print(hdr % (L, len(K), str(MK), "%.5f" % float(floor),
                     "%.5f" % float(Mv), "YES" if ok else "no!!",
                     str(sorted(ks))))

print("\n(2) COARSE BOUND M(v) >= M(K) - A/L holds (verify the inequality directly):")
print("    r can be up to 13; bound is M(v) >= M(K) - A/L regardless.")
viol = 0
for L in (100, 200):
    for trial in range(6):
        r_target = random.randint(2, 13)
        ks = [random.randint(1, r_target) for _ in range(13)]
        A = 10
        a = [random.randint(-A, A) for _ in range(13)]
        v = sorted(set(a[i] + L * ks[i] for i in range(13)))
        if len(v) < 13:
            continue
        K = sorted(set(ks))
        MK = Mw(list(K))
        Mv = Mw(reduce_gcd(v))
        bound = MK - Fraction(A, L)
        holds = Mv >= bound - Fraction(1, 10**6)  # tiny slack for the gcd-reduction shift
        if not holds:
            viol += 1
        print("    L=%-4d r=%-2d  M(K)=%-6s  M(K)-A/L=%.5f   M(v)=%.5f   bound holds? %s"
              % (L, len(K), str(MK), float(bound), float(Mv), "YES" if holds else "NO"))
print("    coarse-bound violations: %d (expect 0)" % viol)

print("\n(3) r=13 DESCENT: all clusters distinct => K is a SMALLER 13-family, M(v)~M(K).")
print("    height(K) = max k_i << height(v) ~ L*max k_i.  Recurse -> single-scale base.")
for L in (200,):
    ks = random.sample(range(1, 40), 13)          # 13 DISTINCT cluster labels
    a = [random.randint(-8, 8) for _ in range(13)]
    v = sorted(set(a[i] + L * ks[i] for i in range(13)))
    K = sorted(ks)
    print("    v height=%d  ->  K=%s height=%d   (ratio ~ L=%d)" %
          (max(v), K, max(K), L))
    print("    M(v)=%.5f   M(K)=%.5f   |M(v)-M(K)|=%.5f  (< A/L=%.5f)" %
          (float(Mw(reduce_gcd(v))), float(Mw(list(K))),
           abs(float(Mw(reduce_gcd(v))) - float(Mw(list(K)))), 8.0 / L))

print("\n(4) RESIDUE = single-scale / compressed families (no scale gap): the hard core.")
print("    The AP {1..13} is the extremal single-scale family: M = 1/14 EXACTLY (tight).")
ap = list(range(1, 14))
print("    M({1..13}) = %s = %.6f   (= 1/14, the LRC(14) threshold)" %
      (str(Mw(ap)), float(Mw(ap))))
print("    A single-scale family has k_i=v_i, a_i=0, L=1 => K=v: coarse bound is vacuous.")
print("    => these need DECORRELATION (Route 1 Fourier floor), NOT the coarse reduction.")

print("\n" + "=" * 78)
print("CONCLUSION.  The coarse reduction, on the CORRECT object (sup M, Route 1):")
print("  * clusters-into-<=12-groups families  =>  M(v) > 1/14  via SETTLED LRC(<=13).")
print("  * distinct-13-cluster families        =>  DESCEND (smaller height).")
print("  * single-scale / compressed families  =>  RESIDUE = the analytic core.")
print("  J-K FREE.  Does NOT close LRC(14); RE-GROUNDS the compression reduction on the")
print("  sup, and isolates the decorrelation core.  Survives opus-S130's Route-2 retire.")
print("=" * 78)
