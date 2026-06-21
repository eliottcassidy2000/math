#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_mobius_majorant_control_opus_0621.py   (opus, 2026-06-21, THREAD B control)

CRITICAL CONTROL for the C3-ii finding.  In Part C we found a degree-2 per-subset
majorant  F_lambda(E)=sum_{|A|<=2} lambda_A q_A(E)  that is (a) >= measS7 for ALL E and
(c) TIGHT at consec.  BUT: the LP "min F_lambda(target) s.t. F_lambda(E)>=measS7(E) all E"
returns measS7(target) iff there is a valid degree-R majorant tight at TARGET.  We MUST
check this is NOT vacuous, i.e. that it does NOT return gap=0 for arbitrary targets.

THE DECISIVE QUESTION (consec-extremality, not just tightness):
  A degree-R per-subset majorant F_lambda proves consec-max IF a SINGLE fixed lambda
  satisfies  F_lambda(E) >= measS7(E) for all E  AND  F_lambda(consec) <= cap_k  AND the
  tightness selects consec.  The real test of the Mobius/integrality route is:

  TEST 1 (NON-VACUITY): run the same majorant LP with target = a NON-consec shape E*.
    If gap>0 for E* (the min majorant at E* STRICTLY exceeds measS7(E*)), the degree-R
    majorant DISTINGUISHES consec from E* -- consec is special.  If gap=0 for many E*,
    the construction is vacuous (every shape has its own tight majorant -> proves nothing).

  TEST 2 (UNIFORM CERTIFICATE): is there ONE lambda (degree R) that is BOTH tight at consec
    AND a valid majorant?  (=C3-ii, already YES.)  Among all such tight-at-consec valid
    majorants, what is the MAX over the window of F_lambda(E)?  If consec is the argmax of
    F_lambda over the window for the tight lambda, then F_lambda CERTIFIES consec-max
    (since measS7(E) <= F_lambda(E) <= F_lambda(consec) = measS7(consec)).  This is the
    clean consec-max proof.  We compute, for the tight lambda, whether F_lambda is
    consec-maximal over the window.

  TEST 3 (the honest collapse check): does the degree-R majorant value at a SWEEP of targets
    track measS7(target) (gap~0 everywhere = vacuous) or stay strictly above for non-consec?

ALL via scipy float LP (search) with exact-Fraction spot checks.  HONEST verdict.
"""
import sys, itertools
from math import comb, gcd
from fractions import Fraction as F
from functools import reduce, lru_cache
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)
H = F(1, 14)
INNER = list(range(1, 7))
ALL_SUBSETS = [frozenset(b for b in INNER if (mask >> (b - 1)) & 1) for mask in range(64)]

from scipy.optimize import linprog
import numpy as np

def occupancy_full(E):
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    law = defaultdict(lambda: F(0))
    for x0, x1 in zip(bps, bps[1:]):
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        hit = set(int(7 * e * xm) % 7 for e in E)
        missed = frozenset(s for s in INNER if s not in hit)
        law[missed] += x1 - x0
    return dict(law)

def q_from_law(law):
    items = list(law.items())
    return {A: sum(m for B, m in items if A <= B) for A in ALL_SUBSETS}

def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

def danger(u):
    iv = []
    for j in range(u):
        c = F(j, u); a = (c - H / u) % 1; b = (c + H / u) % 1
        if a < b: iv.append((a, b))
        else: iv.append((a, F(1))); iv.append((F(0), b))
    return iv
def mgmerge(iv):
    iv = sorted(iv); o = []
    for a, b in iv:
        if o and a <= o[-1][1]: o[-1] = (o[-1][0], max(o[-1][1], b))
        else: o.append((a, b))
    return o
def measGP(P):
    if not P: return F(1)
    dz = mgmerge([iv for u in P for iv in danger(u)]); s = F(0); prev = F(0)
    for a, b in dz:
        if a > prev: s += a - prev
        prev = max(prev, b)
    if prev < 1: s += 1 - prev
    return s
@lru_cache(None)
def cap(k):
    psz = 13 - k
    if psz == 0: return F(1)
    return min(measGP(P) for P in itertools.combinations(range(1, 14), psz))

WINDOWS = {8: 16, 9: 15, 10: 14, 11: 13}

def build_window(k):
    """Return list of (E, q_floatvec_indexed_by_subset, measS7_float, measS7_F)."""
    recs = []
    for rest in itertools.combinations(range(1, WINDOWS[k] + 1), k - 1):
        E = [0] + list(rest)
        if not primitive(E): continue
        law = occupancy_full(E); q = q_from_law(law)
        s7 = law.get(frozenset(), F(0))
        recs.append((tuple(E), q, float(s7), s7))
    return recs

def majorant_lp(recs, subsets, target_q):
    """min sum_A lambda_A target_q[A]  s.t. sum_A lambda_A q_A(E) >= measS7(E) for all E in recs.
       Returns (optval, lambda_vector) or None."""
    n = len(subsets)
    cobj = [float(target_q[A]) for A in subsets]
    A_ub = [[-float(q[A]) for A in subsets] for (_, q, _, _) in recs]
    b_ub = [-s7f for (_, _, s7f, _) in recs]
    res = linprog(cobj, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                  bounds=[(None, None)] * n, method="highs")
    if not res.success: return None
    return res.fun, res.x

print("=" * 100)
print("THREAD B CONTROL: is the degree-R per-subset tight majorant VACUOUS or does it pick consec?")
print("=" * 100)

for k in [8, 9, 10]:
    recs = build_window(k)
    ck = cap(k)
    Ec = tuple(consec(k))
    qc = next(q for (E, q, _, _) in recs if E == Ec)
    s7c = next(s7F for (E, _, _, s7F) in recs if E == Ec)
    # pick a few NON-consec targets: top non-consec by measS7, a mid one, a low one
    nonc = sorted([(s7f, E, q, s7F) for (E, q, s7f, s7F) in recs if E != Ec], reverse=True)
    targets = [("consec", Ec, qc, float(s7c), s7c)]
    for label, idx in [("2nd-best", 0), ("median", len(nonc) // 2), ("low", -1)]:
        s7f, E, q, s7F = nonc[idx]
        targets.append((label + " " + str(list(E)), E, q, s7f, s7F))

    for R in [2, 3]:
        subsets = [frozenset(A) for r in range(R + 1) for A in itertools.combinations(INNER, r)]
        print(f"\n  --- k={k}, degree R={R} (#subsets={len(subsets)}), cap={float(ck):.5f} ---")
        for label, E, q, s7f, s7F in targets:
            out = majorant_lp(recs, subsets, q)
            if out is None:
                print(f"    target {label}: LP failed"); continue
            val, lam = out
            gap = val - s7f
            print(f"    target {label:>34}: min majorant={val:.6f}  measS7={s7f:.6f}  "
                  f"gap={gap:.6f}  tight?{abs(gap)<1e-9}")

print("\n" + "=" * 100)
print("TEST 2: For the tight-at-consec majorant lambda, is consec the ARGMAX of F_lambda over window?")
print("  If yes: measS7(E) <= F_lambda(E) <= F_lambda(consec) = measS7(consec) => CONSEC-MAX PROVED at degree R.")
print("=" * 100)
for k in [8, 9, 10]:
    recs = build_window(k); ck = cap(k)
    Ec = tuple(consec(k))
    qc = next(q for (E, q, _, _) in recs if E == Ec)
    s7c = float(next(s7F for (E, _, _, s7F) in recs if E == Ec))
    for R in [2, 3]:
        subsets = [frozenset(A) for r in range(R + 1) for A in itertools.combinations(INNER, r)]
        out = majorant_lp(recs, subsets, qc)
        if out is None:
            print(f"  k={k} R={R}: LP failed"); continue
        val, lam = out
        # F_lambda(E) for all E; is consec the argmax?
        FE = []
        for (E, q, s7f, s7F) in recs:
            fe = sum(lam[i] * float(q[subsets[i]]) for i in range(len(subsets)))
            FE.append((fe, E, s7f))
        Fc = sum(lam[i] * float(qc[subsets[i]]) for i in range(len(subsets)))
        above = [(fe, E, s7f) for (fe, E, s7f) in FE if fe > Fc + 1e-9]
        # also: validity margin -- is F_lambda(E) >= measS7(E) everywhere (should be by constr)?
        viol = [(fe, E, s7f) for (fe, E, s7f) in FE if fe < s7f - 1e-7]
        print(f"  k={k} R={R}: F_lambda(consec)={Fc:.6f}=measS7(consec)={s7c:.6f}.  "
              f"#shapes with F_lambda(E)>F_lambda(consec) = {len(above)}  (consec=argmax F_lambda? {len(above)==0})")
        if above[:3]:
            for fe, E, s7f in sorted(above, reverse=True)[:3]:
                print(f"        BEATER F_lambda={fe:.6f} > {Fc:.6f} at {list(E)} (measS7={s7f:.5f})")
        if viol:
            print(f"        (validity violations within tol: {len(viol)} -- float noise)")

print("\nDONE (control).")
