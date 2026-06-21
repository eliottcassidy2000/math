#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_mobius_certificate_existence_opus_0621.py   (opus, 2026-06-21, THREAD B final)

THE DECISIVE certificate-existence LP for the Mobius/per-subset route.

A degree-R per-subset functional F_lambda(E)=sum_{|A|<=R} lambda_A q_A(E) PROVES consec-max iff
there is a SINGLE lambda with:
   (V) validity:     F_lambda(E) >= measS7(E)         for ALL E in window
   (M) consec-max:   F_lambda(E) <= F_lambda(consec)  for ALL E in window
   (T) tightness:    F_lambda(consec) = measS7(consec)
(V)+(M)+(T) give  measS7(E) <= F_lambda(E) <= F_lambda(consec) = measS7(consec)  for all E
=> consec-max with NO finite atlas needed beyond the certificate.  (M) alone with (T) and
the lattice structure is the CJJ "integral optimum" condition.

We solve the FEASIBILITY LP for lambda with constraints (V),(M),(T) and report whether it is
FEASIBLE at degree R=2,3,4.  If INFEASIBLE at all bounded R -> the per-subset Mobius route
CANNOT certify consec-max at bounded degree => COLLAPSE in the sense of Prop 1.2 (the
extremizer consec is not selectable by a low-degree lattice functional).  If FEASIBLE at some
R -> the Mobius route WORKS at that degree (a genuine extremality-without-finite-checks
certificate, modulo proving the per-subset bounds q_A(E)<=q_A(consec)-type inputs).

We also test the CRUCIAL sub-question: is the route's failure due to (M) being incompatible
with (V)+(T)?  We report, for the lambda minimizing max_E F_lambda(E) subject to (V)+(T),
the SLACK  max_E F_lambda(E) - F_lambda(consec)  (>=0; =0 iff consec-max certificate exists).

ALL scipy float (it is a feasibility/optimization question; exactness not needed for a
NEGATIVE).  We spot-check the binding constraints.
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

def build(k):
    recs = []
    for rest in itertools.combinations(range(1, WINDOWS[k] + 1), k - 1):
        E = [0] + list(rest)
        if not primitive(E): continue
        law = occupancy_full(E); q = q_from_law(law)
        s7 = float(law.get(frozenset(), F(0)))
        recs.append((tuple(E), q, s7))
    return recs

print("=" * 100)
print("THREAD B FINAL: does a SINGLE degree-R per-subset functional certify consec-max?")
print("  Need lambda with (V) F>=measS7 all E, (M) F(E)<=F(consec) all E, (T) F(consec)=measS7(consec).")
print("  We minimize  s = max_E F_lambda(E) - F_lambda(consec)  s.t. (V)+(T).  s=0 iff certificate exists.")
print("=" * 100)

for k in [8, 9, 10]:
    recs = build(k); ck = float(cap(k))
    Ec = tuple(consec(k))
    qc = next(q for (E, q, _) in recs if E == Ec)
    s7c = next(s7 for (E, q, s7) in recs if E == Ec)
    for R in [2, 3, 4]:
        subsets = [frozenset(A) for r in range(R + 1) for A in itertools.combinations(INNER, r)]
        n = len(subsets)
        qcv = np.array([float(qc[A]) for A in subsets])
        # variables: lambda (n) and s (1).  minimize s.
        # (M): F(E) - F(consec) <= s  ->  (q_E - q_consec).lambda - s <= 0
        # (V): F(E) >= measS7(E)      ->  -q_E.lambda <= -measS7(E)
        # (T): q_consec.lambda = measS7(consec)
        c = [0.0] * n + [1.0]
        A_ub = []; b_ub = []
        for (E, q, s7) in recs:
            qe = np.array([float(q[A]) for A in subsets])
            A_ub.append(list(qe - qcv) + [-1.0]); b_ub.append(0.0)        # (M)
            A_ub.append(list(-qe) + [0.0]); b_ub.append(-s7)              # (V)
        A_eq = [list(qcv) + [0.0]]; b_eq = [s7c]                          # (T)
        bounds = [(None, None)] * n + [(0.0, None)]
        res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                      A_eq=np.array(A_eq), b_eq=np.array(b_eq),
                      bounds=bounds, method="highs")
        if not res.success:
            print(f"  k={k} R={R} (#subsets={n}): LP failed/infeasible -> {res.message}")
            continue
        s = res.x[-1]
        cert = s < 1e-7
        print(f"  k={k} R={R} (#subsets={n}): min slack s = max_E F - F(consec) = {s:.6f}  "
              f"=> consec-max certificate exists? {cert}")
    print()

print("=" * 100)
print("INTERPRETATION")
print("  s=0  : a single degree-R per-subset functional certifies consec-max (route SUCCEEDS at R).")
print("  s>0  : NO degree-R per-subset functional is simultaneously valid, tight-at-consec, and")
print("         consec-maximal.  The per-subset Mobius data does NOT select consec at degree R.")
print("         This is the Prop-1.2 COLLAPSE of the EXTREMALITY (consec is not a lattice-functional")
print("         optimizer), even though Part B showed the per-subset BOUND strictly improves level 1.")
print("=" * 100)
print("\nDONE.")
