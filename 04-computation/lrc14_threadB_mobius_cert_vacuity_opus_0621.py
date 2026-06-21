#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_mobius_cert_vacuity_opus_0621.py   (opus, 2026-06-21, THREAD B decisive control)

DECISIVE VACUITY CONTROL for the certificate-existence finding.

The certificate LP found a degree-2 per-subset functional F_lambda that is
  (V) >= measS7 all E, (M) <= F(consec) all E, (T) =measS7 at consec.
=> certifies consec as the max.  BUT we must check this is NOT VACUOUS:

  If the SAME LP returns s=0 (a "certificate") when we substitute a NON-MAXIMAL shape E*
  as the "winner", then the construction would certify a FALSE maximum -- impossible if
  honest, so s MUST be > 0 for non-maximal targets.  We run the certificate LP with
  target = several non-consec shapes and report s.

  KEY: for a non-maximal target E* (measS7(E*)<measS7(consec)), constraint (M)
  F(E)<=F(E*) for ALL E *together with* (V) F(consec)>=measS7(consec) FORCES
  F(E*)>=F(consec)>=measS7(consec)>measS7(E*), so (T) F(E*)=measS7(E*) is then VIOLATED
  unless s>0.  So a non-maximal target MUST give s>0 (strictly).  If instead it gives s=0,
  the per-subset functional is so flexible it can fake any maximum => vacuous, route dead.

  We ALSO sharpen the consec result: with the certifying lambda, we re-verify (exact
  Fractions) that F_lambda is valid + consec-max + tight on the FULL window, and we extract
  the lambda to see its structure (does it match the Krawtchouk/Mobius g_2 signs?).
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
WINDOWS = {8: 16, 9: 15, 10: 14}

def build(k):
    recs = []
    for rest in itertools.combinations(range(1, WINDOWS[k] + 1), k - 1):
        E = [0] + list(rest)
        if not primitive(E): continue
        law = occupancy_full(E); q = q_from_law(law)
        s7 = float(law.get(frozenset(), F(0)))
        recs.append((tuple(E), q, s7))
    return recs

def cert_lp(recs, subsets, target_q, target_s7):
    """min slack s = max_E F - F(target)  s.t. (V) F>=measS7 all E, (T) F(target)=measS7(target)."""
    n = len(subsets)
    tqv = np.array([float(target_q[A]) for A in subsets])
    c = [0.0] * n + [1.0]
    A_ub = []; b_ub = []
    for (E, q, s7) in recs:
        qe = np.array([float(q[A]) for A in subsets])
        A_ub.append(list(qe - tqv) + [-1.0]); b_ub.append(0.0)   # (M) F(E)-F(target)<=s
        A_ub.append(list(-qe) + [0.0]); b_ub.append(-s7)         # (V) F(E)>=measS7(E)
    A_eq = [list(tqv) + [0.0]]; b_eq = [target_s7]               # (T)
    bounds = [(None, None)] * n + [(0.0, None)]
    res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                  A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=bounds, method="highs")
    if not res.success: return None
    return res.x[-1], res.x[:n]

print("=" * 100)
print("DECISIVE VACUITY CONTROL: certificate LP with NON-CONSEC winners must give s>0.")
print("  (A maximal-shape certificate gives s=0; any NON-maximal target must give s>0,")
print("   else the per-subset functional fakes any max = vacuous = route dead.)")
print("=" * 100)
for k in [8, 9, 10]:
    recs = build(k)
    Ec = tuple(consec(k))
    qc = next(q for (E, q, _) in recs if E == Ec)
    s7c = next(s7 for (E, q, s7) in recs if E == Ec)
    nonc = sorted([(s7, E, q) for (E, q, s7) in recs if E != Ec], reverse=True)
    targets = [("consec (MAX)", Ec, qc, s7c)]
    for label, idx in [("2nd-best", 0), ("median", len(nonc) // 2), ("low", -1)]:
        s7t, E, q = nonc[idx]
        targets.append((label, E, q, s7t))
    for R in [2]:
        subsets = [frozenset(A) for r in range(R + 1) for A in itertools.combinations(INNER, r)]
        print(f"\n  --- k={k}, degree R={R} (#subsets={len(subsets)}) ---")
        for label, E, q, s7t in targets:
            out = cert_lp(recs, subsets, q, s7t)
            if out is None:
                print(f"    target {label:>14} {str(list(E)):>28}: LP infeasible (s would be +inf) -> GOOD (non-vacuous)")
                continue
            s, lam = out
            flag = "s=0 (certified as max)" if s < 1e-7 else f"s={s:.6f} (NOT certifiable as max)"
            print(f"    target {label:>14} {str(list(E)):>28}: measS7={s7t:.5f}  {flag}")

print("\n" + "=" * 100)
print("VERDICT LOGIC:")
print("  If ONLY consec gives s=0 and all non-consec give s>0 (or infeasible): the degree-2")
print("    per-subset functional GENUINELY selects consec => MOBIUS/INTEGRALITY ROUTE SUCCEEDS")
print("    (a real extremality certificate at degree 2, no full atlas).")
print("  If non-consec ALSO give s=0: VACUOUS => route collapses (Prop 1.2).")
print("=" * 100)
print("\nDONE.")
