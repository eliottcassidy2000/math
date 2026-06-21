#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_mobius_symmetric_vs_persubset_opus_0621.py   (opus, 2026-06-21, THREAD B crux)

THE CRUX of the deliverable.  The exact-verify run showed the degree-2 consec-max certifying
functional, when allowed per-subset weights, came out SYMMETRIC-PER-SIZE -- i.e. it is just a
MOMENT functional  F = sum_r y_r S_r  (= level-1 Delsarte), NOT genuinely per-subset.

So the SHARP question for "collapse vs improves":
   Does a SYMMETRIC (per-size = moment = LEVEL-1 Delsarte) degree-R functional ALREADY certify
   consec-max on the FULL window?
     YES  => the per-subset/Mobius (level-l) structure adds NOTHING for the EXTREMALITY
             certificate.  COLLAPSE to level 1 (Prop 1.2), even though Part B showed the
             per-subset BOUND is strictly tighter in self-consistency.  Consec-max IS a
             level-1 Delsarte-LP fact -- which is exactly what HYP-2726/THM-534 already say,
             reproved here as a single moment functional that is valid+tight+consec-maximal.
     NO   => symmetric fails but per-subset succeeds => the Mobius route GENUINELY improves
             the EXTREMALITY certificate (the real CJJ win).

We solve, on the FULL window, the certificate slack-LP
   min s = max_E F(E) - F(consec)   s.t. (V) F>=measS7 all E, (T) F(consec)=measS7(consec)
in TWO variable regimes:
   (S) SYMMETRIC: F = sum_{r=0}^R y_r S_r(E)   (R+1 vars y_r)   [= moment / level-1]
   (P) PER-SUBSET: F = sum_{|A|<=R} lambda_A q_A(E)              [the Mobius route]
and report s_sym vs s_persubset for R=2,3,4.  s=0 <=> certificate exists in that regime.

ALL scipy float (a slack-comparison; exactness not needed for the verdict).  We also report,
for the SYMMETRIC certificate, the y_r (and whether they match the even-Bonferroni / Krawtchouk
dual already on file).
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
def Svec_from_q(q):
    return [sum(q[frozenset(A)] for A in itertools.combinations(INNER, r)) for r in range(7)]
def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))
WINDOWS = {8: 16, 9: 15, 10: 14}

def build(k):
    recs = []
    for rest in itertools.combinations(range(1, WINDOWS[k] + 1), k - 1):
        E = [0] + list(rest)
        if not primitive(E): continue
        law = occupancy_full(E); q = q_from_law(law)
        recs.append((tuple(E), q, float(law.get(frozenset(), F(0)))))
    return recs

def slack_lp(recs, feats_consec, feats_of, s7c):
    """min s s.t. for all E: feat(E).x - feat(consec).x <= s  AND feat(E).x >= measS7(E);
       feat(consec).x = measS7(consec)=s7c.  Returns (s, x) or (None,None)."""
    n = len(feats_consec)
    fc = np.array(feats_consec, dtype=float)
    c = [0.0] * n + [1.0]
    A_ub = []; b_ub = []
    for (E, q, s7) in recs:
        fe = np.array(feats_of(q), dtype=float)
        A_ub.append(list(fe - fc) + [-1.0]); b_ub.append(0.0)
        A_ub.append(list(-fe) + [0.0]); b_ub.append(-s7)
    A_eq = [list(fc) + [0.0]]; b_eq = [s7c]
    res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), A_eq=np.array(A_eq), b_eq=np.array(b_eq),
                  bounds=[(None, None)] * n + [(0.0, None)], method="highs")
    if not res.success: return None, None
    return res.x[-1], res.x[:n]

print("=" * 100)
print("THREAD B CRUX: SYMMETRIC (moment / level-1 Delsarte) vs PER-SUBSET (Mobius) consec-max certificate.")
print("  s_sym = slack of best degree-R MOMENT functional;  s_sub = slack of best per-subset functional.")
print("  s=0 => certificate exists.  s_sym=0 => level-1 Delsarte already certifies (COLLAPSE).")
print("=" * 100)
for k in [8, 9, 10]:
    consec_k = consec(k)
    recs = build(k)
    qc = next(q for (E, q, _) in recs if E == tuple(consec_k))
    s7c = next(s7 for (E, q, s7) in recs if E == tuple(consec_k))
    Sc = Svec_from_q(qc)
    print(f"\n  k={k}  (#shapes={len(recs)})")
    for R in [2, 3, 4]:
        subsets = [frozenset(A) for r in range(R + 1) for A in itertools.combinations(INNER, r)]
        # SYMMETRIC features: S_0..S_R
        feats_consec_sym = [float(Sc[r]) for r in range(R + 1)]
        def feats_sym(q, R=R):
            return [float(sum(q[frozenset(A)] for A in itertools.combinations(INNER, r))) for r in range(R + 1)]
        s_sym, ysym = slack_lp(recs, feats_consec_sym, feats_sym, s7c)
        # PER-SUBSET features: q_A, |A|<=R
        feats_consec_sub = [float(qc[A]) for A in subsets]
        def feats_sub(q, subsets=subsets):
            return [float(q[A]) for A in subsets]
        s_sub, lam = slack_lp(recs, feats_consec_sub, feats_sub, s7c)
        ss = f"{s_sym:.6f}" if s_sym is not None else "infeas"
        sp = f"{s_sub:.6f}" if s_sub is not None else "infeas"
        verdict = ""
        if s_sym is not None and s_sub is not None:
            if s_sym < 1e-7 and s_sub < 1e-7: verdict = "BOTH certify -> COLLAPSE to level-1 (moment suffices)"
            elif s_sym >= 1e-7 and s_sub < 1e-7: verdict = "per-subset certifies, moment does NOT -> Mobius IMPROVES"
            elif s_sym >= 1e-7 and s_sub >= 1e-7: verdict = "NEITHER certifies at this R"
        print(f"    R={R}: s_sym={ss:>10}  s_persubset={sp:>10}   -> {verdict}")
        if R == 2 and s_sym is not None and s_sym < 1e-7:
            yrat = [F(v).limit_denominator(10000) for v in ysym]
            print(f"        symmetric y_r (degree-2 moment cert): {[str(v) for v in yrat]}  "
                  f"(even-Bonferroni would be [1,-1,1])")

print("\n" + "=" * 100)
print("VERDICT: if s_sym=0 at low R for all k, consec-max is a LEVEL-1 (moment/Delsarte) certificate;")
print("  the per-subset Mobius lattice adds NOTHING to the extremality => Prop-1.2 COLLAPSE confirmed.")
print("  (Note: a degree-2 MOMENT certificate that is valid+tight+consec-maximal = a clean, finite,")
print("   single-inequality proof of consec-max on the bounded window -- modulo proving the moment")
print("   inequalities S_r(E) feed it, which is the HYP-2726 Delsarte-dual already on file.)")
print("=" * 100)
print("\nDONE.")
