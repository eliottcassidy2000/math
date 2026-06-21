#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadD_hierarchy_argmax_opus_0621.py  (opus, 2026-06-21, THREAD D)

THE NESTED LP HIERARCHY for measS7, and the ARG-MAX test: is consec the maximizer of the
level-l LP bound?  (CJJ arXiv:2211.01248 -- collapse-vs-improves + does a low level pin consec.)

Two parallel hierarchies on the joint missed-sector law (variables = P(M=A), A subset {1..6}):

  SYMMETRIC level-R  (= THM-534 / Delsarte family): fix the factorial moments S_1..S_R only.
     S_r = sum_{|A|=r} P(... A subset M) = E[C(N,r)].  Shape enters via S_1..S_R (R numbers).
     LP_sym(R; E) = max P(M=emptyset) s.t. sum=1, S_r(law)=S_r(E) r=1..R, law>=0.

  INDIVIDUAL (relation-code) level-l: fix ALL l-subset marginals P(A subset M)=:rho_A, |A|<=l.
     rho_{{s}}=m_s; rho_{{s,t}}=q_{st}; etc.  Shape enters via the INDIVIDUAL marginals, i.e.
     WHICH sectors co-miss (the Lambda(E) support-<=l relation structure), not just their sums.
     LP_ind(l; E) = max P(M=emptyset) s.t. sum=1, all |A|<=l marginals matched, law>=0.

Nesting:  LP_ind(l) <= LP_sym(l) (individual fixes more than the sum)  and both decrease in
their level.  LP_ind(6) = LP_sym(6) = measS7 (all marginals fixed => law determined => tight).

DELIVERABLE (Thread D):
  (A) COLLAPSE-vs-IMPROVES: does each level strictly improve?  (Per CJJ Prop 1.2, improvement
      requires the optimizer to be closed under LINEAR combinations -- Lambda(E) IS linear, so
      we expect improvement; a collapse would be the Prop-1.2 negative.)
  (B) ARG-MAX: as E ranges over primitive shapes, is the level-l LP BOUND maximized by consec?
      (If the level-l bound itself is consec-extremal, then a per-level inequality + monotone
      argument would pin consec.)  We sweep bounded spread and count beaters at each level.
  (C) The LEVEL where LP_ind(l;consec) == measS7(consec) (consec saturates).  If finite and
      small, that level "pins consec" with a finite-information certificate.

All occupancy EXACT (Fractions); LP via scipy highs (floats; gaps reported with tolerance).
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
from collections import defaultdict
import numpy as np
from scipy.optimize import linprog
sys.stdout.reconfigure(line_buffering=True)

SUBSETS = [frozenset(c) for r in range(7) for c in itertools.combinations(range(1, 7), r)]
SUB_IDX = {A: i for i, A in enumerate(SUBSETS)}
EMPTY = SUB_IDX[frozenset()]

def joint_occupancy(E):
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for a in range(7 * e + 1): bps.add(F(a, 7 * e))
    bps = sorted(bps)
    p = [F(0)] * 7
    rho = {A: F(0) for A in SUBSETS}     # rho[A] = P(A subset M)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2; L = x1 - x0
        hit = set(int(7 * e * xm) % 7 for e in E)
        missed = frozenset(s for s in range(1, 7) if s not in hit)
        p[len(missed)] += L
        # add L to rho[A] for every A subset missed
        ml = sorted(missed)
        for r in range(len(ml) + 1):
            for c in itertools.combinations(ml, r):
                rho[frozenset(c)] += L
    return p, rho

def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

# ---- precompute subset-membership masks for marginal rows ----
SUBSET_OF = {}   # SUBSET_OF[A] = list of indices i with A subset SUBSETS[i]
for A in SUBSETS:
    SUBSET_OF[A] = [SUB_IDX[B] for B in SUBSETS if A <= B]

def LP_individual(E, lvl, rho=None):
    """max P(M=emptyset) s.t. sum=1, P(A subset M)=rho[A] for all |A|<=lvl, law>=0."""
    if rho is None: _, rho = joint_occupancy(E)
    nvar = len(SUBSETS)
    c = np.zeros(nvar); c[EMPTY] = -1.0
    A_eq = [[1.0] * nvar]; b_eq = [1.0]
    for A in SUBSETS:
        if 1 <= len(A) <= lvl:
            row = np.zeros(nvar)
            for i in SUBSET_OF[A]: row[i] = 1.0
            A_eq.append(row); b_eq.append(float(rho[A]))
    res = linprog(c, A_eq=np.array(A_eq), b_eq=np.array(b_eq),
                  bounds=[(0, 1)] * nvar, method="highs")
    return -res.fun if res.success else None

def LP_symmetric(E, R, rho=None):
    """max P(M=emptyset) s.t. sum=1, S_r(law)=S_r(E) for r=1..R, law>=0."""
    if rho is None: _, rho = joint_occupancy(E)
    # S_r(E) = sum_{|A|=r} rho[A]
    S = {r: sum(rho[A] for A in SUBSETS if len(A) == r) for r in range(7)}
    nvar = len(SUBSETS)
    c = np.zeros(nvar); c[EMPTY] = -1.0
    A_eq = [[1.0] * nvar]; b_eq = [1.0]
    for r in range(1, R + 1):
        row = np.zeros(nvar)
        for i, B in enumerate(SUBSETS):
            row[i] = float(comb(len(B), r))   # #r-subsets of B = C(|B|,r)
        A_eq.append(row); b_eq.append(float(S[r]))
    res = linprog(c, A_eq=np.array(A_eq), b_eq=np.array(b_eq),
                  bounds=[(0, 1)] * nvar, method="highs")
    return -res.fun if res.success else None

def banner(t): print("\n" + "=" * 92 + f"\n{t}\n" + "=" * 92)

# ============================================================================
if __name__ == "__main__":
    TOL = 1e-7
    banner("(A) COLLAPSE-vs-IMPROVES: nested LP bounds for CONSEC at each level")
    print("  LP_sym(R)  fixes factorial moments S_1..S_R  (the Delsarte/THM-534 family).")
    print("  LP_ind(l)  fixes ALL individual marginals P(A subset M), |A|<=l (relation-code).")
    print("  Both -> measS7 at full level (R=6 / l=6).  IMPROVES iff each level strictly drops.\n")
    for k in [8, 9, 10]:
        C = consec(k); p, rho = joint_occupancy(C); meas = float(p[0])
        print(f"--- k={k} consec, measS7={meas:.6f} ---")
        symvals = [LP_symmetric(C, R, rho) for R in range(1, 7)]
        indvals = [LP_individual(C, l, rho) for l in range(1, 7)]
        print("  LP_sym(R): " + "  ".join(f"R{R}={v:.5f}" for R, v in zip(range(1,7), symvals)))
        print("  LP_ind(l): " + "  ".join(f"l{l}={v:.5f}" for l, v in zip(range(1,7), indvals)))
        # strict improvement flags
        sym_imp = [symvals[i] < symvals[i-1] - TOL for i in range(1, 6)]
        ind_imp = [indvals[i] < indvals[i-1] - TOL for i in range(1, 6)]
        print(f"  sym strictly improves R->R+1: {sym_imp}")
        print(f"  ind strictly improves l->l+1: {ind_imp}")
        # level where ind PINS consec
        pin = next((l for l in range(1, 7) if abs(indvals[l-1] - meas) < TOL), None)
        pin_sym = next((R for R in range(1, 7) if abs(symvals[R-1] - meas) < TOL), None)
        print(f"  level l where LP_ind PINS consec (==measS7): {pin}")
        print(f"  level R where LP_sym PINS consec (==measS7): {pin_sym}\n")

    banner("(B) ARG-MAX TEST: is consec the maximizer of the level-l LP BOUND over shapes?")
    print("  For each shape E the level-l LP bound is a self-upper-bound on measS7(E).")
    print("  If consec maximizes the BOUND, a per-level inequality could pin consec.")
    print("  We sweep bounded spread; report beaters (shapes whose level-l bound exceeds consec's).\n")
    WINDOWS = {8: 14, 9: 12, 10: 12}
    for k in [8, 9]:
        maxE = WINDOWS[k]; C = consec(k)
        pC, rhoC = joint_occupancy(C); measC = float(pC[0])
        # consec bounds at each level
        bC_sym = {R: LP_symmetric(C, R, rhoC) for R in [1, 2, 3]}
        bC_ind = {l: LP_individual(C, l, rhoC) for l in [1, 2, 3]}
        beat_sym = {R: 0 for R in [1, 2, 3]}
        beat_ind = {l: 0 for l in [1, 2, 3]}
        meas_beat = 0
        argmax_ind2 = (measC, C)
        nset = 0
        for rest in itertools.combinations(range(1, maxE + 1), k - 1):
            E = [0] + list(rest)
            if not primitive(E): continue
            nset += 1
            p, rho = joint_occupancy(E); meas = float(p[0])
            if meas > measC + TOL: meas_beat += 1
            for R in [1, 2, 3]:
                if LP_symmetric(E, R, rho) > bC_sym[R] + TOL: beat_sym[R] += 1
            for l in [1, 2, 3]:
                v = LP_individual(E, l, rho)
                if l == 2 and v > argmax_ind2[0]: argmax_ind2 = (v, E)
                if v > bC_ind[l] + TOL: beat_ind[l] += 1
        print(f"--- k={k} (maxE<={maxE}, {nset} primitive shapes) ---")
        print(f"  measS7(consec)={measC:.6f}; #shapes with HIGHER measS7 = {meas_beat} "
              f"(consec-max truth: {meas_beat==0})")
        for R in [1, 2, 3]:
            print(f"  LP_sym(R={R}) consec={bC_sym[R]:.5f}  beaters={beat_sym[R]}  "
                  f"consec-max-of-bound: {beat_sym[R]==0}")
        for l in [1, 2, 3]:
            print(f"  LP_ind(l={l}) consec={bC_ind[l]:.5f}  beaters={beat_ind[l]}  "
                  f"consec-max-of-bound: {beat_ind[l]==0}")
        print(f"  level-2 individual bound argmax: {argmax_ind2[1]} -> {argmax_ind2[0]:.5f}\n")
    print("DONE.")
