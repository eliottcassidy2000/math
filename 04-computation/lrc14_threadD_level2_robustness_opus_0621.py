#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadD_level2_robustness_opus_0621.py  (opus, 2026-06-21, THREAD D)

ROBUSTNESS of the headline Thread-D finding:
  the LEVEL-2 LP bound (pair co-occupancy = support-2 relation-code) is maximized by consec,
  while the LEVEL-1 bound is NOT.

We harden this with:
  (1) WIDER windows + larger k (8,9,10,11,12,13) for the level-2 arg-max test.
  (2) The k=11 BOUNDARY CASE: the THM-534 factorial-moment dual L_y had 1 beater at k=11
      (consec not L_y-max).  Does the FULL level-2 individual LP also have a beater?  This
      separates "symmetric-moment level-2" (S_1,S_2 only) from "full individual level-2"
      (all q_{st}).  If full level-2 has NO beater where symmetric-moment did, the EXTRA
      relation information (which pairs co-miss) is what saves consec.
  (3) The SYMMETRIC-moment level-2 LP_sym(2) (fix S_1,S_2 only) vs FULL individual LP_ind(2):
      is the extra "which-pair" info necessary, or does S_1,S_2 already suffice?

All exact occupancy (Fractions); LP via scipy highs.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
import numpy as np
from scipy.optimize import linprog
sys.stdout.reconfigure(line_buffering=True)

SUBSETS = [frozenset(c) for r in range(7) for c in itertools.combinations(range(1, 7), r)]
SUB_IDX = {A: i for i, A in enumerate(SUBSETS)}
EMPTY = SUB_IDX[frozenset()]
SUBSET_OF = {A: [SUB_IDX[B] for B in SUBSETS if A <= B] for A in SUBSETS}

def joint_occupancy(E):
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for a in range(7 * e + 1): bps.add(F(a, 7 * e))
    bps = sorted(bps)
    p = [F(0)] * 7
    rho = {A: F(0) for A in SUBSETS}
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2; L = x1 - x0
        hit = set(int(7 * e * xm) % 7 for e in E)
        missed = sorted(s for s in range(1, 7) if s not in hit)
        p[len(missed)] += L
        for r in range(len(missed) + 1):
            for c in itertools.combinations(missed, r):
                rho[frozenset(c)] += L
    return p, rho

def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1

def LP_individual(rho, lvl):
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

def LP_symmetric(rho, R):
    S = {r: sum(rho[A] for A in SUBSETS if len(A) == r) for r in range(7)}
    nvar = len(SUBSETS)
    c = np.zeros(nvar); c[EMPTY] = -1.0
    A_eq = [[1.0] * nvar]; b_eq = [1.0]
    for r in range(1, R + 1):
        row = np.array([float(comb(len(B), r)) for B in SUBSETS])
        A_eq.append(row); b_eq.append(float(S[r]))
    res = linprog(c, A_eq=np.array(A_eq), b_eq=np.array(b_eq),
                  bounds=[(0, 1)] * nvar, method="highs")
    return -res.fun if res.success else None

def banner(t): print("\n" + "=" * 92 + f"\n{t}\n" + "=" * 92)

if __name__ == "__main__":
    TOL = 1e-7
    banner("LEVEL-2 ARG-MAX robustness: consec vs all primitive shapes (wider windows)")
    print("  Bnd1 = LP_ind(1) (singleton-only).  Bnd2 = LP_ind(2) (pairs).  Sym2 = LP_sym(2) (S1,S2).")
    print("  Report: #shapes beating consec on each bound, and on measS7 truth.\n")
    WINDOWS = {8: 16, 9: 14, 10: 13, 11: 13, 12: 13, 13: 13}
    for k in sorted(WINDOWS):
        maxE = WINDOWS[k]; C = list(range(k))
        pC, rhoC = joint_occupancy(C); measC = float(pC[0])
        b1C, b2C, s2C = LP_individual(rhoC, 1), LP_individual(rhoC, 2), LP_symmetric(rhoC, 2)
        beat1 = beat2 = beats2 = beatM = nset = 0
        worst2 = (b2C, C)
        for rest in itertools.combinations(range(1, maxE + 1), k - 1):
            E = [0] + list(rest)
            if not primitive(E): continue
            nset += 1
            p, rho = joint_occupancy(E); meas = float(p[0])
            if meas > measC + TOL: beatM += 1
            if LP_individual(rho, 1) > b1C + TOL: beat1 += 1
            v2 = LP_individual(rho, 2)
            if v2 > b2C + TOL: beat2 += 1
            if v2 > worst2[0]: worst2 = (v2, E)
            if LP_symmetric(rho, 2) > s2C + TOL: beats2 += 1
        print(f"k={k} (maxE<={maxE}, {nset} shapes): measS7(consec)={measC:.6f}")
        print(f"   truth   beaters={beatM:4d}  consec-max-meas:{beatM==0}")
        print(f"   Bnd1    consec={b1C:.5f} beaters={beat1:4d}  consec-max:{beat1==0}")
        print(f"   Sym2    consec={s2C:.5f} beaters={beats2:4d}  consec-max:{beats2==0}")
        print(f"   Bnd2    consec={b2C:.5f} beaters={beat2:4d}  consec-max:{beat2==0}  argmax={worst2[1]}")
        print()
    print("DONE.")
