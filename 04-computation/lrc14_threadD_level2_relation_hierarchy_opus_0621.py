#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadD_level2_relation_hierarchy_opus_0621.py  (opus, 2026-06-21, THREAD D)

THE RELATION-CODE LP HIERARCHY for measS7 (Coregliano-Jeronimo-Jones, arXiv:2211.01248).

QUESTION (gap #4, the genuinely-aggregate one): does ADDING the level-2 relation-code
constraints (support-2 / support-3 relation-combination = pair/triple co-occupancy /
MacWilliams positivity on Lambda(E)) to the level-1 Delsarte LP tighten the bound toward
measS7(consec), and at what level does it PIN consec?  Per CJJ Prop 1.2 the hierarchy
improves over level-1 ONLY through LINEARITY of the optimizer.  Lambda(E) IS linear, so
improvement is possible in principle -- we test whether it actually happens HERE.

SETUP.  For a shape E={e_1<...<e_k} (e_1=0), x~Unif[0,1).  Sector s in Z/7 is "missed" if
no e_i x lands in [s/7,(s+1)/7) mod 1.  Sector 0 always hit (e=0).  Inner sectors S={1..6}.
M(x) = set of MISSED inner sectors; N(x)=|M(x)|; measS7(E)=P(N=0)=P(M=emptyset).

THE LP VIEW (the thing the hierarchy operates on).  The DATA at each level:
  Level-0: P(N=0) bounded only by 0<=.<=1.
  Level-1 (DELSARTE): marginal occupancy p_t=P(N=t), t=0..6.  Bound = the moment/Krawtchouk
    LP: max p_0 s.t. moment constraints S_r=sum_t C(t,r)p_t and Krawtchouk positivity, with
    the SHAPE entering only through S_1,...,S_R (R=#binding factorial moments).  This is what
    THM-534 / HYP-2726 / Thread B already do: measS7<=L_y=sum_j c_j M_j.
  Level-2 (RELATION-CODE): the JOINT occupancy -- singleton marginals m_s=P(s in M) AND pair
    co-occupancy q_{st}=P({s,t} subset M).  These are *exactly* the support-1 and support-2
    relation statistics of Lambda(E) (HYP-2719/2723): m_s and q_{st} are Weyl integrals whose
    relation expansion is graded by support.  The level-2 LP maximizes P(M=emptyset) given the
    REALIZABLE (m_s,q_{st}) plus the Bonferroni / covariance / PSD constraints they must satisfy.

HOW LAMBDA(E) ENTERS.  q_{st}=P(sectors s,t BOTH missed) depends on the additive relations
among {e_i mod 7}.  consec={0..k-1} has residues mod 7 = {0,1,2,3,4,5,6}=ALL of Z/7 (full,
anti-MDS): every sector is hit by *some* generator at the lattice level, and the relation code
has MINIMUM distance 2 (the support-2 relation 2*e_1=e_2 etc.).  We compute the EXACT
(m_s,q_{st}) for consec and for competitors, plug them into the level-2 LP, and compare the LP
optimum to measS7(consec).

DELIVERABLE: collapse-vs-improves at level 2, and whether the realizable level-2 data PINS
consec (level needed).  HONEST about negatives (Prop 1.2 collapse is a valuable result).
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
from collections import defaultdict
import numpy as np
from scipy.optimize import linprog
sys.stdout.reconfigure(line_buffering=True)

# ============================================================================
# EXACT occupancy machinery: joint missed-sector statistics for a shape E.
# ============================================================================
def joint_occupancy(E):
    """Return exact statistics of the missed-inner-sector set M(x) over x in [0,1).
       p[t]=P(N=t);  m[s]=P(s in M), s=1..6;  q[(s,t)]=P({s,t} subset M), 1<=s<t<=6.
       Also pset: dict frozenset(M)->measure (full joint law on subsets of {1..6})."""
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for a in range(7 * e + 1):
            bps.add(F(a, 7 * e))
    bps = sorted(bps)
    p = [F(0)] * 7
    m = {s: F(0) for s in range(1, 7)}
    q = {(s, t): F(0) for s in range(1, 7) for t in range(s + 1, 7)}
    pset = defaultdict(lambda: F(0))
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2; L = x1 - x0
        hit = set(int(7 * e * xm) % 7 for e in E)   # includes 0 from e=0
        missed = frozenset(s for s in range(1, 7) if s not in hit)
        p[len(missed)] += L
        for s in missed: m[s] += L
        for s, t in itertools.combinations(sorted(missed), 2): q[(s, t)] += L
        pset[missed] += L
    return p, m, q, dict(pset)

def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

# ============================================================================
# Krawtchouk (binary, 6 inner sectors) for the level-1 dual, exact.
# ============================================================================
def Kraw(j, t, n=6):
    return sum((-1) ** i * comb(t, i) * comb(n - t, j - i) for i in range(j + 1))
KTAB = [[F(Kraw(j, t)) for t in range(7)] for j in range(7)]

DUAL_G = {  # THM-534 level-1 moment-LP duals g(t), t=0..6 (Krawtchouk-nonneg)
    8:  [F((t-1)*(t-2)*(t-4)*(t-5), 40) for t in range(7)],
    9:  [F(-(t-2)*(t-3)*(t-6), 36) for t in range(7)],
    10: [F(-(t-2)*(t-3)*(t-6), 36) for t in range(7)],
    11: [F((t-3)*(t-4), 12) for t in range(7)],
    12: [F((t-3)*(t-4), 12) for t in range(7)],
    13: [F((t-3)*(t-4), 12) for t in range(7)],
}
def L_y_level1(p, k):
    g = DUAL_G[k]; return sum(g[t] * p[t] for t in range(7))

# ============================================================================
# THE LEVEL-1 LP (Delsarte): max p_0 over abstract occupancy p_t given the
# binding factorial moments S_1..S_R of E.  Shape enters only via those moments.
# This is the LP whose dual is L_y.  We solve the PRIMAL to get the LP optimum.
# Constraints: p_t>=0, sum p_t=1, S_r(p)=S_r(E) for r=1..R (the binding moments).
# ============================================================================
def binding_R(k):
    """Highest factorial moment used in the level-1 dual = degree of g."""
    g = DUAL_G[k]
    # express g in factorial basis C(t,r); degree = max r with nonzero coeff
    # g is a polynomial in t of known degree: k=8 deg4, k=9,10 deg3, k=11-13 deg2
    return {8: 4, 9: 3, 10: 3, 11: 2, 12: 2, 13: 2}[k]

def Svec(p, R):
    return [sum(F(comb(t, r)) * p[t] for t in range(7)) for r in range(R + 1)]

def level1_LP_opt(E, k):
    """max p_0 s.t. p>=0, sum p=1, S_r(p)=S_r(E) for r=1..R.  scipy float; report bound."""
    p_E, _, _, _ = joint_occupancy(E)
    R = binding_R(k)
    S_E = [float(sum(comb(t, r) * p_E[t] for t in range(7))) for r in range(R + 1)]
    # variables p_0..p_6, maximize p_0 => minimize -p_0
    c = np.zeros(7); c[0] = -1.0
    A_eq = [[1.0] * 7]; b_eq = [1.0]
    for r in range(1, R + 1):
        A_eq.append([float(comb(t, r)) for t in range(7)]); b_eq.append(S_E[r])
    res = linprog(c, A_eq=np.array(A_eq), b_eq=np.array(b_eq),
                  bounds=[(0, 1)] * 7, method="highs")
    return -res.fun if res.success else None

# ============================================================================
# THE LEVEL-2 LP (relation-code): max P(M=emptyset) given the REALIZABLE
# singleton marginals m_s=P(s in M) and pair co-occupancy q_{st}=P({s,t} subset M).
# This is the genuine "level-2" of CJJ: the data is all 2-subset statistics
# (= support-<=2 relation-combination statistics on Lambda(E)).
#
# Variables: the full joint law on subsets of {1..6}: P(M=A) for A subset {1..6}.
#   (64 variables.)  Objective: max P(M=emptyset).
# Constraints (the LEVEL-2 INFORMATION = what is FIXED by the shape):
#   (E1) sum_A P(M=A) = 1
#   (E2) for each s:        sum_{A: s in A} P(M=A) = m_s        (singleton marginal)
#   (E3) for each s<t:      sum_{A: {s,t} subset A} P(M=A) = q_{st}   (pair marginal)
#   P(M=A) >= 0.
# The LP optimum = the BEST P(M=emptyset) consistent with the fixed 1- and 2-marginals.
# This is exactly: how much does knowing the pair co-occupancy structure (Lambda(E)
# support-2) pin down the covering probability?  We compare LP_opt to actual measS7(E).
#
# *** CRUCIAL: the level-2 LP is shape-dependent ONLY through (m_s, q_{st}). ***
# If LP_opt(consec) < LP_opt(competitor) for some competitor with HIGHER measS7 envelope,
# the hierarchy "sees" the shape only through these and we test whether consec is the
# arg-max of the LEVEL-2 LP BOUND (an upper bound on its own measS7).
# ============================================================================
SUBSETS = [frozenset(c) for r in range(7) for c in itertools.combinations(range(1, 7), r)]
SUB_IDX = {A: i for i, A in enumerate(SUBSETS)}  # 64 subsets

def level2_LP_opt(E):
    """max P(M=emptyset) over joint laws matching (m_s,q_{st}) of E.  scipy float."""
    _, m, q, _ = joint_occupancy(E)
    nvar = len(SUBSETS)
    c = np.zeros(nvar); c[SUB_IDX[frozenset()]] = -1.0   # maximize P(emptyset)
    A_eq = []; b_eq = []
    A_eq.append([1.0] * nvar); b_eq.append(1.0)
    for s in range(1, 7):
        row = [1.0 if s in A else 0.0 for A in SUBSETS]
        A_eq.append(row); b_eq.append(float(m[s]))
    for s in range(1, 7):
        for t in range(s + 1, 7):
            row = [1.0 if (s in A and t in A) else 0.0 for A in SUBSETS]
            A_eq.append(row); b_eq.append(float(q[(s, t)]))
    res = linprog(c, A_eq=np.array(A_eq), b_eq=np.array(b_eq),
                  bounds=[(0, 1)] * nvar, method="highs")
    return -res.fun if res.success else None

def level1_marginal_LP_opt(E):
    """LEVEL-1 ANALOG with the SAME variable space (joint law on subsets) but ONLY the
       SINGLETON marginals m_s fixed (no pair info).  Isolates the marginal-only bound."""
    _, m, q, _ = joint_occupancy(E)
    nvar = len(SUBSETS)
    c = np.zeros(nvar); c[SUB_IDX[frozenset()]] = -1.0
    A_eq = [[1.0] * nvar]; b_eq = [1.0]
    for s in range(1, 7):
        row = [1.0 if s in A else 0.0 for A in SUBSETS]
        A_eq.append(row); b_eq.append(float(m[s]))
    res = linprog(c, A_eq=np.array(A_eq), b_eq=np.array(b_eq),
                  bounds=[(0, 1)] * nvar, method="highs")
    return -res.fun if res.success else None

def banner(t): print("\n" + "=" * 92 + f"\n{t}\n" + "=" * 92)

# ============================================================================
if __name__ == "__main__":
    banner("THREAD D: the relation-code LP HIERARCHY for measS7 (level-1 -> level-2)")
    print("Per shape E we report:")
    print("  measS7 = actual P(M=emptyset)               (the truth)")
    print("  L1bnd  = level-1 LP optimum (singleton-marginal-only, joint-law LP)")
    print("  L2bnd  = level-2 LP optimum (singleton + PAIR marginals fixed)")
    print("A bound IMPROVES (collapses->improves) at level 2 iff L2bnd < L1bnd strictly.")
    print("The hierarchy PINS the shape iff L2bnd == measS7 (LP optimum equals truth).")

    for k in [8, 9]:
        C = consec(k)
        measC = joint_occupancy(C)[0][0]
        L1C = level1_marginal_LP_opt(C)
        L2C = level2_LP_opt(C)
        print(f"\n--- k={k}  CONSEC E={C} ---")
        print(f"  measS7(consec) = {float(measC):.6f}")
        print(f"  L1bnd(consec)  = {L1C:.6f}   gap1 = {L1C-float(measC):+.6f}")
        print(f"  L2bnd(consec)  = {L2C:.6f}   gap2 = {L2C-float(measC):+.6f}")
        print(f"  level-2 IMPROVES over level-1: {L2C < L1C - 1e-9}")
        print(f"  level-2 PINS consec (L2bnd==measS7): {abs(L2C-float(measC)) < 1e-9}")
        # also the THM-534 moment-LP bound for reference
        if k in DUAL_G:
            Ly = float(L_y_level1(joint_occupancy(C)[0], k))
            print(f"  [ref] THM-534 L_y(consec) = {Ly:.6f}  (the deg-{binding_R(k)} factorial-moment dual)")
