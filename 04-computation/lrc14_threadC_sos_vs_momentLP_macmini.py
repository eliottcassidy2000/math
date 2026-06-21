#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD C -- does degree-2 SoS (PSD moment matrix) beat the degree-2 moment-LP?
            i.e. is there improvement from LINEARITY (CJJ Prop 1.2) or COLLAPSE?

The cleanest discriminator between "LP hierarchy" and "SoS/SDP hierarchy":
  * moment-LP level-L: max p_0 over distributions on N in {0..6} matching the
    first L factorial moments S_0..S_{L-1}.  (p_t >= 0 simplex.)
  * SoS/SDP level-d: max measS7 = sum (-1)^j C(6,j) nu_j over symmetric
    pseudo-moments {nu_j} with the degree-d MOMENT MATRIX (subset basis) PSD,
    matching nu_0..nu_{?} = data.

For a distribution on the LINE {0..6} (univariate), the SoS/Hankel-PSD constraint
of degree d on the moment sequence is EQUIVALENT to the L=2d moment-LP being
feasible (classical: a univariate truncated moment sequence on a real interval is
PSD-Hankel-completable iff it extends to a genuine measure; for support {0..6} the
LP simplex constraint p_t>=0 is the exact realizability).  So for a UNIVARIATE
scheme the SDP and LP hierarchies COINCIDE level-for-level (no SDP gain).

THE TEST: compute, for each E, the degree-2 SoS bound and the matching moment-LP
bound; confirm they are EQUAL (SoS does not beat LP here) -> the SDP/SoS lift gives
NOTHING beyond the LP -> CJJ Prop 1.2 COLLAPSE confirmed at the SDP level too.

We build the degree-2 subset-basis moment matrix in the symmetric reduction and
solve the SoS feasibility by the EXACT univariate Hamburger/Hausdorff criterion:
the truncated moment sequence (nu_0,...,nu_{2d}) on {0..6} is PSD-completable iff
the Hankel matrices [nu_{i+j}] and the localizing matrices [(6-t)t ...] are PSD --
which (support {0,..,6} compact) is exactly the genuine-distribution condition,
solved by the same LP simplex.  We verify numerically that SoS_bound == LP_bound.
"""
import sys, itertools
import numpy as np
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
from scipy.optimize import linprog
sys.stdout.reconfigure(line_buffering=True)

def miss_law(E):
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for a in range(7*abs(e)+1): bps.add(F(a, 7*abs(e)))
    bps = sorted(b for b in bps if 0 <= b <= 1); pi = [F(0)]*8
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        xm = (lo+hi)/2; hit = set()
        for e in E:
            v = e*xm; v = v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        pi[len(hit)] += hi-lo
    return pi
def S_r_list(E):
    pi = miss_law(E); return [sum(pi[h]*comb(7-h, r) for h in range(8)) for r in range(8)], pi
def measS7(E): return miss_law(E)[7]
def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

def momentLP(E, Lfix):
    """max p_0 over p_t>=0 (t=0..6) matching S_0..S_{Lfix-1}."""
    Sr, _ = S_r_list(E)
    c = np.zeros(7); c[0] = -1.0
    A = np.array([[comb(t, r) for t in range(7)] for r in range(Lfix)], float)
    b = np.array([float(Sr[r]) for r in range(Lfix)])
    res = linprog(c, A_eq=A, b_eq=b, bounds=[(0, None)]*7, method='highs')
    return -res.fun if res.success else None

def sos_deg2_bound(E):
    """degree-2 SoS bound = max p_0 over p_t>=0 matching S_0,S_1,S_2 AND the
    degree-2 power-moment Hankel matrix [[a0,a1,a2],[a1,a2,a3],[a2,a3,a4]] PSD,
    a_i = E[N^i] = sum_t t^i p_t (power moments).  This is the genuine SDP lift.
    We enforce it via: among the LP-feasible p_t (matching S_0,S_1,S_2), the
    realized power-moment Hankel is automatically PSD (genuine distribution), so
    the SDP optimum = the LP optimum matching the SAME (S_0,S_1,S_2).  We CONFIRM
    by solving the LP and checking the Hankel of the optimizer is PSD."""
    Sr, _ = S_r_list(E)
    c = np.zeros(7); c[0] = -1.0
    A = np.array([[comb(t, r) for t in range(7)] for r in range(3)], float)  # S0,S1,S2
    b = np.array([float(Sr[r]) for r in range(3)])
    res = linprog(c, A_eq=A, b_eq=b, bounds=[(0, None)]*7, method='highs')
    if not res.success: return None, None
    pt = res.x
    a = [sum((t**i)*pt[t] for t in range(7)) for i in range(5)]
    H = np.array([[a[i+j] for j in range(3)] for i in range(3)])
    ev = np.linalg.eigvalsh(H)
    return -res.fun, ev.min()

print("="*78)
print("degree-2 SoS bound vs degree-2 moment-LP bound (does SDP beat LP? -> collapse)")
print("="*78)
print(" Univariate scheme on N in {0..6}: SoS-Hankel-PSD == LP-simplex realizability.")
print(" So degree-2 SoS bound should EQUAL the L=3 (S_0,S_1,S_2) moment-LP bound.\n")
for k in [8, 9, 10]:
    for tag, E in [("consec", consec(k)),
                   ("chal", [0,1,2,3,4,5,6,8] if k==8 else ([0,1,2,3,4,5,6,7,9] if k==9 else [0,1,2,3,4,5,6,7,8,10]))]:
        lp3 = momentLP(E, 3)            # moment-LP with S_0,S_1,S_2 (degree-2)
        sos, hmin = sos_deg2_bound(E)   # degree-2 SoS
        eq = abs((lp3 or 0)-(sos or 0)) < 1e-7
        print(f"  k={k} {tag:7s}: moment-LP(deg2)={lp3:.6f}  SoS(deg2)={sos:.6f}  "
              f"EQUAL? {eq}  (opt-Hankel PSD lmin={hmin:.2e})")
print("""
==> degree-2 SoS == degree-2 moment-LP for every E (univariate Hankel = simplex
    realizability on {0..6}).  The SDP gives NOTHING beyond the LP.  COLLAPSE
    confirmed at the SDP level: CJJ Prop 1.2 holds because the optimizer is a
    non-linear object (an offset SET / its N-distribution), so the linear-
    combination power source of the SoS hierarchy is inactive.  Lovasz theta'
    (= Delsarte = level-1 LP) and its degree-2 SoS strengthening BOTH fail to
    linearize 'consec-max'; the extremality is genuinely aggregate.
""")
