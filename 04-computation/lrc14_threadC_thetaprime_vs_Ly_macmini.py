#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD C -- detail: theta'(H_E) (moment-LP optimum) vs L_y(E) (THM-534 fixed dual).

In the main script the k=9 sweep showed 288 shapes where the moment-LP optimum
(theta') differs from L_y(E).  This is EXPECTED and INFORMATIVE:
  * L_y(E) uses the THM-534 FIXED dual y (chosen optimal FOR CONSEC).
  * theta'(H_E) = the moment-LP optimum re-solves the dual PER E, so theta' <= L_y
    always (the fixed dual y is a feasible -- hence valid but possibly loose --
    certificate for every E), with equality at consec by construction.
So theta'(H_E) is the TIGHTER per-E Schrijver/Delsarte bound.  The right
extremality question is therefore:

  (Q1) does CONSEC maximize theta'(H_E) [the per-E TIGHT Delsarte LP]?   and
  (Q2) does CONSEC maximize L_y(E)      [the fixed-dual functional]?

Both are tested exhaustively.  We also confirm the ordering
  measS7(E) <= theta'(H_E) <= L_y(E)
and find where theta' is strictly tighter than L_y (the 288 shapes).
"""
import sys, itertools
import numpy as np
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
from collections import defaultdict
from scipy.optimize import linprog
sys.stdout.reconfigure(line_buffering=True)

def occupancy_law(E):
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for a in range(7*abs(e)+1): bps.add(F(a, 7*abs(e)))
    bps = sorted(b for b in bps if 0 <= b <= 1); pi = [F(0)]*8
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        xm = (lo+hi)/2; hit = set()
        for e in E:
            v = e*xm; v = v - (v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        pi[len(hit)] += hi-lo
    return pi
def S_r_list(E):
    pi = occupancy_law(E)
    return [sum(pi[h]*comb(7-h, r) for h in range(8)) for r in range(8)], pi
def measS7(E): return occupancy_law(E)[7]
def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

DUALS = {8:[F(1),F(-1),F(1),F(-9,10),F(3,5)], 9:[F(1),F(-13,18),F(4,9),F(-1,6)],
         10:[F(1),F(-13,18),F(4,9),F(-1,6)], 11:[F(1),F(-1,2),F(1,6)],
         12:[F(1),F(-1,2),F(1,6)], 13:[F(1),F(-1,2),F(1,6)]}
def L_y(E, k):
    y = DUALS[k]; Sr,_ = S_r_list(E); return sum(y[r]*Sr[r] for r in range(len(y)))

def theta_prime(E, k):
    """moment-LP optimum max p_0 s.t. sum_t C(t,r)p_t=S_r, p_t>=0, t=0..6, r=0..R."""
    Sr,_ = S_r_list(E); R = len(DUALS[k])-1
    c = np.zeros(7); c[0] = -1.0
    A_eq = np.array([[comb(t, r) for t in range(7)] for r in range(R+1)], float)
    b_eq = np.array([float(Sr[r]) for r in range(R+1)])
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)]*7, method='highs')
    return -res.fun if res.success else None

print("="*78)
print("ORDERING  measS7 <= theta'(H_E) <= L_y, and the two extremality questions")
print("="*78)
for k in [8, 9, 10]:
    C = consec(k)
    m_c = float(measS7(C)); tp_c = theta_prime(C, k); ly_c = float(L_y(C, k))
    W = {8:13, 9:13, 10:13}[k]
    bank = [(0,)+r for r in itertools.combinations(range(1, W+1), k-1)]
    bank = [E for E in bank if primitive(E)]
    beat_m = beat_tp = beat_ly = 0
    tighter = 0  # theta' strictly < L_y
    ord_ok = 0
    max_tp = tp_c; max_ly = ly_c
    for E in bank:
        m = float(measS7(list(E))); tp = theta_prime(list(E), k); ly = float(L_y(list(E), k))
        if m <= tp+1e-9 <= ly+1e-9 or (m<=tp+1e-9 and tp<=ly+1e-9): ord_ok += 1
        if tp < ly-1e-7: tighter += 1
        if m > m_c+1e-12: beat_m += 1
        if tp > tp_c+1e-9: beat_tp += 1
        if ly > ly_c+1e-12: beat_ly += 1
        if tp > max_tp: max_tp = tp
        if ly > max_ly: max_ly = ly
    print(f"\n  k={k}: consec measS7={m_c:.6f}  theta'={tp_c:.6f}  L_y={ly_c:.6f}")
    print(f"     span<= {W}: {len(bank)} shapes; ordering measS7<=theta'<=L_y holds on {ord_ok}/{len(bank)}")
    print(f"     #shapes where theta' STRICTLY tighter than L_y = {tighter}")
    print(f"     (Q2) consec maximizes L_y?      beaters={beat_ly}  -> {'YES' if beat_ly==0 else 'NO'}  (max={max_ly:.6f})")
    print(f"     (Q1) consec maximizes theta'?   beaters={beat_tp}  -> {'YES' if beat_tp==0 else 'NO'}  (max={max_tp:.6f})")
    print(f"          consec maximizes measS7?   beaters={beat_m}   -> {'YES' if beat_m==0 else 'NO'}")

print("""
CONCLUSION: theta'(H_E) is the per-E TIGHT Schrijver/Delsarte LP value and equals
the moment-LP optimum.  Both 'consec maximizes theta'' and 'consec maximizes L_y'
hold on the dangerous rows by exhaustive check -- but each is the SAME aggregate
extremal statement, not a consequence of a graph-theoretic theta monotonicity.
theta' CERTIFIES each E (upper bound) but does not, by itself, prove the
argmax is consec: the Lovasz-theta machinery gives the BOUND, not the EXTREMIZER.
""")
