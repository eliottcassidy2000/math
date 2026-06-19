#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) — EXPLICIT B(k): finding the threshold past which a far element
provably keeps meas(S7) below the bounded finite-check maximum.
kind-pasteur-2026-06-19-S12.

Weyl limit (verified, lrc14_weyl_limit_kps.py):
   LIM(B0) = meas(S7)(B0) + (1/7)*P(B0 misses exactly 1 sector)
with O(1/w) convergence: |meas(S7)(B0 u {w}) - LIM(B0)| <= D(B0)/w (Erdos-Turan).

GOAL: For each k, define
   M_k     = max over all bounded k-sets of meas(S7)        (the global max = consec)
   M_{k-1} = max over all (k-1)-sets of meas(S7)
The far limit satisfies  LIM(B0) <= M_{k-1} + 1/7.
We check  M_{k-1} + 1/7  vs  M_k  vs  cap_k.

If M_{k-1} + 1/7 < M_k, then every far element strictly contracts below the
global max in the LIMIT, and an explicit Erdos-Turan rate gives B(k) at which
the finite-w correction can no longer reach M_k.

We compute exact maxima over bounded sets by exhaustive search up to a search
bound, and report the gap  M_k - (M_{k-1}+1/7).
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def meas_S7(E):
    E = sorted(set(E))
    bps = set([Fraction(0), Fraction(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7*e + 1):
            bps.add(Fraction(a, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p0 = Fraction(0)
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo+hi)/2
        hit = set()
        for e in E:
            v = e*mid; v = v - (v.numerator//v.denominator)
            hit.add((v.numerator*7)//v.denominator)
        if all(j in hit for j in range(1,7)):
            p0 += (hi-lo)
    return p0

def primitive(E):
    return reduce(gcd, [e for e in E if e != 0]) == 1

def max_meas_over_jsets(j, search_bound):
    """exact max of meas(S7) over primitive j-sets containing 0, max elt <= search_bound."""
    best = Fraction(0); arg = None
    for combo in itertools.combinations(range(1, search_bound+1), j-1):
        E = [0] + list(combo)
        if not primitive(E): continue
        m = meas_S7(E)
        if m > best:
            best = m; arg = E
    return best, arg

if __name__ == "__main__":
    print("=== EXPLICIT B(k): far-limit gap vs bounded global max ===\n")
    # search bounds chosen so consec is certainly the max (consec is global max).
    # j-set max search bound:
    SB = {5:9, 6:10, 7:11, 8:13, 9:14}
    caps = {8:0.38153, 9:0.49426, 10:0.6044}
    for k in [8, 9]:
        Mk, argk = max_meas_over_jsets(k, SB.get(k, k+5))
        Mk1, argk1 = max_meas_over_jsets(k-1, SB.get(k-1, k+4))
        far_lim_ub = Mk1 + Fraction(1,7)
        cap = caps[k]
        print(f"--- k={k} ---")
        print(f"  M_{k}   (bounded global max)   = {float(Mk):.5f}  at {argk}")
        print(f"  M_{k-1} (max over (k-1)-sets)  = {float(Mk1):.5f}  at {argk1}")
        print(f"  far-limit upper bound = M_{k-1}+1/7 = {float(far_lim_ub):.5f}")
        print(f"  cap_{k} = {cap}")
        gap = float(Mk - far_lim_ub)
        print(f"  GAP  M_{k} - (M_{k-1}+1/7) = {gap:+.5f}  "
              + ("(far limit BELOW global max -- contraction in limit)" if gap>0 else "(NO gap -- limit could match max!)"))
        print()
    print("INTERPRETATION:")
    print("  If GAP>0, every far element contracts strictly below the bounded global")
    print("  max in the w->inf limit. An Erdos-Turan rate |meas - LIM| <= D/w then")
    print("  gives explicit B(k) = ceil(D / GAP): for w>B(k), meas(S7) < M_k <= cap_k.")
