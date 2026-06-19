#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) — does the far-limit ever exceed cap_k?  The decisive comparison.
kind-pasteur-2026-06-19-S12.

The crude bound LIM(B0) <= M_{k-1} + 1/7 is too loose at k=9 (gives 0.470,
exceeds M_9 but NOT cap_9=0.494). The ACTUAL far-limit is
   LIM(B0) = meas(S7)(B0) + (1/7)*P1(B0),  P1=P(exactly 1 sector missed).
We maximize LIM(B0) over ALL bounded (k-1)-bases B0 and compare to cap_k AND to
the global bounded max M_k.

If max_B0 LIM(B0) < M_k, then asymptotically every far element is dominated by
the bounded finite-check max -> reduction closes (with an Erdos-Turan rate for
the finite-w tail).  Even if max LIM(B0) is between M_k and cap_k, that still
shows far elements never violate the cap in the limit.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def dist_missed_count(E):
    E = sorted(set(E))
    bps = set([Fraction(0), Fraction(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7*e + 1):
            bps.add(Fraction(a, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [Fraction(0)]*7
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo+hi)/2
        hit = set()
        for e in E:
            v = e*mid; v = v - (v.numerator//v.denominator)
            hit.add((v.numerator*7)//v.denominator)
        t = sum(1 for j in range(1,7) if j not in hit)
        p[t] += (hi-lo)
    return p

def primitive(E):
    return reduce(gcd, [e for e in E if e != 0]) == 1

def max_LIM_over_bases(km1, search_bound):
    """max over (km1)-bases B0 of LIM(B0)=p0 + p1/7, and also report max meas_S7(B0)."""
    best_lim = Fraction(0); arg_lim = None
    best_m = Fraction(0); arg_m = None
    for combo in itertools.combinations(range(1, search_bound+1), km1-1):
        B0 = [0]+list(combo)
        if not primitive(B0): continue
        p = dist_missed_count(B0)
        lim = p[0] + p[1]*Fraction(1,7)
        if lim > best_lim:
            best_lim = lim; arg_lim = B0
        if p[0] > best_m:
            best_m = p[0]; arg_m = B0
    return best_lim, arg_lim, best_m, arg_m

if __name__ == "__main__":
    print("=== max far-limit LIM(B0) vs M_k and cap_k ===\n")
    caps = {8:Fraction(38153,100000), 9:Fraction(49426,100000), 10:Fraction(6044,10000)}
    SB = {7:11, 8:13, 9:13}
    Mk_map = {8:Fraction(0), 9:Fraction(0)}
    for k in [8, 9]:
        km1 = k-1
        sb = SB[km1]
        best_lim, arg_lim, best_m_km1, arg_m_km1 = max_LIM_over_bases(km1, sb)
        # global bounded max M_k = consec value (known); compute it
        from itertools import combinations
        # M_k via consec (global max, proven elsewhere)
        consec = list(range(k))
        # meas of consec:
        pc = dist_missed_count(consec); Mk = pc[0]
        cap = caps[k]
        print(f"--- k={k} ---")
        print(f"  max far-limit LIM(B0) over (k-1)-bases = {float(best_lim):.5f}  at B0={arg_lim}")
        print(f"  M_{k} (consec bounded global max)        = {float(Mk):.5f}")
        print(f"  cap_{k}                                  = {float(cap):.5f}")
        d_Mk = float(Mk - best_lim)
        d_cap = float(cap - best_lim)
        print(f"  M_{k}  - maxLIM = {d_Mk:+.5f}   {'(limit BELOW global max: reduction closes!)' if d_Mk>0 else '(limit ABOVE M_k)'}")
        print(f"  cap_{k}- maxLIM = {d_cap:+.5f}   {'(limit BELOW cap: far elt never violates cap)' if d_cap>0 else '(LIMIT EXCEEDS CAP!)'}")
        print()
