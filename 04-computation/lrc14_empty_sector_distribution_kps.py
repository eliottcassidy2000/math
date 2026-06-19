#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The empty-sector count distribution p_t = meas{x: N(x)=t}, t=0..6, for LRC(14)-S3.
N(x) = #{ j in 1..6 : sector [j/7,(j+1)/7) is MISSED by the orbit {frac(e_i x): e in E} }.
meas(S7)=p_0. L_y(E)=E[g(N)] (THM-534 dual). HYP-2607: consec maximizes L_y.
kind-pasteur-2026-06-19-S11.  Goal: SEE why consec is extremal (convex-order/coupling intuition).

Exact via breakpoints x = a/(7 e_i) (where some e_i x crosses a sector wall j/7).
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def sector(xfrac):
    """which 1/7-sector (0..6) does a Fraction in [0,1) lie in."""
    return (xfrac.numerator * 7) // xfrac.denominator   # floor(7x)

def N_at(E, x):
    """number of sectors among {1..6} missed by {frac(e x): e in E} at point x (Fraction)."""
    hit = set()
    for e in E:
        v = e * x
        v = v - (v.numerator // v.denominator)   # frac
        s = (v.numerator * 7) // v.denominator
        hit.add(s)
    return sum(1 for j in range(1, 7) if j not in hit)

def dist_p(E):
    """exact distribution p_t (Fraction), t=0..6, of N over x in [0,1)."""
    E = sorted(set(E))
    # breakpoints: x = a/(7 e) for e in E (e>0), a=0..7e ; plus 0 and 1
    bps = set([Fraction(0), Fraction(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7 * e + 1):
            bps.add(Fraction(a, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [Fraction(0)] * 7
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo + hi) / 2
        t = N_at(E, mid)
        p[t] += (hi - lo)
    return p

def g_poly(k):
    """g(t) for t=0..6 (THM-534 integer-root dual), as Fractions, for cluster size k."""
    g = []
    for t in range(7):
        if k == 8:
            val = Fraction((t-1)*(t-2)*(t-4)*(t-5), 40)
        elif k in (9, 10):
            val = Fraction(-(t-2)*(t-3)*(t-6), 36)
        else:  # 11,12,13
            val = Fraction((t-3)*(t-4), 12)
        g.append(val)
    return g

def L_y(E, k):
    p = dist_p(E); g = g_poly(k)
    return sum(p[t] * g[t] for t in range(7)), p

def consec(k):  # 0 in E, |E|=k
    return list(range(k))

if __name__ == "__main__":
    print("=== distribution of N (p_0..p_6) and L_y for consec vs competitors ===\n")
    for k in [8, 9, 10]:
        g = g_poly(k)
        print(f"--- k={k}  g(t)={[str(x) for x in g]} (support {[t for t in range(7) if g[t]!=0]}) ---")
        C = consec(k)
        Lc, pc = L_y(C, k)
        print(f"  consec={C}")
        print(f"    p = {[f'{float(x):.4f}' for x in pc]}   (p0=meas S7={float(pc[0]):.4f})")
        print(f"    L_y(consec) = {Lc} = {float(Lc):.5f}")
        print(f"    E[N]={float(sum(t*pc[t] for t in range(7))):.4f}  Var[N]={float(sum(t*t*pc[t] for t in range(7))-(sum(t*pc[t] for t in range(7)))**2):.4f}")
        # a few competitors: small perturbations of consec
        comps = []
        if k == 8:
            comps = [[0,1,2,3,4,5,6,8],[0,1,2,3,4,5,7,8],[0,1,2,3,4,6,7,8],[0,2,3,4,5,6,7,9],[0,1,2,3,5,7,9,11]]
        elif k == 9:
            comps = [[0,1,2,3,4,5,6,7,9],[0,1,2,3,4,5,6,8,9],[0,1,2,3,4,5,7,8,10]]
        elif k == 10:
            comps = [[0,1,2,3,4,5,6,7,8,10],[0,1,2,3,4,5,6,7,9,10]]
        for Cc in comps:
            if not (reduce(gcd,Cc)==1): continue
            Lcc, pcc = L_y(Cc, k)
            mark = "  <-- beats consec!" if Lcc > Lc else ""
            print(f"    {Cc}: L_y={float(Lcc):.5f} p={[f'{float(x):.3f}' for x in pcc]}{mark}")
        print()
