#!/usr/bin/env python3
"""
Complete opus THM-615's general case: M(2U u {w1,w2}) >= 1/12 for LOOSE U (mac-mini-2026-07-04-S37).
opus proved the HARD point-only AP case (U=c*{1..11}); handed off the loose case (M(U)>1/12, where
{g_E>=1/12} is an INTERVAL, more room). Folding identity (opus): M(2U u {w1,w2}) = max_t min(g_U(2t), Psi(t)),
Psi(t)=max(min(a,b), 1/2-max(a,b)), a=||w1 t||, b=||w2 t||. Non-extremity (Psi>=1/12) <=> NOT (one of a,b<1/12
AND other>5/12). Want: some t in {g_U(2t)>=1/12} is non-extremity.
 (1) verify: for LOOSE 11-runner U, M(2U u {w1,w2})>=1/12 for ALL odd w1,w2 (small).
 (2) measure U's hiding-interval width delta vs the extremity-box scale 1/(6 min(w1,w2)); test the dodge.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import numpy as np

def gcd_all(xs): return reduce(gcd, xs)
def nd(x):
    x = x % 1
    return x if x <= 1-x else 1-x
def M_exact(sp):
    vs = sorted(set(sp)); Q = set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1, len(vs)): Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    best = F(0)
    for q in Q:
        if q < 2: continue
        for a in range(1, q):
            m = min(nd(v*F(a,q)) for v in sp)
            if m > best: best = m
    return best
def MU(U):  # M(U) exact
    return M_exact(U)
def hiding_interval_width(U):
    """width of the largest interval of {t: g_U(t) >= 1/12} (in t)."""
    pts = {F(0), F(1)}
    for v in U:
        for k in range(2*v+1):
            for s in (1,-1):
                t = F(12*k+s, 12*v)   # ||v t||=1/12 crossings
                if 0 <= t <= 1: pts.add(t)
    pts = sorted(pts); best = F(0)
    for i in range(len(pts)-1):
        a, b = pts[i], pts[i+1]
        if min(nd(v*((a+b)/2)) for v in U) >= F(1,12): best = max(best, b-a)
    return best

if __name__ == "__main__":
    import sys
    def out(*a): print(*a); sys.stdout.flush()
    # LOOSE 11-runner U (M(U) > 1/12)
    Us = [list(range(1,11))+[12], list(range(1,11))+[13], [1,2,3,4,5,6,7,8,9,11,13],
          list(range(1,11))+[15], [1,2,3,4,5,6,7,8,10,11,12]]
    out("(1) LOOSE 11-runner U: M(U), hiding-interval width delta, and min M(2U u {w1,w2}) over odd w1,w2<=39:")
    out(f"{'U (first 6)':>22} {'M(U)':>9} {'delta(1/12-interval)':>20} {'min M(2Uu2odd)':>15} {'>=1/12?':>8}")
    for U in Us:
        U = sorted(set(U))[:11]
        if len(U) != 11: continue
        mu = MU(U); dU = hiding_interval_width(U)
        E = [2*u for u in U]
        worst = F(1)
        for w1 in range(1,40,2):
            for w2 in range(w1+2,40,2):
                if w1 in E or w2 in E: continue
                S = E+[w1,w2]
                if len(set(S)) != 13: continue
                m = M_exact(S)
                if m < worst: worst = m
        out(f"{str(U[:6]):>22} {float(mu):>9.5f} {str(dU)+' '+f'{float(dU):.4f}':>20} {float(worst):>15.5f} {str(worst>=F(1,12)):>8}")
    out("\n(2) the dodge: on a hiding-interval J (width delta) in t, the extremity set X={a<1/12,b>5/12}∪sym is")
    out("    a union of boxes of width <= 1/(6 min(w1,w2)); if delta exceeds the box+gap scale, J has a")
    out("    non-extremity point => M>=1/12. LARGE tighteners: boxes tiny, trivial. SMALL tighteners: need delta big.")
    out("    (opus proved the delta=0 point-only AP case; loose U has delta>0 => interval room.)")
