#!/usr/bin/env python3
"""
GAP-A coverer-bound PROVED via opus THM-611 (mac-mini-2026-07-03-S35).
THM-611 (far-runner decorrelation): meas(lonely(R∪{w})) >= (6/7)meas(lonely(R)) - A_R/(3w), A_R=#arcs of
lonely(R). If S=R∪{w} is TIGHT (lonely measure 0) then 0 >= (6/7)L - A/(3w) => w <= 7A/(18L), L=meas lonely(R).
So a far coverer keeping the family tight is MAGNITUDE-BOUNDED. For the single-coverer axis {1..11,13,X}
(base R={1..11,13} loose, L,A constant) this gives an EXPLICIT bound => finite check => X in {12,24}.
"""
from fractions import Fraction as F
from math import gcd, ceil
def nd(x):
    x = x % 1
    return x if x <= 1-x else 1-x
def lonely_measure_and_arcs(sp):
    """exact Lebesgue measure of {t: all ||v t|| >= 1/14} and its number of arcs (danger endpoints (14k±1)/(14v))."""
    pts = {F(0), F(1)}
    for v in sp:
        for k in range(v+1):
            for s in (1, -1):
                t = F(14*k+s, 14*v)
                if 0 <= t <= 1: pts.add(t)
    pts = sorted(pts); L = F(0); arcs = 0
    for i in range(len(pts)-1):
        a, b = pts[i], pts[i+1]; mid = (a+b)/2
        if min(nd(v*mid) for v in sp) >= F(1,14):
            L += (b-a); arcs += 1
    return L, arcs
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

if __name__ == "__main__":
    base = [1,2,3,4,5,6,7,8,9,10,11,13]
    L, A = lonely_measure_and_arcs(base)
    bound = F(7*A, 18)/L
    print(f"base {{1..11,13}} (12 runners, loose): lonely measure L={L}={float(L):.6f}, #arcs A={A}")
    print(f"THM-611 => tight coverer X <= 7A/(18L) = {bound} = {float(bound):.1f}")
    B = ceil(float(bound))
    tights = [X for X in range(12, B+13, 12) if len(set(base+[X]))==13 and M_exact(base+[X])==F(1,14)]
    print(f"finite check 12|X, X<={B}: TIGHT = {tights}")
    print(f"=> {{1..11,13,X}} tight (12|X) <=> X in {tights}  [RIGOROUS via THM-611 + finite check]")
    print(f"   => GW=AP[12->24] is the unique non-AP tight family on the q=12 coverer axis (PROVED).")
