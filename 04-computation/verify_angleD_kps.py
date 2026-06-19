#!/usr/bin/env python3
# Independent adversarial verification of ANGLE D claims (LRC(14) endgame).
# Own breakpoint engine, Fraction exact.
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd

def N_at(E, x):
    """count sectors j in 1..6 with NO orbit point of {frac(e*x)} in [j/7,(j+1)/7)."""
    hit = set()
    for e in E:
        v = e * x
        v = v - (v.numerator // v.denominator)  # frac in [0,1)
        s = (v.numerator * 7) // v.denominator   # floor(7*frac)
        hit.add(s)
    return sum(1 for j in range(1, 7) if j not in hit)

def maxgap(E, x):
    """max circular gap between consecutive orbit points (as Fraction)."""
    pts = set()
    for e in E:
        v = e * x
        v = v - (v.numerator // v.denominator)
        pts.add(v)
    pts = sorted(pts)
    if len(pts) == 1:
        return Fraction(1)
    g = Fraction(0)
    for i in range(len(pts)):
        nxt = pts[(i+1) % len(pts)]
        d = (nxt - pts[i]) % 1 if i < len(pts)-1 else (pts[0] + 1 - pts[i])
        if i < len(pts)-1:
            d = pts[i+1] - pts[i]
        else:
            d = pts[0] + 1 - pts[i]
        if d > g: g = d
    return g

def breakpoints(E):
    bps = set([Fraction(0), Fraction(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7*e + 1):
            bps.add(Fraction(a, 7*e))
    return sorted(b for b in bps if 0 <= b <= 1)

def dist_p(E):
    E = sorted(set(E))
    bps = breakpoints(E)
    p = [Fraction(0)]*7
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo+hi)/2
        p[N_at(E, mid)] += (hi-lo)
    return p

def meas_maxgap_small(E):
    """meas{x : maxgap < 1/7}."""
    E = sorted(set(E))
    bps = breakpoints(E)
    m = Fraction(0)
    one7 = Fraction(1,7)
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid=(lo+hi)/2
        if maxgap(E, mid) < one7:
            m += (hi-lo)
    return m

def consec(k): return list(range(k))

print("=== (1) consec distributions ===")
expected = {
 8: ["481/1470","359/1470","25/147","26/245","17/210","5/98","1/49"],
 9: ["2447/5880","653/2940","27/196","23/245","13/196","9/196","1/56"],
 10:["8899/17640","551/2940","127/1176","145/1764","55/882","5/126","1/63"],
}
for k in [8,9,10]:
    p = dist_p(consec(k))
    print(f"k={k}: p = {[str(x) for x in p]}")
    exp = [Fraction(s) for s in expected[k]]
    ok = (p == exp)
    print(f"   matches claimed: {ok}  ; sum={sum(p)} ; p0={p[0]}={float(p[0]):.6f}")
    if not ok:
        for t in range(7):
            if p[t]!=exp[t]: print(f"     MISMATCH t={t}: got {p[t]} exp {exp[t]}")

print("\n=== (2) dilated AP {0,2,...,2(k-1)} identical to consec ===")
for k in [8,9,10]:
    pc = dist_p(consec(k))
    pd = dist_p([2*i for i in range(k)])
    print(f"k={k}: identical = {pc==pd}")

print("\n=== (4) maxgap<1/7 => N=0, and meas{maxgap<1/7} ===")
exp_mg = {8:"44/735", 9:"47/294", 10:"11/49"}
for k in [8,9,10]:
    m = meas_maxgap_small(consec(k))
    p0 = dist_p(consec(k))[0]
    print(f"k={k}: meas(maxgap<1/7)={m}={float(m):.6f}  claimed {exp_mg[k]}={float(Fraction(exp_mg[k])):.6f}  match={m==Fraction(exp_mg[k])}  <=p0? {m<=p0}")

print("\n=== (3) consec uniquely maximizes p_0 over bounded-spread (k=8, spread<=11) ===")
k=8
C=consec(k); p0c=dist_p(C)[0]
beat=0; tie=0; ties=[]
maxspread=11
# all sets with 0, |E|=k, max element <= maxspread
elems=list(range(1,maxspread+1))
cnt=0
for combo in itertools.combinations(elems, k-1):
    E=[0]+list(combo)
    if reduce(gcd,E)!=1: continue
    cnt+=1
    p0=dist_p(E)[0]
    if p0>p0c: beat+=1
    elif p0==p0c and E!=C: tie+=1; ties.append(E)
print(f"k=8 spread<=11: tested {cnt} sets (gcd=1), beat={beat}, tie={tie}")
print(f"   tie sets (should be only AP class {{0,d,2d,...}}): {ties}")
