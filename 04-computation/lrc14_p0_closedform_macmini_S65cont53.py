#!/usr/bin/env python3
"""cont.53: hunt p0(consec_k) = P(all 7 sectors hit) closed form. p0=0 for k<=6, positive
k>=7. Test structural forms; the three-gap orbit {0,x,..,(k-1)x} hits all 7 sectors."""
from fractions import Fraction as F
def p0_consec(k):
    E=list(range(k)); pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(x for x in pts if 0<=x<=1); p0=F(0)
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        if len(hit)==7: p0+=b-a
    return p0
print("p0(consec_k) exact + numerator over 7-smooth denominators:")
for k in range(7,18):
    p0=p0_consec(k)
    # denominator structure
    print(f"  k={k}: p0 = {p0} = {float(p0):.6f},  den={p0.denominator}")
print("\ntry p0 * 7*(k-1) (echo the p6/p5 denominators):")
for k in range(7,18):
    p0=p0_consec(k)
    v=p0*7*(k-1)
    print(f"  k={k} ({'even' if k%2==0 else 'odd '}): p0*7(k-1) = {v} = {float(v):.5f}")
print("\ntry 1-p0 = T1 = meas(S7), and T1*something:")
for k in range(7,18):
    T1=1-p0_consec(k)
    print(f"  k={k}: T1 = {T1} = {float(T1):.6f}, T1*7(k-1)={float(T1*7*(k-1)):.4f}")
