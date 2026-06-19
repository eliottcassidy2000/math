#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL part 5: nail down the wide-set maximum near consec-plus-outlier,
and confirm NOTHING wide reaches cap. Exact.
"""
import itertools
from fractions import Fraction as F
from functools import reduce
from math import gcd

def N_at(E,x):
    hit=set()
    for e in E:
        v=e*x; v=v-(v.numerator//v.denominator)
        hit.add((v.numerator*7)//v.denominator)
    return sum(1 for j in range(1,7) if j not in hit)

def dist_p(E):
    E=sorted(set(E))
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(F(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1)
    p=[F(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        p[N_at(E,(lo+hi)/2)]+=(hi-lo)
    return p

def g_poly(k):
    return [F((t-1)*(t-2)*(t-4)*(t-5),40) for t in range(7)]

def L_y(E,k=8):
    p=dist_p(E); g=g_poly(k)
    return sum(p[t]*g[t] for t in range(7))

# exact value of the found wide set
E=[0,1,2,3,4,5,6,17]
L=L_y(E)
print(f"E={E}: L_y={L}={float(L):.6f} spread={max(E)}")

# systematically: consec block of 7 + one outlier d, for d=8..40
print("\nconsec-7 + one outlier (0..6, then d):")
best=F(0); bestd=None
for d in range(8,41):
    E=[0,1,2,3,4,5,6,d]
    if reduce(gcd,E)!=1: continue
    L=L_y(E)
    if L>best: best=L; bestd=d
    if d<=25: print(f"  d={d}: L_y={float(L):.6f} spread={d}")
print(f"  best over d in 8..40: d={bestd} L_y={float(best):.6f}")

# consec-6 + two outliers, push spread up
print("\nconsec-6 (0..5) + two outliers:")
best2=F(0); bestpair=None
for d1 in range(6,30):
    for d2 in range(d1+1,40):
        E=[0,1,2,3,4,5,d1,d2]
        if reduce(gcd,E)!=1: continue
        L=L_y(E)
        if L>best2: best2=L; bestpair=(d1,d2)
print(f"  best: {bestpair} L_y={float(best2):.6f}")

# cap for k=8
print(f"\ncap_8=0.3815.  Any wide config above approaches but stays under? max found = {float(max(best,best2)):.6f}")
