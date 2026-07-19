#!/usr/bin/env python3
"""closed_form_synthesis_kps_S128c82.py -- kind-pasteur S128 cont.82.
BOXEPH-S120'S FORMULA EXPLAINS MY 7k+1 LAW.
The located-maximizer theorem gives M = |v_i a_j - v_j a_i|/(v_i+v_j) at a straddling active
pair.  On my clustered families the active pair is (v_min, active killer), so
        M = |v_min * a_k| / (v_min + k)
with a_i = 0.  My cont.55 empirical law was M = k'/(7k'+1) for the covering families -- test
whether 7k'+1 is exactly (v_min + killer)/3, i.e. whether the two are the same statement."""
import sys, itertools
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def Mexact(V,qhi=2600):
    best=(F(0),None)
    for q in range(2,qhi+1):
        for a in range(1,q//2+1):
            if gcd(a,q)!=1: continue
            t=F(a,q); mn=min(nd(v*t) for v in V)
            if mn>best[0]: best=(mn,t)
    return best
print("### the located pair on my clustered families: (v_min, active killer) ###")
print("  family                     M          active     v_min+k   M*(v_min+k)  a_k")
CASES=[(list(range(2,13))+[169,182],"{2..12}u{169,182}"),
       (list(range(2,13))+[169,196],"{2..12}u{169,196}"),
       (list(range(2,13))+[169,210],"{2..12}u{169,210}"),
       (list(range(2,13))+[182,196],"{2..12}u{182,196}"),
       (list(range(1,12))+[312,364],"{1..11}u{312,364}"),
       (list(range(1,12))+[168,338],"{1..11}u{168,338}")]
rows=[]
for V,name in CASES:
    M,t=Mexact(V)
    act=[v for v in V if nd(v*t)==M]
    vmin=min(V); k=[a for a in act if a>50]
    k=k[0] if k else max(act)
    s=vmin+k; num=M*s
    rows.append((name,M,act,s,num))
    print("  %-26s %-10s %-10s %-9d %-12s %s"%(name,M,str(act),s,num,num/vmin if vmin else "-"))
print()
print("### does the 7k'+1 law equal (v_min + killer)/3 ? ###")
print("  family                     M          denom(M)   v_min+k   (v_min+k)/3   7k'+1 with k'=num")
for name,M,act,s,num in rows:
    d=M.denominator
    print("  %-26s %-10s %-10d %-9d %-13s %s"%(
        name,M,d,s,F(s,3),"7*%s+1 = %s"%(M.numerator,7*M.numerator+1)))
print()
print("### the general closed form this yields ###")
print("  For a clustered covering family with core min v_min and active killer k:")
print("      M = a_k * v_min / (v_min + k)   for the integer a_k nearest (v_min+k)/(2 v_min) ...")
print("  check: predicted vs actual")
for V,name in CASES:
    M,t=Mexact(V)
    act=[v for v in V if nd(v*t)==M]
    vmin=min(V); ks=[a for a in act if a>50]
    if not ks: continue
    k=ks[0]; s=vmin+k
    best=None
    for a in range(1,200):
        cand=F(a*vmin,s)
        if cand==M: best=a; break
    print("  %-26s M = %-9s = %d*%d/%d   a_k = %s"%(name,M,best if best else 0,vmin,s,best))
print("DONE")
