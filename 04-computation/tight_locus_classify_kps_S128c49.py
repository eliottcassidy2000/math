#!/usr/bin/env python3
"""tight_locus_classify_kps_S128c49.py -- classify parametrized tight-candidate families by
escape route (Layer1 dilation / Layer2 sieve / small / comparable / residual). Exact."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def M_exact(V):
    cand=set()
    for v in V:
        for a in range(1,2*v): cand.add(F(a,2*v))
    for i in range(len(V)):
        for j in range(i+1,len(V)):
            for s in (V[i]+V[j],abs(V[i]-V[j])):
                if s:
                    for a in range(1,s): cand.add(F(a,s))
    b=F(0)
    for t in cand:
        if 0<t<1:
            m=min(nd(v*t) for v in V)
            if m>b: b=m
    return b
def route(V):
    A=[abs(v) for v in V]
    if reduce(gcd,A)>=2: return "non-primitive(Layer1:common-res)"
    if max(A)<23: return "small(max<23)"
    if not any(a>13*b for a in A for b in A): return "comparable(no gap)"
    for q in range(2,15):
        if min(nd(v*F(1,q)) for v in V)>=F(1,14): return "sieve-covered(Layer2:q=%d)"%q
    return "*** PRIMITIVE COVERING GAP max>=23 (residual!) ***"
if __name__=="__main__":
    fams=[]
    for c in [24,25,27,36,48]: fams.append(("{1..11,13,%d}"%c, list(range(1,12))+[13,c]))
    for mult in [84,98,182,196]: fams.append(("{1..12,%d}"%mult, list(range(1,13))+[mult]))
    fams.append(("{1..13}", list(range(1,14)))); fams.append(("{1..11,14,27}", list(range(1,12))+[14,27]))
    for name,V in fams:
        M=M_exact(V); r=route(V)
        print("  %-18s %s  -> %s"%(name, "TIGHT" if M==F(1,14) else "M=%s"%M, r))
