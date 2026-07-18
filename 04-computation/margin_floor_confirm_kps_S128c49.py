#!/usr/bin/env python3
"""margin_floor_confirm_kps_S128c49.py -- light re-confirm: 25 independent starts, coarse grid,
does the trapped minimizer beat 1/7?  (independent seed from cont.49 v1)"""
import sys, random
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
random.seed(112358)
def Mf(V,Q=14000):
    best=0.0; iq=1.0/Q
    for p in range(1,Q):
        t=p*iq; m=2.0
        for v in V:
            d=abs(v*t-round(v*t))
            if d<m:
                m=d
                if m<=best: break
        if m>best: best=m
    return best
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def M_exact(V):
    cand=set()
    for v in V:
        for a in range(1,2*v): cand.add(F(a,2*v))
    for i in range(13):
        for j in range(i+1,13):
            for s in (V[i]+V[j],abs(V[i]-V[j])):
                if s:
                    for a in range(1,s): cand.add(F(a,s))
    b=F(0)
    for t in cand:
        if 0<t<1:
            m=min(nd(v*t) for v in V)
            if m>b: b=m
    return b
def trapped(V):
    if len(set(V))!=13 or any(v<=0 for v in V) or max(V)<23: return False
    if not any(a>13*b for a in V for b in V): return False
    if not all(any(j!=i and V[i]<=13*V[j] for j in range(13)) for i in range(13)): return False
    for i0 in range(13):
        g=reduce(gcd,[V[j] for j in range(13) if j!=i0],0)
        if g>=2 and V[i0]%g!=0: return False
    for q in range(2,15):
        if min(nd(v*F(1,q)) for v in V)>=F(1,14): return False
    for d in range(2,max(V)+1):
        for a in range(d):
            if all((v-a)%d==0 for v in V): return False
    return True
best=(2.0,None)
for s in range(25):
    lo=random.choice([20,28,38,50])
    V=None
    for _ in range(3000):
        c=sorted(random.sample(range(lo,lo*15),13))
        if trapped(c): V=c; break
    if not V: continue
    cur=Mf(V)
    for _ in range(150):
        moved=False
        for i in range(13):
            for dv in (-2,-1,1,2):
                W=sorted(V[:i]+[V[i]+dv]+V[i+1:])
                if trapped(W):
                    m=Mf(W)
                    if m<cur-1e-9: V,cur=W,m; moved=True; break
            if moved: break
        if not moved: break
    if cur<best[0]: best=(cur,tuple(V))
c,Vb=best
Mex=M_exact(list(Vb))
print("confirm min: M=%s=%.6f (%.3fx of 1/14) V=%s"%(Mex,float(Mex),float(Mex)*14,list(Vb)))
print("beat 1/7?", float(Mex)<1/7-1e-9)
print("DONE")
