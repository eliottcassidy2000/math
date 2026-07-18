#!/usr/bin/env python3
"""margin_floor_hard_kps_S128c49.py -- push HARDER: can any trapped family beat M = 1/7?
Seed from the 1/7 minimizer + low-M samples; aggressive descent; larger moves; 200 starts."""
import sys, random
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
random.seed(2718)
def Mf(V, Q=30000):
    best=0.0; invQ=1.0/Q
    for p in range(1,Q):
        t=p*invQ; m=2.0
        for v in V:
            x=v*t; d=abs(x-round(x))
            if d<m:
                m=d
                if m<=best: break
        if m>best: best=m
    return best
def nd(x):
    fx=x-int(x)
    if fx<0: fx+=1
    return min(fx,1-fx)
def M_exact(V):
    cand=set(); n=len(V)
    for v in V:
        for a in range(1,2*v): cand.add(F(a,2*v))
    for i in range(n):
        for j in range(i+1,n):
            for s in (V[i]+V[j],abs(V[i]-V[j])):
                if s==0: continue
                for a in range(1,s): cand.add(F(a,s))
    best=F(0)
    for t in cand:
        if 0<t<1:
            m=min(nd(v*t) for v in V)
            if m>best: best=m
    return best
def gap(V): A=V; return any(a>13*b for a in A for b in A)
def compressed(V): return all(any(j!=i and V[i]<=13*V[j] for j in range(13)) for i in range(13))
def distinct(V): return len(set(V))==13
def noncluster(V):
    for i0 in range(13):
        g=reduce(gcd,[V[j] for j in range(13) if j!=i0],0)
        if g>=2 and V[i0]%g!=0: return False
    return True
def covering(V):
    for q in range(2,15):
        if min(nd(v*F(1,q)) for v in V)>=F(1,14): return False
    return True
def comm_res(V):
    for d in range(2,max(V)+1):
        for a in range(d):
            if all((v-a)%d==0 for v in V): return False
    return True
def trapped(V):
    return (distinct(V) and all(v>0 for v in V) and gap(V) and compressed(V)
            and max(V)>=23 and noncluster(V) and covering(V) and comm_res(V))
def descend(V):
    V=sorted(V); cur=Mf(V)
    for _ in range(400):
        moved=False
        order=list(range(13)); random.shuffle(order)
        for i in order:
            for dv in (-3,-2,-1,1,2,3):
                W=sorted(V[:i]+[V[i]+dv]+V[i+1:])
                if not trapped(W): continue
                mW=Mf(W)
                if mW<cur-1e-9: V,cur=W,mW; moved=True; break
            if moved: break
        if not moved: break
    return V,cur
seeds=[[42,66,96,108,150,228,229,247,375,377,396,414,552],
       [21,35,66,72,91,93,95,141,187,210,219,230,291],
       [21,38,78,116,133,140,144,187,201,206,216,225,277]]
best=(2.0,None)
for sd in seeds:
    if trapped(sd):
        V,c=descend(sd)
        if c<best[0]: best=(c,tuple(V))
for s in range(200):
    lo=random.choice([18,22,26,32,40,55,70])
    V=None
    for _ in range(3000):
        cand=sorted(random.sample(range(lo,lo*15),13))
        if trapped(cand): V=cand; break
    if V is None: continue
    V,c=descend(V)
    if c<best[0]:
        best=(c,tuple(V)); print("  new min ~ %.5f"%c)
c,Vb=best; Vb=list(Vb)
Mex=M_exact(Vb)
print("\nHARDEST trapped minimum: V=%s"%Vb)
print("  EXACT M = %s = %.6f ; margin over 1/14 = %.6f ; trapped=%s"%(Mex,float(Mex),float(Mex)-1/14,trapped(Vb)))
print("  M as multiple of 1/14: %.4f x"%(float(Mex)*14))
print("  broke below 1/7?", float(Mex) < 1/7 - 1e-9)
print("DONE")
