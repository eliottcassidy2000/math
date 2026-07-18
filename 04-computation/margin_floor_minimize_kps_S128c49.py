#!/usr/bin/env python3
"""margin_floor_minimize_kps_S128c49.py -- kind-pasteur S128 cont.49.
MARGIN FLOOR: minimize M(V) over the TRAPPED cut via multi-start local descent.
If even the optimizer cannot drive M below 1/14 + delta, delta is the empirical margin floor
and 'trapped => M > 1/14' is stress-tested from the hardest direction.
Fast float M for search; exact rational M for the best families found.
"""
import sys, random
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
random.seed(31415)
LAMf = 1.0/14

def Mf(V, Q=24000):
    """fast float lower-estimate of M = max_t min_i ||vi t||."""
    best = 0.0
    invQ = 1.0/Q
    for p in range(1, Q):
        t = p*invQ
        m = 2.0
        for v in V:
            x = v*t
            d = abs(x - round(x))
            if d < m:
                m = d
                if m <= best: break
        if m > best: best = m
    return best

def nd(x):
    fx = x-int(x)
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

def gap(V): A=[abs(v) for v in V]; return any(a>13*b for a in A for b in A)
def compressed(V): A=[abs(v) for v in V]; return all(any(j!=i and A[i]<=13*A[j] for j in range(len(A))) for i in range(len(A)))
def distinct(V): return len(set(V))==len(V)
def maxge23(V): return max(V)>=23
def noncluster(V):
    n=len(V)
    for i0 in range(n):
        g=reduce(gcd,[V[j] for j in range(n) if j!=i0],0)
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
            and maxge23(V) and noncluster(V) and covering(V) and comm_res(V))

def rand_trapped(lo):
    for _ in range(4000):
        V=sorted(random.sample(range(lo, lo*15), 13))
        if trapped(V): return V
    return None

best_global=(2.0,None)
STARTS=60
for s in range(STARTS):
    lo=random.choice([20,25,30,40,50])
    V=rand_trapped(lo)
    if V is None: continue
    curM=Mf(V)
    improved=True
    while improved:
        improved=False
        for i in range(13):
            for dv in (-2,-1,1,2):
                W=V[:]; W[i]+=dv; W=sorted(W)
                if not trapped(W): continue
                mW=Mf(W)
                if mW < curM - 1e-9:
                    V,curM=W,mW; improved=True; break
            if improved: break
    if curM < best_global[0]:
        best_global=(curM, tuple(V))
        print("  start %2d: new min M ~ %.5f  V=%s"%(s,curM,list(V)))
print()
mf, Vbest = best_global
Vb = list(Vbest)
Mex = M_exact(Vb)
print("MINIMIZED trapped family: V = %s"%Vb)
print("  float M ~ %.6f ; EXACT M = %s = %.6f"%(mf, Mex, float(Mex)))
print("  margin over 1/14 = %s = %.6f"%(Mex-F(1,14), float(Mex)-1/14))
print("  trapped verified:", trapped(Vb), " primitive:", reduce(gcd,Vb)==1)
print()
if float(Mex) > 1/14:
    print(">>> Even minimized, the trapped family has M = %.5f > 1/14 (margin %.4f)."%(float(Mex), float(Mex)-1/14))
    print(">>> Empirical margin floor delta >~ %.4f on the trapped cut (multi-start descent, %d starts)."%(float(Mex)-1/14, STARTS))
else:
    print(">>> ALERT: minimizer drove a trapped family to M <= 1/14 -- investigate:", Vb, Mex)
print("DONE")
