#!/usr/bin/env python3
"""covering_tight_test_kps_S128c50.py -- THE SHARP CONJECTURE: covering => M > 1/14 (not tight).
'covering' = every modulus q in {2..14} divides some speed (operational sieve def).
If covering families never get tight, the residual (primitive COVERING tight) is EMPTY.
Sample + minimize M over covering families; report the covering-family M-minimum."""
import sys, random
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
random.seed(1729)
LAM=F(1,14)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def Mgrid(V,Q=16000):
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
def is_covering(V):
    for q in range(2,15):
        if not any(v%q==0 for v in V): return False
    return True
def rand_covering():
    # ensure each q in 2..14 hit: greedily assign multiples, fill to 13 distinct
    for _ in range(2000):
        need=list(range(2,15)); V=set()
        random.shuffle(need)
        for q in need:
            if any(v%q==0 for v in V): continue
            m=q*random.randint(1,6)
            V.add(m)
        while len(V)<13:
            V.add(random.randint(2,90))
        V=sorted(V)[:13]
        if len(V)==13 and is_covering(V) and len(set(V))==13:
            return V
    return None
best=(2.0,None)
samples=0
for s in range(3000):
    V=rand_covering()
    if V is None: continue
    samples+=1
    m=Mgrid(V)
    if m<best[0]: best=(m,tuple(V))
# local descent from best
c,Vb=best; Vb=list(Vb)
for _ in range(200):
    moved=False
    for i in range(13):
        for dv in (-2,-1,1,2):
            W=sorted(Vb[:i]+[Vb[i]+dv]+Vb[i+1:])
            if len(set(W))==13 and all(v>0 for v in W) and is_covering(W):
                m=Mgrid(W)
                if m<c-1e-9: Vb,c=W,m; moved=True; break
        if moved: break
    if not moved: break
Mex=M_exact(Vb)
print("== covering-family M-minimum (%d covering samples + descent) =="%samples)
print("  min covering family V=%s"%Vb)
print("  EXACT M=%s=%.6f ; > 1/14? %s ; margin=%.5f"%(Mex,float(Mex),Mex>LAM,float(Mex)-1/14))
print("  any covering family TIGHT (M=1/14)?", Mex==LAM)
if Mex>LAM:
    print("  >>> covering => M > 1/14 SUPPORTED: min covering M = %.5f (%.2fx threshold)"%(float(Mex),float(Mex)*14))
    print("  >>> the residual (primitive COVERING tight) is EMPTY on the sample -- equality horn closes conditionally.")
else:
    print("  *** a covering tight family exists:", Vb, "-- residual nonempty ***")
print("DONE")
