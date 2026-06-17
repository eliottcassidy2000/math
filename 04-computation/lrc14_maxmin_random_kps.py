#!/usr/bin/env python3
"""
lrc14_maxmin_random_kps  (kind-pasteur 2026-06-17) -- randomized counterexample hunt.
LRC(14): every 13-set has max-min >= 1/14.  Hunt for max-min < 1/14 over:
 - random 13-sets (various ranges)
 - "view-obstruction" style sets (multiples, geometric, mixed)
 - local optimization: start from random, greedily replace one element to LOWER max-min.
Use a fast accurate float max-min (dense grid + kink refinement); exact-verify the best.
"""
import sys, io, random
try: sys.stdout=io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
except Exception: pass
from fractions import Fraction as Fr
from functools import reduce
from math import gcd

def fn_f(x):
    x=x-int(x)
    if x<0: x+=1.0
    return x if x<=1-x else 1-x

def maxmin_float(S, refine=True):
    ks=set()
    for v in S:
        for j in range(v+1):
            ks.add(j/v)
    ks=sorted(ks); cand=list(ks)
    for i in range(len(ks)-1):
        a,b=ks[i],ks[i+1]
        cand.append(0.5*(a+b))
    best=0.0; bt=0.0
    for t in cand:
        m=1.0
        for v in S:
            d=fn_f(v*t)
            if d<m: m=d
            if m<best: break
        if m>best: best=m; bt=t
    if refine:
        # local refine around bt by bisection on the lower envelope
        lo,hi=bt-1e-4,bt+1e-4
        for _ in range(60):
            t1=lo+(hi-lo)/3; t2=hi-(hi-lo)/3
            m1=min(fn_f(v*t1) for v in S); m2=min(fn_f(v*t2) for v in S)
            if m1<m2: lo=t1
            else: hi=t2
        tm=(lo+hi)/2; mm=min(fn_f(v*tm) for v in S)
        if mm>best: best=mm; bt=tm
    return best,bt

def fn(x):
    f=x-(x.numerator//x.denominator); return f if f<=1-f else 1-f
def maxmin_exact(S):
    S=sorted(set(S)); cands=set()
    for v in S:
        for j in range(v+1): cands.add(Fr(j,v))
    n=len(S)
    for i in range(n):
        a=S[i]
        for jx in range(i+1,n):
            b=S[jx]
            for p in range(a+1):
                for q in range(b+1):
                    if a!=b:
                        t=Fr(p-q,a-b)
                        if 0<=t<=1: cands.add(t)
                    t2=Fr(p+q,a+b)
                    if 0<=t2<=1: cands.add(t2)
    best=Fr(-1); bt=None
    for t in cands:
        m=min(fn(v*t) for v in S)
        if m>best: best=m; bt=t
    return best,bt

THRESHf=1/14.0; THRESH=Fr(1,14)
random.seed(20260617)
best=(1.0,None)
N=200000
for it in range(N):
    rng=random.choice([14,20,30,50,80,120,200])
    S=random.sample(range(1,rng+1),13)
    mm,_=maxmin_float(S, refine=False)
    if mm<best[0]: best=(mm,tuple(sorted(S)))
print("random search (%d trials): lowest float max-min = %.7f (1/14=%.7f)"%(N,best[0],THRESHf))
print("  config:",best[1])
mm,t=maxmin_exact(list(best[1]))
print("  EXACT max-min:",mm,"=",float(mm),"  < 1/14?",mm<THRESH)

# local descent from several random starts to MINIMIZE max-min
print("\nlocal descent (minimize max-min):")
glob=(1.0,None)
for start in range(40):
    rng=random.choice([14,20,40,80,150])
    S=set(random.sample(range(1,rng+1),13))
    cur,_=maxmin_float(list(S),refine=False)
    improved=True
    while improved:
        improved=False
        Sl=list(S)
        for e in Sl:
            for w in range(1,201):
                if w in S: continue
                T=(S-{e})|{w}
                mm,_=maxmin_float(list(T),refine=False)
                if mm<cur-1e-12:
                    S=T; cur=mm; improved=True; break
            if improved: break
    if cur<glob[0]: glob=(cur,tuple(sorted(S)))
print("  best descent float max-min = %.7f"%glob[0],"config",glob[1])
mm,t=maxmin_exact(list(glob[1]))
print("  EXACT:",mm,"=",float(mm)," < 1/14?",mm<THRESH," at tau=",t)
print("\nCONCLUSION: smallest max-min found is", min(best[0],glob[0]),
      "vs 1/14 =",THRESHf,"-- counterexample found?", min(best[0],glob[0])<THRESHf-1e-9)
