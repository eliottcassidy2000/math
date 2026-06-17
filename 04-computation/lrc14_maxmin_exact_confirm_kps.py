#!/usr/bin/env python3
"""
lrc14_maxmin_exact_confirm_kps  (kind-pasteur 2026-06-17)
The random float screen produced FALSE positives (float maxmin UNDERestimates true
maxmin by missing crossing peaks). Here we EXACT-verify a large batch of the
lowest-float candidates to confirm: is there ANY 13-set with exact max-min < 1/14?
We collect the 400 lowest-float configs from a fresh random+structured sweep,
exact-verify all, and report the true minimum.
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
def maxmin_float(S):  # conservative screen (underestimates)
    ks=set()
    for v in S:
        for j in range(v+1): ks.add(j/v)
    ks=sorted(ks); cand=list(ks)
    for i in range(len(ks)-1): cand.append(0.5*(ks[i]+ks[i+1]))
    best=0.0
    for t in cand:
        m=1.0
        for v in S:
            d=fn_f(v*t)
            if d<m: m=d
        if m>best: best=m
    return best
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
THRESH=Fr(1,14)
random.seed(424242)
pool=[]
# random sweep, keep speeds modest so exact maxmin is tractable
for _ in range(120000):
    rng=random.choice([14,18,24,30,40,60])
    S=tuple(sorted(random.sample(range(1,rng+1),13)))
    pool.append((maxmin_float(S),S))
# structured: AP and all single-doublings & small perturbations
AP=list(range(1,14))
for e in AP:
    for w in range(14,61):
        if w in AP and w!=e: continue
        S=tuple(sorted([x for x in AP if x!=e]+[w]))
        pool.append((maxmin_float(S),S))
pool.sort(key=lambda x:x[0])
seen=set(); verified=[]
for mf,S in pool:
    if S in seen: continue
    seen.add(S)
    mm,t=maxmin_exact(S)
    verified.append((mm,S,t,mf))
    if len(verified)>=400: break
verified.sort(key=lambda x:x[0])
print("Exact-verified %d lowest-float candidates."%len(verified))
print("Smallest 12 by EXACT max-min:")
for mm,S,t,mf in verified[:12]:
    flag=" <<< COUNTEREXAMPLE" if mm<THRESH else (" (=1/14 tight)" if mm==THRESH else "")
    print("  exact=%s=%.7f (float screen said %.5f)%s  S=%s"%(mm,float(mm),mf,flag,S))
mn=min(v[0] for v in verified)
print("\nTRUE smallest exact max-min over all verified:",mn,"=",float(mn))
print("1/14 =",float(THRESH),"  Any below 1/14?",sum(1 for v in verified if v[0]<THRESH))
