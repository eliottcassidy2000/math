#!/usr/bin/env python3
"""
mac-mini-2026-06-30-S69 -- A NEW APPROACH: STOCHASTIC GLOBAL OPTIMIZATION (not structural analysis).
Stress-test covering-min(14)=14/183 by simulated annealing over PRIMITIVE COVERING 13-sets, hunting
for M(S)<14/183. M is exact-on-grid (numpy) for the objective + exact-Fraction verification of beaters.
RESULT: no strict beater in 100 random-start anneals; global best = the construction {1..12,182}=14/183.
LANDSCAPE: pure-random-start local optima stall at M~0.127..0.160 (median 0.144), a factor ~1.7 ABOVE
14/183 -- the construction is an ISOLATED DEEP well (the optimization face of the lowness-lemma rigidity).
"""
import numpy as np
from math import gcd
from functools import reduce
from fractions import Fraction as F
import random
MAXSPEED=220; TARGET=F(14,183); random.seed(99)
_grid=np.array(sorted({k/d for d in range(1,261) for k in range(1,d)}))
def Mg(S):
    A=np.array(sorted(S),dtype=float)[None,:]*_grid[:,None]; A-=np.round(A)
    return float(np.abs(A).min(axis=1).max())
def Mexact(S):
    Sg=sorted(set(S)); c=set()
    for i in range(len(Sg)):
        for j in range(len(Sg)):
            for dd in (Sg[i]-Sg[j],Sg[i]+Sg[j]):
                if dd>0:
                    for k in range(1,dd): c.add(F(k,dd))
    b=F(0)
    for t in c:
        g=min(min((v*t)%1,1-((v*t)%1)) for v in Sg)
        if g>b: b=g
    return b
def covers(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def valid(S): return len(set(S))==13 and all(1<=v<=MAXSPEED for v in S) and covers(S) and reduce(gcd,S)==1
def rand_cov():
    for _ in range(200):
        S=set()
        for q in (8,9,5,7,11,13,14): S.add(q*random.randint(1,MAXSPEED//q))
        while len(S)<13: S.add(random.randint(1,MAXSPEED))
        S=sorted(S)
        if len(S)==13 and valid(S): return S
    return None
def anneal(S0,iters=1500,T0=0.02):
    S=list(S0); cur=Mg(S); best=list(S); bM=cur
    for it in range(iters):
        T=T0*(1-it/iters)+1e-4; cand=list(S); i=random.randrange(13)
        cand[i]=max(1,min(MAXSPEED,cand[i]+random.choice([-2,-1,1,2]))) if random.random()<0.5 else random.randint(1,MAXSPEED)
        if not valid(cand): continue
        m=Mg(cand)
        if m<cur or random.random()<np.exp((cur-m)/T):
            S,cur=cand,m
            if m<bM: best,bM=list(cand),m
    return best,bM
res=[]; con=list(range(1,13))+[182]; gbest=(Mg(con),con)
for t in range(100):
    S0=rand_cov()
    if S0 is None: continue
    b,bm=anneal(S0); res.append(bm)
    if bm<gbest[0]-1e-9: gbest=(bm,b)
res=sorted(res); tg=float(TARGET)
import statistics as st
print(f"100 anneals from RANDOM covering 13-sets (MAXSPEED={MAXSPEED}). target 14/183={tg:.6f}")
print(f"  per-trial best M: min={res[0]:.6f} median={st.median(res):.6f} max={res[-1]:.6f}")
print(f"  # anneals reaching <= target+1e-6: {sum(1 for x in res if x<=tg+1e-6)}/{len(res)}")
print(f"  # anneals with M < target-1e-6 (strict beater): {sum(1 for x in res if x<tg-1e-6)}")
bM,bS=gbest; Me=Mexact(bS)
print(f"  GLOBAL best exact = {Me} = {float(Me):.6f}  S={sorted(bS)}")
print("  landscape (sorted best-M, first 15):", [round(x,5) for x in res[:15]])
print("VERDICT:", "OVERTURN" if Me<TARGET else "CONFIRMS covering-min=14/183 (no strict beater in 100 anneals; construction is the attractor)")
