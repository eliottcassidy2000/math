#!/usr/bin/env python3
"""
covering_min_hillclimb_klein.py  --  klein-2026-06-30-S53

Does a covering set BEAT the construction n/Phi6(n) for n=12,13,14?  (The open edge.)
Random sampling misses the SPREAD-structured beaters (verified real at n=7,8,9: e.g.
{1,3,4,5,7,11,18,32} at n=9 has M=4/33 < 9/73 and does NOT contain the core -- so the construction
is the loosest rung, not the covering-min).  HILL-CLIMBING (local swap search, minimizing M over
covering sets) is far better at finding them.  Start from the construction + random covering sets +
scaled beaters; swap one speed at a time to lower M while staying a primitive covering set.

Reports, per n: construction value, best M found, and whether it BEATS the construction (rung < n).
"""
from fractions import Fraction as F
from math import gcd
import random

def Phi6(n): return n*n-n+1
def dist0(x,D):
    x%=D; return min(x,D-x)
def M(S):
    Dmax=2*max(S)+2; best=F(0)
    for D in range(2,Dmax+1):
        for a in range(1,D):
            m=D
            for s in S:
                d=dist0(s*a,D); m=d if d<m else m
                if m==0: break
            if F(m,D)>best: best=F(m,D)
    return best
def is_cov(S,n): return all(any(s%q==0 for s in S) for q in range(2,n+1))
def prim(S):
    g=0
    for s in S: g=gcd(g,s)
    return g==1
def ok(S,n): return len(set(S))==n-1 and prim(S) and is_cov(S,n)

def rung_of(Mv, n):
    # M = r/(r(n-1)+1) => r = M/(1 - M(n-1)); report nearest simple rung if it fits
    for r in range(1, 3*n):
        if Mv == F(r, r*(n-1)+1): return r
    return None

def random_cov(n, Bmax, rng):
    S=set()
    for q in sorted(range(2,n+1), key=lambda _:rng.random()):
        if any(s%q==0 for s in S): continue
        if len(S)>=n-1: return None
        S.add(q*rng.randint(1,max(1,Bmax//q)))
    t=0
    while len(S)<n-1 and t<300:
        S.add(rng.randint(1,Bmax)); t+=1
    S=sorted(S)
    return S if ok(S,n) else None

def hillclimb(S0, n, Bmax, rng, iters=400):
    S=list(S0); bestM=M(S)
    for _ in range(iters):
        i=rng.randrange(len(S))
        cand=rng.randint(1,Bmax)
        T=sorted(set(S[:i]+S[i+1:]+[cand]))
        if len(T)!=n-1 or not ok(T,n): continue
        mt=M(T)
        if mt<bestM:
            S=T; bestM=mt
    return S,bestM

if __name__=="__main__":
    print(f"{'n':>3} {'construction':>14} {'best found (hillclimb)':>24} {'rung':>5}  verdict")
    print("-"*74)
    for n in [12,13,14]:
        thr=F(n,Phi6(n)); rng=random.Random(777+n)
        constr=list(range(1,n-1))+[n*(n-1)]
        best=(M(constr), constr)
        Bmax=4*n*(n-1)
        restarts = 60 if n<=13 else 40
        # seed with construction, random covering sets
        seeds=[constr]
        for _ in range(restarts):
            s=random_cov(n,Bmax,rng)
            if s: seeds.append(s)
        for s in seeds:
            Sc,Mc=hillclimb(s,n,Bmax,rng, iters=250)
            if Mc<best[0]: best=(Mc,Sc)
        bM,bS=best
        r=rung_of(bM,n)
        verdict = ("= construction (rung n)" if bM==thr else
                   f"BEATS construction!" if bM<thr else "worse")
        print(f"{n:>3} {str(thr):>7}={float(thr):.5f} {str(bM):>10}={float(bM):.5f} {str(r):>5}  {verdict}")
        if bM<thr: print(f"      beater set: {bS}")
    print()
    print("hillclimb finds spread beaters random sampling misses; if best==construction the construction")
    print("survives local search (evidence, not proof, that it is the covering-min at that n).")
