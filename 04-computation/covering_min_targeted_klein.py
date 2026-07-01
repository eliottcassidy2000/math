#!/usr/bin/env python3
"""
covering_min_targeted_klein.py  --  klein-2026-06-30-S53

Targeted structured beater search at n=12,13,14 (the open edge).  The found beaters are NEARLY DENSE:
n=10 covering-min = {1,2,3,5,6,7,8,9,30} = {1..9}\{4} + killer 30 (M=4/37, rung 4).  Structure = drop a
REDUNDANT small speed (its resonance still covered) + add a tuned killer, landing at a LOW rung.
So: enumerate {1,..,n-1} minus a small dropped subset, + killers (<=4n), keep primitive covering sets
of size n-1, compute the FULL exact M, and report any that beats the construction n/Phi6.
Motivation: random/hillclimb MISS beaters (missed 4/33 at n=9, 3/31 at n=11) -- a structured sweep is
the honest test of whether the construction is really the covering-min for n>=12.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

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
def rung_of(Mv,n):
    for r in range(1,4*n):
        if Mv==F(r,r*(n-1)+1): return r
    return None

def search(n, drop_max=2, Kbound=None, Kadd_max=2):
    Kbound = Kbound or 5*n
    thr=F(n,Phi6(n))
    base=list(range(1,n))            # {1,..,n-1}
    small=list(range(2,n))           # droppable (keep 1)
    killers=list(range(n, Kbound+1)) # candidate large speeds
    best=(F(n,Phi6(n)), list(range(1,n-1))+[n*(n-1)])
    beaters=[]
    for dsz in range(0, drop_max+1):
        for drp in combinations(small, dsz):
            kept=[x for x in base if x not in drp]
            need = (n-1) - len(kept)          # killers to add
            if need < 0: continue
            if need==0:
                cand=[kept]
            elif need==1:
                cand=[kept+[k] for k in killers]
            elif need==2:
                cand=[kept+[k1,k2] for i,k1 in enumerate(killers) for k2 in killers[i+1:]]
            else:
                continue
            if need>Kadd_max: continue
            for S in cand:
                S=sorted(set(S))
                if len(S)!=n-1: continue
                if not prim(S) or not is_cov(S,n): continue
                mv=M(S)
                if mv<best[0]:
                    best=(mv,S)
                if mv<thr:
                    beaters.append((mv,S))
    return thr,best,beaters

if __name__=="__main__":
    for n in [12,13,14]:
        # drop up to 2, add up to 2 killers (need == drop size + (n-1 - (n-1)) ... kept=n-1-dsz, need=dsz)
        thr,best,beaters=search(n, drop_max=2, Kbound=5*n, Kadd_max=2)
        bM,bS=best; r=rung_of(bM,n)
        print(f"n={n}: construction {thr}={float(thr):.5f}; BEST found {bM}={float(bM):.5f} rung={r} "
              f"{'BEATS!' if bM<thr else '= construction'}")
        if beaters:
            beaters.sort()
            print(f"   {len(beaters)} beaters; tightest: {beaters[0][1]} M={beaters[0][0]}={float(beaters[0][0]):.5f} rung={rung_of(beaters[0][0],n)}")
        else:
            print(f"   0 beaters in the drop<=2 + <=2-killer(<=5n) structured family")
