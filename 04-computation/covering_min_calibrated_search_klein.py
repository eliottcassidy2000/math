#!/usr/bin/env python3
"""
covering_min_calibrated_search_klein.py  --  klein-2026-06-30-S60

NEW APPROACH to the covering-min OPEN EDGE (HYP-3764: is the construction n/Phi6 the covering-min for
n>=12?). Previous searches (random, hillclimb, drop<=2) DEMONSTRABLY MISSED known beaters (n=9,11), so
their "no beater at n>=12" was worthless. This search is CALIBRATED to the KNOWN beater structure:
  every verified beater = {1..n-1} minus a few small speeds + a few killers (<= ~3.5n), e.g.
    n=9  {1,3,4,5,7,11,18,32}  = drop {2,6,8} + killers {11,18,32}   (d=3, killer<=32<5n)
    n=10 {1,2,3,5,6,7,8,9,30}  = drop {4}     + killer {30}          (d=1)
So: enumerate drop-sets R (|R|<=DMAX from {2..n-1}, keep speed 1) + |R| killers in [n, B], keep primitive
COVERING sets (a multiple of every q in {2..n}), and test M < construction via fast early-exit M_lt.

CRUCIAL: VALIDATE on n=9,10,11 (known covering-min 4/33, 4/37, 3/31) BEFORE trusting n=12,13,14. If the
search recovers the known beaters, its verdict at n=12,13,14 is meaningful; if not, it is not (honest).
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

def Phi6(n): return n*n-n+1
def dist0(x,D):
    x%=D; return min(x,D-x)

def M_lt(S, thr):
    """True iff M(S) < thr (all gaps < thr), early-exit at first gap >= thr."""
    mx=max(S)
    for D in range(2, 2*mx+2):
        for a in range(1,D):
            m=D
            for s in S:
                d=dist0(s*a,D)
                if d<m: m=d
                if m==0: break
            if F(m,D)>=thr:
                return False
    return True

def M_exact(S):
    mx=max(S); best=F(0)
    for D in range(2,2*mx+2):
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

def search(n, DMAX, B):
    """min covering-set M over the calibrated family; return (best_M, best_set, #beaters, sample_beaters)."""
    thr=F(n,Phi6(n))
    small=list(range(2,n))                 # droppable (keep speed 1 always)
    killers=list(range(n, B+1))
    best=(thr, list(range(1,n-1))+[n*(n-1)])  # construction baseline
    beaters=[]
    for d in range(0, DMAX+1):
        for R in combinations(small, d):
            core=[1]+[x for x in small if x not in R]
            need = (n-1)-len(core)          # = d killers
            if need<0: continue
            if need==0:
                cands=[tuple()]
            else:
                cands=combinations(killers, need)
            for K in cands:
                S=sorted(set(core)|set(K))
                if len(S)!=n-1: continue
                if not prim(S) or not is_cov(S,n): continue
                if M_lt(S, best[0]):          # strictly beats current best
                    m=M_exact(S)
                    if m<best[0]: best=(m,S)
                    if m<thr: beaters.append((m,tuple(S)))
    beaters.sort()
    return best[0], best[1], len(beaters), beaters[:3]

if __name__=="__main__":
    print("VALIDATION (must recover known beaters): n=9->4/33, n=10->4/37, n=11->3/31")
    for n in [9,10,11]:
        thr=F(n,Phi6(n))
        bM,bS,nb,samp=search(n, DMAX=3, B=5*n)
        ok = bM<thr
        print(f"  n={n}: construction {thr}={float(thr):.4f}; SEARCH best {bM}={float(bM):.4f} "
              f"{'BEATS (validated)' if ok else 'FAILED to beat -- method invalid here'}; #beaters={nb}")
        if ok: print(f"        tightest: {list(bS)}  M={bM}")
    print()
    print("THE OPEN EDGE: n=12,13,14 (calibrated drop<=3 + killers<=5n):")
    for n in [12,13,14]:
        thr=F(n,Phi6(n))
        bM,bS,nb,samp=search(n, DMAX=3, B=5*n)
        verdict = "BEATS construction -- CONSTRUCTION IS NOT THE COVERING-MIN!" if bM<thr else "= construction (no beater in the calibrated family)"
        print(f"  n={n}: construction {thr}={float(thr):.5f}; SEARCH best {bM}={float(bM):.5f}  -> {verdict}")
        if bM<thr: print(f"        BEATER: {list(bS)}  M={bM}")
