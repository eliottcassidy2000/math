#!/usr/bin/env python3
"""
covering_min_probe_klein.py  --  klein-2026-06-30-S53

Two-stage beater probe (fast): does a covering set beat the construction n/Phi6(n) at n=12,13,14?
Low-rung beaters bind at a SMALL modulus D0 = a(n-1)+1 (rung a small).  Stage 1: hillclimb minimizing
a CAPPED-D surrogate M_cap (D<=Dcap) -- fast, catches low-rung structure.  Stage 2: VERIFY each
promising candidate with the FULL exact M (uncapped) -- a set that is tight up to Dcap but has a hole
at larger D is rejected here.  So a reported beater is a genuine covering set with full M < n/Phi6.

Also re-confirms the known beaters n=7,8,9 as a sanity check that the method finds spread beaters.
"""
from fractions import Fraction as F
from math import gcd
import random

def Phi6(n): return n*n-n+1
def dist0(x,D):
    x%=D; return min(x,D-x)

def M_cap(S, Dcap):
    best=F(0)
    for D in range(2,Dcap+1):
        for a in range(1,D):
            m=D
            for s in S:
                d=dist0(s*a,D); m=d if d<m else m
                if m==0: break
            if F(m,D)>best: best=F(m,D)
    return best

def M_full(S):
    return M_cap(S, 2*max(S)+1)

def is_cov(S,n): return all(any(s%q==0 for s in S) for q in range(2,n+1))
def prim(S):
    g=0
    for s in S: g=gcd(g,s)
    return g==1
def ok(S,n): return len(set(S))==n-1 and prim(S) and is_cov(S,n)

def rung_of(Mv,n):
    for r in range(1,3*n):
        if Mv==F(r,r*(n-1)+1): return r
    return None

def rand_cov(n,B,rng):
    S=set()
    for q in sorted(range(2,n+1),key=lambda _:rng.random()):
        if any(s%q==0 for s in S): continue
        if len(S)>=n-1: return None
        S.add(q*rng.randint(1,max(1,B//q)))
    t=0
    while len(S)<n-1 and t<300: S.add(rng.randint(1,B)); t+=1
    S=sorted(S); return S if ok(S,n) else None

def climb(S,n,B,rng,Dcap,iters):
    best=M_cap(S,Dcap)
    for _ in range(iters):
        i=rng.randrange(len(S))
        T=sorted(set(S[:i]+S[i+1:]+[rng.randint(1,B)]))
        if len(T)!=n-1 or not ok(T,n): continue
        mt=M_cap(T,Dcap)
        if mt<best: S=T; best=mt
    return S,best

def probe(n, rng, restarts, iters):
    thr=F(n,Phi6(n)); B=3*n; Dcap=4*n  # low speeds/moduli where low-rung beaters live
    constr=list(range(1,n-1))+[n*(n-1)]
    cands=[]
    seeds=[rand_cov(n,B,rng) for _ in range(restarts)]
    seeds=[s for s in seeds if s]
    for s in seeds:
        Sc,_=climb(s,n,B,rng,Dcap,iters)
        cands.append(Sc)
    # stage 2: full exact M on the distinct candidates
    best=(M_full(constr),tuple(constr))
    seen=set()
    for S in cands:
        key=tuple(S)
        if key in seen: continue
        seen.add(key)
        mf=M_full(S)
        if mf<best[0]: best=(mf,key)
    return thr,best

if __name__=="__main__":
    # sanity: known beaters n=7,8,9
    print("SANITY (known beaters, full M):")
    for n,S in [(7,[1,2,5,6,7,8]),(8,[1,4,5,6,7,11,16]),(9,[1,3,4,5,7,11,18,32])]:
        print(f"  n={n} M={M_full(S)} rung={rung_of(M_full(S),n)} < constr {F(n,Phi6(n))}? {M_full(S)<F(n,Phi6(n))}")
    print()
    print(f"{'n':>3} {'construction':>14} {'best found (2-stage probe)':>27} {'rung':>5}  verdict")
    print("-"*76)
    for n in [12,13,14]:
        rng=random.Random(4200+n)
        thr,(bM,bS)=probe(n,rng, restarts=(120 if n<=13 else 90), iters=180)
        r=rung_of(bM,n)
        v = "= construction" if bM==thr else ("BEATS!" if bM<thr else "worse")
        print(f"{n:>3} {str(thr):>7}={float(thr):.5f} {str(bM):>11}={float(bM):.5f} {str(r):>5}  {v}")
        if bM<thr: print(f"      beater: {list(bS)}")
