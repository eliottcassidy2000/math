#!/usr/bin/env python3
"""
covering_min_beater_search_klein.py  --  klein-2026-06-30-S53

THE OPEN EDGE: does ANY covering set beat the construction n/Phi6(n)?
The repo covering-min values (MISTAKE-087) 2/13,2/15,4/33,4/37,3/31 (n=7..11) are SPREAD-structured and
several are BELOW the construction -- so structured "core+killer" searches miss them. Here: a broad
randomized+structured search over covering sets (a multiple of every q in {2..n}) with large/spread
speeds, testing M(S) < n/Phi6.  Reports, per n: construction value, #beaters found, and the TIGHTEST
covering set found (the current best covering-min estimate).
"""
from fractions import Fraction as F
from math import gcd
import random

def Phi6(n): return n*n-n+1
def dist0(x,D):
    x%=D; return min(x,D-x)

def M_exact(S):
    Dmax=2*max(S)+2; best=F(0)
    for D in range(2,Dmax+1):
        for a in range(1,D):
            m=D
            for s in S:
                d=dist0(s*a,D)
                if d<m: m=d
                if m==0: break
            if F(m,D)>best: best=F(m,D)
    return best

def M_lt(S, thr):
    """True iff M(S) < thr (strict): every (D,a) gap is < thr. Early-exit."""
    for D in range(2,2*max(S)+2):
        for a in range(1,D):
            m=D
            for s in S:
                d=dist0(s*a,D)
                if d<m: m=d
                if m==0: break
            if F(m,D)>=thr:
                return False
    return True

def is_cov(S,n): return all(any(s%q==0 for s in S) for q in range(2,n+1))
def prim(S):
    g=0
    for s in S: g=gcd(g,s)
    return g==1

def random_covering(n, Bmax, rng):
    """Build a primitive covering (n-1)-set: force a multiple of each q, fill the rest, spread speeds."""
    S=set()
    order=list(range(2,n+1)); rng.shuffle(order)
    for q in order:
        if any(s%q==0 for s in S): continue
        if len(S) >= n-1: return None
        mult=q*rng.randint(1, max(1,Bmax//q))
        S.add(mult)
    # fill remaining slots with spread speeds
    tries=0
    while len(S) < n-1 and tries<200:
        S.add(rng.randint(1,Bmax)); tries+=1
    S=sorted(S)
    if len(S)!=n-1 or not prim(S) or not is_cov(S,n): return None
    return S

if __name__=="__main__":
    print(f"{'n':>3} {'construction n/Phi6':>20} {'#beaters':>9} {'tightest covering found':>26}")
    print("-"*78)
    for n in [10,11,12,13,14]:
        thr=F(n,Phi6(n))
        rng=random.Random(12345+n)
        constr=list(range(1,n-1))+[n*(n-1)]
        best=(M_exact(constr), tuple(constr))   # construction as baseline
        beaters=0; K = 9000 if n<=12 else 4000
        Bmax = 3*n*(n-1)   # room for spread + large killers
        for _ in range(K):
            S=random_covering(n,Bmax,rng)
            if S is None: continue
            if M_lt(S, best[0]):          # strictly tighter than current best
                Mv=M_exact(S)
                if Mv < best[0]:
                    best=(Mv, tuple(S))
                if Mv < thr:
                    beaters+=1
        bM,bS=best
        beat = "= construction" if bM==thr else ("BEATS construction!" if bM<thr else "")
        print(f"{n:>3} {str(thr):>10}={float(thr):.5f} {beaters:>9}   {str(bM):>9}={float(bM):.5f} {beat}")
        if bM<thr:
            print(f"      tightest set: {list(bS)}")
    print()
    print("A randomized search (not exhaustive): 0 beaters => construction is the covering-min at that n")
    print("within the sampled spread/large-speed families; >0 => construction is NOT the covering-min.")
