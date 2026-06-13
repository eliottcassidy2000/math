#!/usr/bin/env python3
"""Round 2 (Lemma B structure): hardness vs 3-term richness, and the bridge to the
all-shielded worry-set. opus-2026-06-03-S577."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random
def dist(x): x%=1; return min(x,1-x)
def Mexact(V):
    cands=set()
    for i in range(len(V)):
        for j in range(len(V)):
            if i==j: continue
            for D in (V[i]+V[j], abs(V[i]-V[j])):
                if D==0: continue
                for m in range(1,D): cands.add(F(m,D))
        cands.add(F(1,2*V[i]))
    best=F(0)
    for t in cands:
        mn=min(dist(v*t) for v in V)
        if mn>best: best=mn
    return best
def three_terms(V):
    s=set(V); return [(a,b,a+b) for a,b in combinations(sorted(V),2) if a+b in s]
def small_pairs(V,k):  # reduced sum <= k+1
    return [(a,b) for a,b in combinations(sorted(V),2) if (a+b)//gcd(a,b)<=k+1]
def shielded(V,a,b,k):  # pair (a,b) blocked by a SHIELD: some c with (a+b)|c (incl c=a+b)
    D=a+b
    return any(c%D==0 for c in V if c not in (a,b))
def all_small_shielded(V,k):
    sp=small_pairs(V,k)
    return bool(sp) and all(shielded(V,a,b,k) for (a,b) in sp)
def prim(V):
    g=0
    for v in V: g=gcd(g,v)
    return tuple(sorted(v//g for v in V))
def main():
    rng=random.Random(13)
    print("(a) M-delta binned by #3-term relations (hard configs are 3-term-rich?)")
    for k in [6,8,10]:
        delta=F(1,k+1); B=3*k+2; bins={}
        for _ in range(8000):
            V=prim(tuple(sorted(rng.sample(range(1,B+1),k))))
            if len(V)!=k: continue
            t3=len(three_terms(V)); M=Mexact(V)
            key=min(t3,5)
            bins.setdefault(key,[]).append(float(M-delta))
        print(f"  k={k}: "+"; ".join(f"#3t={key}{'+' if key==5 else ''}: minMargin={min(v):+.4f}(n={len(v)})" for key,v in sorted(bins.items())))
    print()
    print("(b) near-tight configs (M-delta < 0.02): are all small pairs shielded? 3-term count?")
    for k in [6,8,10]:
        delta=F(1,k+1); B=3*k+2; near=0; shielded_all=0; mean3=0
        for _ in range(12000):
            V=prim(tuple(sorted(rng.sample(range(1,B+1),k))))
            if len(V)!=k: continue
            M=Mexact(V)
            if float(M-delta)<0.02:
                near+=1; mean3+=len(three_terms(V))
                if all_small_shielded(V,k): shielded_all+=1
        if near: print(f"  k={k}: near-tight configs={near}; all-small-pairs-shielded={shielded_all}/{near} "
                       f"({100*shielded_all/near:.0f}%); mean #3-term={mean3/near:.1f}")
if __name__=='__main__': main()
