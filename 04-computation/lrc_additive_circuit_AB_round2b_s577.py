#!/usr/bin/env python3
"""Round 2b (light): the bridge -- near-tight configs are 3-term-rich AND all small
pairs shielded (=> global-clock worry-set = ker Phi = C'). opus-2026-06-03-S577."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random
def dist(x): x%=1; return min(x,1-x)
def Mfloat(V,k):
    # fast float M via crossing times t=m/(vi+-vj); good enough to flag near-tight
    cs=set()
    for i in range(len(V)):
        for j in range(i+1,len(V)):
            for D in (V[i]+V[j],abs(V[i]-V[j])):
                if D: 
                    for m in range(1,D): cs.add(m/D)
    best=0.0
    for t in cs:
        mn=min(min((v*t)%1,1-((v*t)%1)) for v in V)
        if mn>best: best=mn
    return best
def three_terms(V):
    s=set(V); return [(a,b,a+b) for a,b in combinations(sorted(V),2) if a+b in s]
def small_pairs(V,k):
    return [(a,b) for a,b in combinations(sorted(V),2) if (a+b)//gcd(a,b)<=k+1]
def all_small_shielded(V,k):
    sp=small_pairs(V,k)
    return bool(sp) and all(any(c%(a+b)==0 for c in V if c not in (a,b)) for (a,b) in sp)
def prim(V):
    g=0
    for v in V: g=gcd(g,v)
    return tuple(sorted(v//g for v in V))
def main():
    rng=random.Random(99)
    print("Bridge: near-tight (M-delta<0.02) configs -- 3-term-rich AND all-small-shielded?")
    for k in [6,8,10]:
        delta=1/(k+1); B=2*k+6; near=sh=cf=0; m3=0
        for _ in range(20000):
            V=prim(tuple(sorted(rng.sample(range(1,B+1),k))))
            if len(V)!=k: continue
            M=Mfloat(V,k)
            if M-delta<0.02:
                near+=1; t3=len(three_terms(V)); m3+=t3
                if t3==0: cf+=1
                if all_small_shielded(V,k): sh+=1
        if near:
            print(f"  k={k:2d}: near-tight={near}; mean#3term={m3/near:.1f}; circuit-free among them={cf} "
                  f"({100*cf/near:.1f}%); all-small-shielded={sh} ({100*sh/near:.0f}%)")
if __name__=='__main__': main()
