#!/usr/bin/env python3
"""Round 3: the convergence. Lemma B (3-term-rich) hard configs are lonely via the
straddle pinch t=j/(k+1) (the delta-clock), blocked only by a multiple of (k+1) = C'.
Verify near-tight 3-term-rich configs: (a) no multiple of k+1, (b) delta-clock witness.
opus-2026-06-03-S577."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random
def dist(x): x%=1; return min(x,1-x)
def Mfloat(V):
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
    s=set(V); return sum(1 for a,b in combinations(sorted(V),2) if a+b in s)
def delta_clock_witness(V,k):  # some j in 1..k with all ||v_i j/(k+1)|| >= 1/(k+1)
    n=k+1
    return any(all(dist(v*F(j,n))>=F(1,n) for v in V) for j in range(1,n))
def prim(V):
    g=0
    for v in V: g=gcd(g,v)
    return tuple(sorted(v//g for v in V))
def main():
    rng=random.Random(202)
    print("Convergence: near-tight 3-term-rich configs -> no mult of (k+1) AND delta-clock witness?")
    for k in [6,8,10,12]:
        n=k+1; delta=1/n; B=2*k+6
        near=0; nomult=0; clockwit=0; hasmult=0
        for _ in range(25000):
            V=prim(tuple(sorted(rng.sample(range(1,B+1),k))))
            if len(V)!=k: continue
            if Mfloat(V)-delta<0.02:
                near+=1
                hasm=any(v%n==0 for v in V)
                if not hasm: nomult+=1
                else: hasmult+=1
                if delta_clock_witness(V,k): clockwit+=1
        if near:
            print(f"  k={k:2d} (n={n}): near-tight={near}; NO multiple of {n}={nomult} ({100*nomult/near:.0f}%); "
                  f"delta-clock j/{n} witness={clockwit} ({100*clockwit/near:.0f}%); had-multiple={hasmult}")
if __name__=='__main__': main()
