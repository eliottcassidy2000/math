#!/usr/bin/env python3
"""exp turns +(line) into x(circle); runners are helices. LRC hardness grades INVERSELY
with the hyperoperation level of speed growth: additive(AP)=hardest, x2/geometric=easy,
exponential/tetration=trivial. The clock witnesses = roots of unity; worry-set =
PRIMITIVE roots (cyclotomic). opus-2026-06-03-S588."""
from fractions import Fraction as F
from math import gcd
def dist(x): x%=1; return min(x,1-x)
def Mval(V):
    V=[v for v in V if v]; 
    if not V: return F(1,2)
    from itertools import combinations
    cands=set()
    for a,b in combinations(V,2):
        for D in (a+b,abs(a-b)):
            if D:
                for m in range(1,D): cands.add(F(m,D))
    for v in V: cands.add(F(1,2*v))
    return max(min(dist(v*t) for v in V) for t in cands)
def main():
    print("Hardness vs HYPEROPERATION level of speed growth (k speeds, delta=1/(k+1)):")
    for k in [5,6,7]:
        delta=F(1,k+1)
        AP=tuple(range(1,k+1))                         # level 1: ADDITION (dense)
        geo=tuple(2**i for i in range(k))              # level 2->3: x2 / geometric (lacunary)
        fast=tuple([1]+[3**i for i in range(1,k)])     # faster (more lacunary)
        for name,V in [("additive AP",AP),("geometric 2^i",geo),("lacunary 3^i",fast)]:
            print(f"   k={k} {name:14s} {V}: M-delta = {float(Mval(V)-delta):+.4f}")
        print()
    print("exp realizes the clock as ROOTS OF UNITY; AP witnesses = PRIMITIVE roots (cyclotomic):")
    for n in [7,12,14]:
        AP=tuple(range(1,n))
        W=[j for j in range(1,n) if all(dist(v*F(j,n))>=F(1,n) for v in AP)]
        prim=[j for j in range(1,n) if gcd(j,n)==1]
        print(f"   n={n}: AP clock-witnesses {W} = primitive n-th roots e^(2pi i j/n), j in (Z/n)* {prim}: {W==prim}")
if __name__=='__main__': main()
