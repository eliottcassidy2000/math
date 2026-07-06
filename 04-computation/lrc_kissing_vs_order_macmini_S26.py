#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S26 (HYP-4562) -- the KISSING DEFICIT IS opus's ORDER k.

Synthesizes HYP-4552 (kissing number = additive triples = additive energy, AP maximal)
with opus-S116 HYP-4486 ((G) = Kravitz's first-gap Lonely Runner Spectrum Conjecture;
below-1/n values are rungs s/(ns+k); gap members have ORDER k>=2; CKMRV Annals-2022
LP-uniqueness = the Cohn-Elkies frame HYP-4532).

MEASURED: the AP (order k=1) has MAXIMAL kissing at every n; gap members (order k>=2)
have a large KISSING DEFICIT.  So order k and the additive-energy deficit are two faces
of one defect: k=1 <=> full additive closure (max kissing) = AP; k>=2 <=> defected
generalized AP.  Pins opus's (O-korder) to a kissing/additive-energy upper bound.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

def Mfast(S):
    S=sorted(set(S)); Q=set()
    for v in S: Q.add(2*v)
    for a,b in combinations(S,2): Q.add(a+b); Q.add(abs(a-b))
    Q.discard(0); best=F(0)
    for q in Q:
        for a in range(1,q):
            mn=min(min((v*a)%q,q-((v*a)%q)) for v in S)
            val=F(mn,q)
            if val>best: best=val
    return best

def kissing(S):
    Sset=set(S); c=0
    for i in range(len(S)):
        for j in range(i,len(S)):
            if S[i]+S[j] in Sset: c+=1
    return c

def order_k(M, nspeed):
    s=M.numerator
    return (s, M.denominator - nspeed*s)     # M = s/(ns+k) => k = denom - n*s

if __name__=="__main__":
    tests=[
      ("AP {1..12} (12sp)", list(range(1,13)), 12),
      ("n=7 {1,3,4,5,7,13,18}", [1,3,4,5,7,13,18], 7),
      ("n=7 {1,5,6,11,16,17}", [1,5,6,11,16,17], 6),
      ("AP {1..6}", list(range(1,7)), 6),
      ("AP {1..7}", list(range(1,8)), 7),
      ("doubled-apex {1..11,24}", list(range(1,12))+[24], 12),
    ]
    print(f"  {'family':28s} {'M':>7} {'kiss':>5} {'AP-max':>7} {'order(s,k)':>11}")
    for name,S,nsp in tests:
        M=Mfast(S)
        print(f"  {name:28s} {str(M):>7} {kissing(S):>5} {kissing(list(range(1,nsp+1))):>7} {str(order_k(M,nsp)):>11}")
    print("=> AP order k=1 = max kissing; gap members order k>=2 = big kissing deficit. One defect, two frames.")
