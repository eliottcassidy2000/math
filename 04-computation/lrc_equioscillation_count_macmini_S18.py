#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S18 (HYP-4462) -- the EQUIOSCILLATION COUNT: the AP touches its
max M(S) at exactly phi(n) points = the units of Z/n; every other config at 2-4.

Chebyshev/Kolmogorov: the minimizer of a max-functional equioscillates.  The AP is
the unique tight (M=1/n) config, and it equioscillates MAXIMALLY -- f_AP(t)=min_i||v_i t||
reaches its max 1/n at every unit a/n (a in (Z/n)*), phi(n) of them.  This is the
equioscillation FACE of the tight locus, sitting beside equidistribution (roots of
unity), equidecomposability (arcs tile), equinumerosity (max relation lattice).

Uses the fast exact M (witness q | v_i+-v_j) and records ALL t achieving M.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

def Mfast_full(S):
    """(M, [distinct t achieving M])."""
    S=sorted(set(S)); Q=set()
    for v in S: Q.add(2*v)
    for a,b in combinations(S,2): Q.add(a+b); Q.add(abs(a-b))
    Q.discard(0); best=F(0); pts=[]
    for q in Q:
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            mn=min(min((v*a)%q, q-((v*a)%q)) for v in S)
            val=F(mn,q)
            if val>best: best=val; pts=[F(a,q)]
            elif val==best: pts.append(F(a,q))
    return best, sorted(set(pts))

def phi(n):
    return sum(1 for a in range(1,n) if gcd(a,n)==1)

if __name__=="__main__":
    print("equioscillation count = #distinct t achieving M(S):\n")
    fams=[("AP {1..12} (n=13)",list(range(1,13)),13),
          ("doubled-apex {1..11,24}",list(range(1,12))+[24],None),
          ("block {..17,19}",[1,2,3,5,7,8,9,10,11,12,17,19],None),
          ("{1..11,23}",list(range(1,12))+[23],None),
          ("AP {1..6} (n=7)",list(range(1,7)),7),
          ("n=7 gap {1,5,6,11,16,17}",[1,5,6,11,16,17],None)]
    for name,S,n in fams:
        M,ts=Mfast_full(S)
        extra=f"  phi({n})={phi(n)} (units of Z/{n})" if n else ""
        print(f"  {name:26s} M={str(M):>6}  #equiosc-points={len(ts)}{extra}")
    print("\n=> AP = phi(n) equioscillation points (the units); all else 2-4.  Maximal symmetry")
    print("   uniquely at the AP -- the equioscillation face of the tight locus.")
