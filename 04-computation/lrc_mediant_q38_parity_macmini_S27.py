#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S27 (HYP-4572) -- CROSS-N leverage on opus-S117's crux:
the mediant 3/(3N+2) always fits the window (=> (G)=non-achievability), and the
N=12 crux 3/38 = q=38=2*19 descends by PARITY to a mod-19 clearance-2 feasibility.
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

if __name__=="__main__":
    print("KNOWN gap members, corrected order k = q - (#speeds)*s:")
    for name,S in [("6sp {1,5,6,11,16,17}",[1,5,6,11,16,17]),("7sp {1,3,4,5,7,13,18}",[1,3,4,5,7,13,18])]:
        M=Mfast(S); nsp=len(S); k=M.denominator-nsp*M.numerator
        lo,hi=F(1,nsp+1),F(2,2*nsp+1)
        print(f"  {name}: M={M} in gap({lo},{hi})={lo<M<hi}; s={M.numerator}, order k={k}")
    print("\nMEDIANT 3/(3N+2) across N -- scale-invariant ~0.65 of window (value always fits):")
    for N in [6,7,8,10,12]:
        med=F(3,3*N+2); lo,hi=F(1,N+1),F(2,2*N+1)
        print(f"  N={N:>2}: {med}={float(med):.5f}; rise {float(med-lo):.5f} = {float((med-lo)/(hi-lo)):.2f} of window")
    print("\nq=38=2*19 mediant 3/38 hole {0,1,2,36,37} by CRT (mod2 x mod19):")
    for r in [0,1,2,36,37]:
        print(f"  {r} = (mod2 {r%2}, mod19 {r%19})")
    print("  => EVEN speeds 2w HALVE to clearance-2 mod 19 (wa avoids {0,+-1} mod 19); ODD avoid {+-1} mod 19.")
    print("  The mediant reduces to a mod-19 clearance-2 COVERING feasibility (opus O-mediant + my E_p/O_p seed).")
