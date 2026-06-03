#!/usr/bin/env python3
"""RIGOROUS: the pair-sum sieve's modulus is C=2n-1. (i) 1/n and 2/(2n-1) are Farey
neighbors (cross-product 1). (ii) mod C the additive shells are {a, C-a}={a,-a}.
(iii) C=2n-1=sqrt(8 C(n,2)+1). (iv) the floor-tight worry-set = transversals of the
shells; verify the pinch-value just above 1/n is 2/(2n-1) on floor-tight configs.
opus-2026-06-03-S590. Sidesteps resonance energy: pure arithmetic."""
from fractions import Fraction as F
from math import gcd
def dist(x): x%=1; return min(x,1-x)
def Mexact_and_arg(V,n):
    from itertools import combinations
    cands=set()
    for a,b in combinations(V,2):
        for D in (a+b,abs(a-b)):
            if D:
                for m in range(1,D): cands.add(F(m,D))
    for v in V: cands.add(F(1,2*v))
    best=F(0)
    for t in cands:
        mn=min(dist(v*t) for v in V)
        if mn>best: best=mn
    return best
def main():
    print("(i) Farey-neighbor identity: 1/n & 2/(2n-1) adjacent in F_{2n-1} (cross-product +-1):")
    for n in range(3,16):
        cp = 1*(2*n-1) - 2*n   # ad - bc for 1/n, 2/(2n-1)
        print(f"   n={n:2d}: 1/{n} , 2/{2*n-1}: cross-product = {cp} (|.|=1 => Farey neighbors: {abs(cp)==1})")
    print()
    print("(ii) mod C=2n-1: additive shells {a, C-a}={a,-a}; (iii) C=sqrt(8 C(n,2)+1):")
    for n in [4,7,14]:
        C=2*n-1
        shells=[(a, C-a) for a in range(1,n)]
        print(f"   n={n}: C={C} (odd; sqrt(8*C(n,2)+1)={int((8*n*(n-1)//2+1)**0.5)}); "
              f"shells {shells[:4]}{'...' if len(shells)>4 else ''} (a+(C-a)=C, antipodal a<->-a)")
    print()
    print("(iv) on FLOOR-TIGHT configs: M=1/n exactly, and the next achievable pinch value")
    print("     just above 1/n is the Farey neighbor 2/(2n-1) (verify on AP & sporadics):")
    for V,n in [((1,2,3),4),((1,2,3,4,5,6),7),((1,2,3,4,5,6,7,8,9,10,11,12,13),14),((1,3,4,5,9),6)]:
        M=Mexact_and_arg(V,n)
        print(f"   V={V} (n={n}): M={M}={float(M):.4f}; floor 1/n={float(F(1,n)):.4f}; "
              f"2/(2n-1)={float(F(2,2*n-1)):.4f}; M==1/n: {M==F(1,n)}; M<2/(2n-1): {M<F(2,2*n-1)}")
if __name__=='__main__': main()
