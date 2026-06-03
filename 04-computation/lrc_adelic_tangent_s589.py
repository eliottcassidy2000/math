#!/usr/bin/env python3
"""Out-of-box tangent: LRC loneliness as ADELIC -- archimedean (cyclotomic, the static
unit-orbit rigidity) x 2-adic (doubling/Frobenius, the dynamical rigidity). The n=2q
obstruction = the 'bad place' p=2. Also a Mahler/height probe. opus-2026-06-03-S589."""
from fractions import Fraction as F
from math import gcd
def dist(x): x%=1; return min(x,1-x)
def main():
    print("ADELIC split of the n=14 AP worry-set (n=2*7):")
    n=14; AP=list(range(1,n))
    arch = [j for j in range(1,n) if gcd(j,n)==1]      # archimedean: primitive n-th roots = witnesses
    sevenths = [j for j in range(1,n) if gcd(j,n)==7]   # the q-place: 7th roots (even/14)
    twos = [j for j in range(1,n) if gcd(j,n)==2 and gcd(j,n)!=14]
    print(f"  clock j=1..13 stratified by gcd(j,14):")
    print(f"   archimedean WITNESSES (gcd=1, primitive 14th roots, lonely): {arch}")
    print(f"   7-place (gcd=7 -> j=7, the apex, NOT lonely): {[j for j in range(1,n) if gcd(j,n)==7]}")
    print(f"   2-place (gcd=2, the 7th-root shadow, NOT lonely): {[j for j in range(1,n) if gcd(j,n)==2]}")
    # 2-adic: Frobenius-at-2 = doubling; on the witnesses it FRAGMENTS (each odd->even)
    print(f"  Frobenius-at-2 (x2 mod 14) on the witnesses: ", {j:(2*j)%n for j in arch}, "(all land OUTSIDE the witnesses => 2 ramifies/fragments)")
    print(f"  cyclotomic field: Q(zeta_14)=Q(zeta_7) (since zeta_2=-1 rational) -- SAME as the prime 7!")
    print()
    print("Mahler/HEIGHT probe: the AP margin vs a product-of-phases at the witness t=1/n")
    for n in [7,12,14]:
        AP=list(range(1,n)); t=F(1,n)
        phases=[dist(v*t) for v in AP]
        mn=min(phases); prod=1.0
        for p in phases: prod*=float(p)*n   # normalized
        print(f"  n={n}: min phase (margin)={float(mn):.4f}=1/n; normalized product Prod(n*||v/n||)={prod:.3e}")
if __name__=='__main__': main()
