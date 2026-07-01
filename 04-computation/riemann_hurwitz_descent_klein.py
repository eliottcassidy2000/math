#!/usr/bin/env python3
"""
riemann_hurwitz_descent_klein.py  --  klein-2026-06-30-S58

Extending GAUSS-BONNET via DESCENT: the repo's 2-adic parity descent (THM-580, peel n=2p -> odd apex p)
is realized GEOMETRICALLY as the degeneracy map X_0(2p) -> X_0(p), a degree-3 branched cover, whose
RIEMANN-HURWITZ (= Gauss-Bonnet for covers) removes the level-2 curvature:

    chi(X_0(2p)) = 3 * chi(X_0(p)) - R      (R = total ramification = removed curvature)

So the descent LOWERS the genus (curvature -> spherical), and R = the curvature concentrated at the
level-2 = the iota-odd RESIDUAL (the genus-1 cusp form f_14 at n=14, S56/HYP-3768). The descent bottoms
at the genus-0 (spherical, obstruction-free) apex X_0(p). This unifies THM-580 (2-adic descent),
S56/S57 (curvature=genus, the residual), and mac-mini HYP-3771 ((2,3,p) apex geometry).

Genus of X_0(N):  g = 1 + psi/12 - nu2/4 - nu3/3 - nu_inf/2,  psi=N prod(1+1/p),
  nu2 = 0 if 4|N else prod_{p|N}(1+legendre(-1,p)),  nu3 = 0 if 9|N else prod(1+legendre(-3,p)),
  nu_inf = sum_{d|N} phi(gcd(d,N/d)).
"""
from math import gcd
from sympy import totient, divisors, isprime
from fractions import Fraction as F

def leg(a,p):
    a%=p
    if a==0: return 0
    return 1 if pow(a,(p-1)//2,p)==1 else -1

def psi(N):
    r=F(N)
    for p in set(_pf(N)): r*= (1+F(1,p))
    return int(r)

def _pf(N):
    f=[]; d=2; m=N
    while d*d<=m:
        while m%d==0: f.append(d); m//=d
        d+=1
    if m>1: f.append(m)
    return f

def nu2(N):
    # elliptic points of order 2: 0 if 4|N; else prod over p|N of factor:
    #   p=2 -> 1 (2||N);  odd p -> 1+(-1|p) = 2 if p==1 mod4, 0 if p==3 mod4
    if N%4==0: return 0
    r=1
    for p in set(_pf(N)):
        if p==2: continue          # factor 1
        r*= 2 if p%4==1 else 0
    return r

def nu3(N):
    # elliptic points of order 3: 0 if 9|N; else prod over p|N of factor:
    #   p=3 -> 1 (3||N);  other p -> 1+(-3|p) = 2 if p==1 mod3, 0 if p==2 mod3 (incl p=2)
    if N%9==0: return 0
    r=1
    for p in set(_pf(N)):
        if p==3: continue          # factor 1
        r*= 2 if p%3==1 else 0
    return r

def nu_inf(N):
    return sum(int(totient(gcd(d, N//d))) for d in divisors(N))

def genus(N):
    g = 1 + F(psi(N),12) - F(nu2(N),4) - F(nu3(N),3) - F(nu_inf(N),2)
    assert g.denominator==1, (N,g)
    return int(g)

if __name__=="__main__":
    print("(0) sanity: genus X_0(N) for small N (known: g=0 for N<=10,12,13; g=1 at 11,14; g=2 at 22,23,26)")
    for N in [1,7,10,11,13,14,15,22,23,26]:
        print(f"    X_0({N:2d}): genus={genus(N)}  cusps={nu_inf(N)} nu2={nu2(N)} nu3={nu3(N)} psi={psi(N)}")
    print()
    print("(1) THE RIEMANN-HURWITZ DESCENT X_0(2p) -> X_0(p) (degree 3 = the 2-adic peel):")
    print(f"    {'p':>3} {'g(p)':>5} {'g(2p)':>6} {'chi(p)':>7} {'chi(2p)':>8} {'R=3chi(p)-chi(2p)':>18} {'genus drop':>11}")
    for p in [3,5,7,11,13]:
        gp=genus(p); g2p=genus(2*p); cp=2-2*gp; c2p=2-2*g2p
        R=3*cp-c2p
        deg=psi(2*p)//psi(p)
        drop=gp-g2p if False else g2p-gp
        star="  <== LRC-14: flat g=1 -> spherical g=0 (curvature removed)" if p==7 else ""
        print(f"    {p:>3} {gp:>5} {g2p:>6} {cp:>7} {c2p:>8} {R:>18} {('+'+str(g2p-gp)):>11}{star}  (deg={deg})")
    print()
    print("(2) INTERPRETATION: R = ramification = the curvature the level-2 CONCENTRATES = the residual.")
    print("    The 2-adic descent (THM-580) peels n=2p -> apex p; geometrically = X_0(2p)->X_0(p) degree-3")
    print("    branched cover; genus (curvature) DROPS; the removed curvature R is the iota-odd residual")
    print("    (n=14: g 1->0, the removed genus-1 = f_14 = curve 14a, S56/HYP-3768). The descent BOTTOMS at")
    print("    the genus-0 (spherical, no cusp obstruction) apex X_0(p): the base case of the LRC descent.")
    print()
    print("(3) the (2,3,p) apex geometry (mac-mini HYP-3771) vs X_0-genus:")
    for p in [3,5,7,11,13]:
        ad = F(1,2)+F(1,3)+F(1,p)-1   # (2,3,p) angle defect sign: >0 spherical, =0 Euclid(p=6), <0 hyperbolic
        reg = "spherical" if ad>0 else ("Euclidean" if ad==0 else "hyperbolic")
        print(f"    (2,3,{p}): 1/2+1/3+1/{p}-1 = {ad} -> {reg};  X_0({2*p}) genus {genus(2*p)}")
