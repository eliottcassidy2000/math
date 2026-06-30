#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THE LRC->HEXAGONAL BRIDGE + SUBCONDITIONS (klein-S26, fixed).

M(S) = lonely-runner GAP = max over rational t=a/D of min_s ||s a/D|| = j/D at the BINDING modulus D
(HYP-3551). Covering-min 14/183 = n/Phi_6(n) at {1..12,182}; D=Phi_6(n)=Eisenstein norm; Z/D=hexagonal
quotient. Subconditions: dense cores (where the min lives), antipodal binding (1,-1), the lcm killer,
the hexagonal embedding, discrete Kershner (most-spread minimizes covering radius).
"""
from math import gcd
from fractions import Fraction as F

def gap_exact(S, Dmax):
    """max over (D,a), 2<=D<=Dmax, gcd(a,D)=1, of min_{s in S} dist(s*a mod D, 0)/D. Returns (M, D, a, j)."""
    best=(F(0),0,0,0)
    Sl=sorted(set(s % 1 if False else s for s in S))
    for D in range(2, Dmax+1):
        bj=0
        for a in range(1, D):
            m=D
            for s in Sl:
                r=(s*a)%D; d=r if r<=D-r else D-r
                if d<m: m=d
                if m<=bj: break
            if m>bj: bj=m; ba=a
        val=F(bj,D)
        if val>best[0]: best=(val,D,ba,bj)
    return best

print("="*82); print(" (1) M = gap = covering radius j/D: reproduce 14/183 and 7/89"); print("="*82)
for name,S in [("{1..12,182}", list(range(1,13))+[182]),
               ("{1..11,13,84}", list(range(1,12))+[13,84])]:
    M,D,a,j = gap_exact(S, 2*max(S)+2)
    print(f"   {name:<16} M = {j}/{D} = {float(M):.5f}  (binding rot a={a});  D=Phi_6(14)=183? {D==183}")

print("\n"+"="*82); print(" (2) DENSEST-CORE-WINS: {1..13}\\{skip} + minimal killer lcm(skip,14); the covering-min"); print("="*82)
def lcm(a,b): return a*b//gcd(a,b)
rows=[]
for skip in range(1,14):
    core=[x for x in range(1,14) if x!=skip]; killer=lcm(skip,14); S=core+[killer]
    M,D,a,j=gap_exact(S, 2*max(S)+2); rows.append((float(M),skip,j,D,killer))
for Mf,skip,j,D,killer in sorted(rows):
    tag=" <== covering-min" if abs(Mf-14/183)<1e-6 else ""
    print(f"   skip {skip:>2} (killer {killer:>3}): M={j}/{D}={Mf:.5f}{tag}")
print("   => densest coverable core {1..12} (skip 13) + killer 182=lcm(13,14) -> 14/183 (HYP-3551).")

print("\n"+"="*82); print(" (3) HEXAGONAL embedding + antipodal binding"); print("="*82)
D=183
print(f"   zeta_6=14 mod 183: 14^6 mod 183={pow(14,6,D)} (order 6); 183=Phi_6(14)=13*14+1 (Eisenstein norm).")
print(f"   binding pair (1,182): 182=-1 mod 183 (ANTIPODAL); the gap 14/183 = hex-quotient covering radius/D.")

print("\n"+"="*82); print(" (4) DISCRETE KERSHNER: min covering radius for k speeds in Z/D (tightest covering)"); print("="*82)
import itertools
def cov_rad(S,D): return max(min(min((s*a)%D,D-(s*a)%D) for s in S) for a in range(1,D))
for D,k in [(13,4),(21,5),(31,6)]:
    best=(10**9,None)
    for S in itertools.combinations(range(1,D),k):
        j=cov_rad(S,D)
        if j<best[0]: best=(j,S)
    print(f"   Z/{D}, k={k}: min covering radius={best[0]} (~D/2k={D/(2*k):.1f}); a minimizer={best[1]}")
print("   => the tightest covering uses MAXIMALLY-SPREAD speeds (the difference-set/hexagonal pole). For the")
print("      LRC covering-min the spread is forced into Z/Phi_6(n) (hexagonal), where Kershner is optimal.")
