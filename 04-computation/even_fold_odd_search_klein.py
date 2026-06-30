#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""EVEN-FOLD + smarter ODD-n covering-min search + connections (klein-S35)."""
import math, itertools
from fractions import Fraction as F

def gap_exact(S, Dmax):
    best=(F(0),0,0)
    for D in range(2,Dmax+1):
        bj=0; ba=0
        for a in range(1,D):
            m=D
            for s in S:
                r=(s*a)%D; d=r if r<=D-r else D-r
                if d<m: m=d
                if m<=bj: break
            if m>bj: bj=m; ba=a
        if F(bj,D)>best[0]: best=(F(bj,D),D,ba)
    return best

print("="*84); print(" (1) VERIFY the odd-n spread beaters + examine structure (binding modulus)"); print("="*84)
beaters={7:[1,2,5,6,7,8], 9:[1,3,4,5,7,11,18,32]}
for n,S in beaters.items():
    M,D,a=gap_exact(S, 2*max(S)+2)
    img=sorted((s*a)%D for s in S)
    print(f"  n={n}: S={S} -> M={M}={float(M):.5f}, binding D={D}, rotation a={a}")
    print(f"     speeds*a mod D = {img}  (avoid the {M.numerator}-nbhd of 0; covering radius {M.numerator})")
    # structure: is S*a mod D an interval/spread? differences?
    print(f"     gaps of speeds*a on Z/{D}: {sorted((img[(i+1)%len(img)]-img[i])%D for i in range(len(img)))}")

print("\n"+"="*84); print(" (2) THE EVEN-FOLD: descent of the equally-spaced {1..n-1} for even n=2p -> {1..p-1}"); print("="*84)
def descend(S):
    O=[v for v in S if v%2==1]; E=[v for v in S if v%2==0]
    return O, sorted(set(v//2 for v in E))
for n in [6,8,10,14]:
    S=list(range(1,n)); O,Sp=descend(S)
    M,_,_=gap_exact(S,2*n)
    Mp,_,_=gap_exact(Sp,2*max(Sp)+2) if Sp else (F(0),0,0)
    print(f"  n={n}: equally-spaced {{1..{n-1}}} -> odd part O={O}, E/2={Sp} (= {{1..{n//2-1}}}, the p={n//2} equally-spaced!); M(n)={M}={float(M):.4f}, M(p)={Mp}={float(Mp):.4f}")
print("  => the even-fold descends the equally-spaced {1..2p-1} to the equally-spaced {1..p-1}: the worst-case")
print("     STRUCTURE is preserved by the 2-adic descent. Even n=2p folds to p (Mode B); p<=7 is the proven LRC range.")

print("\n"+"="*84); print(" (3) EVEN vs ODD cycle: bipartite (apex gap=0) vs non-bipartite (apex gap>0)"); print("="*84)
import cmath
for p in [6,7,8,9]:
    w=cmath.exp(2j*math.pi/p)
    # signless Laplacian Q(C_p)=2I+A(C_p), lambda_min = 2-2cos(pi/p) if p odd, 0 if p even (bipartite)
    lam_min=min(2+2*math.cos(2*math.pi*k/p) for k in range(p))
    print(f"  C_{p}: bipartite={p%2==0}; lambda_min(2I+A(C_p))={lam_min:.5f} ({'=0 (even cycle, apex gap VANISHES)' if p%2==0 else '>0 (odd cycle, apex gap = the obstruction)'})")
print("  => EVEN n: the even cycle C_n is BIPARTITE -> apex gap 0 -> the equally-spaced is the degenerate")
print("     cusp (M=1/n, measure 0). ODD n: C_n non-bipartite -> apex gap>0 -> the spread/realizability frontier.")
print("     This is the same even/odd = bipartite/non-bipartite = sigma-even/sigma-odd split (HYP-3604/3728);")
print("     the odd cycle is the tournament's Condorcet 3-cycle atom (HYP-3602).")
