#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THE LRC->2D REDUCTION: torus lift + three-distance gaps (klein-S28).

(A) CF(n/Phi_6(n)) = [0; n-1, n] (since n^2-n+1 = (n-1)n+1).
(B) THREE-DISTANCE: the binding config (runner positions at t*=a/D) has <=3 distinct gaps; M=(largest gap)/2/D.
(C) TORUS LIFT: Eisenstein-embed the speeds into the hexagonal lattice Z[w]/(n-w) (w=zeta_6, |a+bw|^2=a^2+ab+b^2);
    see the 2D geometry (line? hexagonal?) and the deep hole.
"""
import math
from fractions import Fraction as F

def cf(p,q):
    out=[]
    while q: out.append(p//q); p,q=q,p%q
    return out
def phi6(n): return n*n-n+1

print("="*78); print(" (A) CF(n/Phi_6(n)) = [0; n-1, n]"); print("="*78)
for n in [3,5,7,14]:
    P=phi6(n); print(f"   n={n}: {n}/{P} -> CF {cf(n,P)}   (= [0, n-1={n-1}, n={n}])")

print("\n"+"="*78); print(" (B) THREE-DISTANCE: binding config {1..12,182}*14 mod 183 -- the gaps"); print("="*78)
D=183; a=14; S=list(range(1,13))+[182]
pos=sorted((s*a)%D for s in S)
gaps=[]
for i in range(len(pos)):
    nxt = pos[(i+1)%len(pos)]
    g = (nxt - pos[i]) % D
    gaps.append(g)
from collections import Counter
gc=Counter(gaps)
print(f"   runner positions: {pos}")
print(f"   gaps (Counter): {dict(sorted(gc.items()))}  -> {len(gc)} distinct (three-distance: <=3)")
maxgap=max(gaps); print(f"   largest gap = {maxgap}; observer 0 sits in it; M = (largest gap)/2 / D = {maxgap//2}/{D} = {float(F(maxgap//2,D)):.5f}")
print(f"   (= 14/183, the covering-min; the deep hole is the half of the size-{maxgap} gap around 0)")

print("\n"+"="*78); print(" (C) TORUS LIFT: Eisenstein reps of the speeds in Z[w]/(14-w), |a+bw|^2=a^2+ab+b^2"); print("="*78)
def eis_rep(r, D=183, k=14, R=20):
    best=None
    for b in range(-R,R+1):
        for aa in range(-R,R+1):
            if (aa + k*b) % D == r%D:
                nrm=aa*aa+aa*b+b*b
                if best is None or nrm<best[0]: best=(nrm,aa,b)
    return best
print("   speed -> Eisenstein rep (a,b), norm:")
online=0
for s in S:
    nrm,aa,bb=eis_rep(s)
    if bb==0: online+=1
    print(f"     {s:>3} -> ({aa:>3},{bb:>2})  norm {nrm:>3}{'  [on the a-axis line, b=0]' if bb==0 else ''}")
print(f"   => {online}/{len(S)} speeds lie on the a-axis (b=0): the speed set is a LINE in the hexagonal lattice.")
print(f"   so the covering-min is a 1D LINE in the 2D hexagonal lattice; the three-gap (B) governs the line,")
print(f"   and the 3 gaps are the line's view of the 3 hexagonal nearest-neighbor directions (the A2 3-fold).")
