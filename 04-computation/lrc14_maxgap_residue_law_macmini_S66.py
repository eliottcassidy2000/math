#!/usr/bin/env python3
"""mac-mini-S66 (corrected): THE MAX-GAP RESIDUE LAW.

THEOREM (exact, elementary). For co-offsets E and x=a/b in lowest terms,
    maxgap{frac(e_i x) : i} = G_E(a,b) / b,
where G_E(a,b) = the max circular gap between consecutive elements of the occupied residue
set R = {a*e_i mod b} in Z/b (gap = residue-distance to next occupied residue). So:
    x=a/b is a GOOD period center  <=>  G_E(a,b) > 2b/7.
This reframes good-period existence (THM-527-A) as an ARITHMETIC residue-gap search over
moduli b -- exact, finite per b, no analysis. The arc around a/b has half-width
    w(a,b) = (G_E(a,b)/b - 2/7) / (2 * s(a,b)),  s(a,b)=max drift = max_i|e_i - e_center|,
so FAT arcs (a ruler integer j can land in one <=> w > 1/(2 Vmax)) come from SMALL b / small
drift. GOAL: (1) verify the residue law exactly; (2) for each cluster, find the modulus b
giving the FATTEST good arc, and confirm it is small/bounded."""
from fractions import Fraction as F
from math import gcd

def maxgap_exact(E,x):
    ph=sorted(set((F(e)*x)%1 for e in E)); n=len(ph)
    if n<=1: return F(1)
    return max((ph[(i+1)%n]-ph[i]) if i<n-1 else (ph[0]+1-ph[n-1]) for i in range(n))

def G_residue(E,a,b):
    R=sorted(set((a*e)%b for e in E)); m=len(R)
    if m==0: return b
    return max((R[(i+1)%m]-R[i]) if i<m-1 else (R[0]+b-R[m-1]) for i in range(m))

# (1) verify the law on random (E,a,b)
import random; random.seed(3)
print("(1) VERIFY maxgap{frac(e x)} == G_E(a,b)/b at x=a/b (exact):")
bad=0
for _ in range(3000):
    k=random.randint(3,13); E=random.sample(range(0,40),k)
    b=random.randint(2,40); a=random.randint(1,b-1)
    if gcd(a,b)!=1: continue
    lhs=maxgap_exact(E,F(a,b)); rhs=F(G_residue(E,a,b),b)
    if lhs!=rhs: bad+=1
print(f"    3000 random trials: mismatches = {bad}  ({'RESIDUE LAW CONFIRMED' if bad==0 else 'FAILED'})")

# (2) fattest good-arc modulus per cluster
def fattest_good(E, Bmax=80):
    D=max(E)-min(E)
    best=None  # (halfwidth, a, b, G)
    for b in range(2,Bmax+1):
        for a in range(1,b):
            if gcd(a,b)!=1: continue
            G=G_residue(E,a,b)
            if F(G,b)>F(2,7):
                # drift bound: as x moves off a/b, phase i moves at rate e_i; the binding
                # gap edges move at rate = spread of the two edge speeds <= D. half-width:
                hw=(F(G,b)-F(2,7))/(2*max(D,1))
                if best is None or hw>best[0]: best=(hw,a,b,G)
    return best,D

clusters = {
 "consec k=9 {0..8}": list(range(9)),
 "consec k=13 {0..12} (AP,OUT OF SCOPE)": list(range(13)),
 "perforated k=7 {0..8}\\{1,7}": [0,2,3,4,5,6,8],
 "spread-21 k=11": [0,2,4,6,8,10,12,14,16,18,21],
 "deep-well cluster {0..12}@core": list(range(13)),
 "k=8 block {0..7}": list(range(8)),
}
print("\n(2) fattest good arc per cluster (small b => fat => ruler-hittable at moderate Vmax):")
print(f"{'cluster':40s} | spread | fattest (a/b, G) | half-width | Vmax needed (1/2hw)")
print("-"*104)
for nm,E in clusters.items():
    best,D=fattest_good(E)
    if best:
        hw,a,b,G=best
        vneed=int(1/(2*float(hw)))+1 if hw>0 else -1
        print(f"{nm:40s} | {D:6d} | {a}/{b} G={G} ({G}/{b}={float(F(G,b)):.3f}) | {float(hw):.5f} | ~{vneed}")
    else:
        print(f"{nm:40s} | {D:6d} | NONE within Bmax  <-- no good period at all")
print("\nINTERPRETATION: the fattest good arc uses SMALL b (b<=~2*spread). Its half-width")
print("hw=(G/b-2/7)/(2*spread) is bounded below when b is small and G/b comfortably >2/7, so a")
print("ruler period j lands in it once Vmax > 1/(2 hw) -- an EXPLICIT V0 from the residue law.")
print("The good-period existence (THM-527-A) = 'some modulus b has residue-gap G>2b/7 with hw>0'.")
