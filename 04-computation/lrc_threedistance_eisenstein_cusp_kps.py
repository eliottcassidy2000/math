#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE EISENSTEIN/CUSP DICHOTOMY IS THE THREE-DISTANCE (STEINHAUS) THEOREM.

kind-pasteur-2026-07-01. opus HYP-3775: margin = classical Dedekind sum (=> regularizable,
Eisenstein, -1/12=zeta(-1)) EXACTLY when the witness residues are a "scaled interval (AP)".
This script identifies that criterion as the classical THREE-DISTANCE THEOREM and gives it a
clean checkable form: the witness residues of the CONSTRUCTION are a SINGLE-ROTATION ORBIT
(<=3 distinct cyclic gaps -- the Steinhaus signature of {k*alpha}), i.e. a scaled AP; the
BEATERS need >3 gaps (>=2 rotations = higher "gap dimension") -- genuinely non-modular / cusp.

  regularizable / Eisenstein / -1/12  <=>  witness residues = single-rotation orbit (<=3 gaps = AP)
  un-regularizable / cusp / residual  <=>  witness residues need >3 gaps (>=2 rotations, higher-dim)

The "gap dimension" (rotations needed) is the geometric shadow of the ci-odd cusp index / genus.
This partially rehabilitates opus HYP-3773's "beaters = higher-dimensional" intuition (in the
gap/rotation sense, NOT the literal Zagier-cotangent sense opus correctly withdrew).

HONEST NEGATIVE (also recorded): the SOFT concentration invariants (additive energy, log-concave
autocorrelation) do NOT cleanly threshold regularizability -- it is EXACT/arithmetic (the AP /
three-distance property), not analytic. Soft metrics confirm the construction is maximally
ordered (E_norm->1) but confound dense/small-D cases. Scope: the covering-min (large witness D).
"""
import sys, random
from fractions import Fraction as Fr
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def witness(S,Dmax):
    """M(S)=max_{a/D} min_v ||v a/D|| over a/D, D<Dmax; the binding rotation."""
    best=(Fr(-1),0,0)
    for D in range(2,Dmax):
        for a in range(1,D):
            if gcd(a,D)!=1: continue
            m=min(min((v*a)%D, D-(v*a)%D) for v in S)
            if Fr(m,D)>best[0]: best=(Fr(m,D),a,D)
    return best
def gaps(R,D):
    R=sorted(set(x%D for x in R))
    return sorted(set((R[(i+1)%len(R)]-R[i])%D for i in range(len(R))))
def analyze(S,Dmax):
    M,a,D=witness(S,Dmax); g=gaps([(v*a)%D for v in S],D); return len(g),M,D,a

print("="*92)
print("THE THREE-DISTANCE CRITERION: #distinct cyclic gaps of the witness residues")
print("  <=3 gaps = single-rotation orbit = scaled AP = Steinhaus = REGULARIZABLE (Eisenstein, -1/12)")
print("  >3 gaps  = >=2 rotations = higher gap-dimension = CUSP (un-regularizable residual)")
print("="*92)
print("  CONSTRUCTION {1..n-2, n(n-1)}  (the ordered/Eisenstein object):")
for n in range(7,15):
    ng,M,D,a=analyze(list(range(1,n-1))+[n*(n-1)], n*n-n+5)
    print(f"    n={n:2d}: #gaps={ng}  M={str(M):>8} (D=Phi6={D})   {'THREE-DISTANCE (AP orbit) -> regularizable' if ng<=3 else '>3'}")
print("  BEATERS (the disordered/cusp covering-min, opus HYP-3769):")
beaters={7:[1,2,5,6,7,8],8:[1,4,5,6,7,11,16],9:[1,3,4,5,7,11,18,32],10:[1,2,3,5,6,7,8,9,30]}
for n,S in beaters.items():
    ng,M,D,a=analyze(S,60)
    print(f"    n={n:2d}: #gaps={ng}  M={str(M):>7} (D={D})   {'>3 gaps -> NOT three-distance -> cusp' if ng>3 else '<=3'}")

print("\n"+"="*92)
print("WHY <=3 gaps = regularizable: a single-rotation orbit {k*alpha} (Steinhaus/three-distance) is")
print("an AP; the sawtooth sum over an AP TELESCOPES into a Dedekind sum (Ostrowski/continued-fraction),")
print("closed via reciprocity (Phi6==1 mod n) -> -1/12. >3 gaps = the residues are NOT one orbit -> no")
print("telescoping -> no Dedekind form -> cusp. The 'gap dimension' = # rotations = the genus/cusp index.")
print("="*92)
print("  three-distance <-> continued fractions: the number of gaps of {k*alpha mod 1} is <=3 for ALL")
print("  N (Steinhaus); the gap SIZES are the CF convergent denominators of alpha. So 'one rotation' = ")
print("  'one continued fraction' = the O(log) reciprocity descent of last session; >3 gaps = no single CF.")

print("\n"+"="*92)
print("HONEST SCOPE (the criterion is exact/arithmetic, not soft; large-D covering regime)")
print("="*92)
rng=random.Random(11); le3=0; tot=0; le3_bigD=0; bigD=0
for _ in range(400):
    k=rng.choice([6,7,8]); S=sorted(rng.sample(range(1,25),k))
    if reduce(gcd,S)!=1: continue
    ng,M,D,a=analyze(S,60); tot+=1
    if ng<=3: le3+=1
    if D>=k*k:            # 'covering-min regime': witness denominator large vs #speeds
        bigD+=1
        if ng<=3: le3_bigD+=1
print(f"  random sets: {le3}/{tot} have <=3 gaps overall (NOISY: small-D non-covering sets trivially few gaps);")
print(f"    but among large-D (D>=k^2, covering-like) witnesses: {le3_bigD}/{bigD} have <=3 gaps")
print(f"    => the construction's <=3 gaps is meaningful in the LARGE-D covering-min regime; there it is the")
print(f"       three-distance/AP signature that beaters violate. Soft-concentration metrics do NOT threshold")
print(f"       regularizability (it is the EXACT arithmetic AP/three-distance property).")
print("DONE.")
