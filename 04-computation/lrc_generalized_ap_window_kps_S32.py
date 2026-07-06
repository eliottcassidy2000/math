#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S32: the STRUCTURE x WIDTH residual at n=13 -- does the
shrinking window 1/325 admit any generalized-AP gap member?

opus HYP-4456: Freiman is necessary-not-sufficient.  The gap member at n=7 is the
GENERALIZED AP {1,5,6,11,16,17} (base AP {1,6,11,16} d=5, plus 5=6-1, 17=16+1),
which TILES (M=5/33 in (1/7,2/13)) with sub-max energy.  STRUCTURE narrows a gap
member to a generalized AP; the window WIDTH 1/((k+1)(2k+1)) ~ 1/(2k^2) decides
survival: k=6 window 1/91 ADMITS {1,5,6,11,16,17}; the OPEN residual is whether
the k=12 window 1/325 admits ANYTHING.

This censuses generalized-AP families on 12 runners (base AP + few perturbations,
= the n=7 gap-member SHAPE lifted) and checks whether any lands in the n=13 gap
(1/13, 2/25).  It also measures the M-RISE of each deficit family vs the window
width, quantifying opus's residual: M-rise < window <=> survives.
"""
from fractions import Fraction
from itertools import combinations, product
import numpy as np

LO, HI = Fraction(1, 13), Fraction(2, 25)
WINDOW = HI - LO   # = 1/325
def M_exact(v):
    v = [x for x in v if x != 0]
    S=int(sum(abs(x) for x in v)); Q=min(4*S,2*max(abs(x) for x in v)+2); va=np.array(v,dtype=np.int64); bn,bd=0,1
    for q in range(2,Q+1):
        a=np.arange(1,q,dtype=np.int64); r=np.outer(va,a)%q; d=np.minimum(r,q-r); bq=int(d.min(axis=0).max())
        if bq*bd>bn*q: bn,bd=bq,q
    return Fraction(bn,bd)
from math import gcd
from functools import reduce
def prim(v):
    g=reduce(gcd,[abs(x) for x in v]); return tuple(sorted(x//g for x in v))

print(f"gap = (1/13, 2/25); window width = {WINDOW} = 1/{1//WINDOW if WINDOW.numerator==1 else WINDOW.denominator} = {float(WINDOW):.5f}", flush=True)
print(f"(n=7 window 1/((7)(13)) = 1/91 = {1/91:.5f}; n=13 window ~3.6x narrower)", flush=True)

# ---- (0) verify the n=7 gap member ----
print("\n=== (0) verify the n=7 generalized-AP gap member {1,5,6,11,16,17} ===", flush=True)
n7 = [1,5,6,11,16,17]
M7 = M_exact(n7)
print(f"  M({n7}) = {M7} = {float(M7):.5f}; in (1/7,2/13)=({1/7:.4f},{2/13:.4f})? {Fraction(1,7)<M7<Fraction(2,13)}", flush=True)

# ---- (1) n=13 generalized-AP census: base AP (12 runners, various d) + perturbations ----
print("\n=== (1) n=13 generalized APs: base AP d + end perturbations; any in the gap (1/13,2/25)? ===", flush=True)
ingap=[]; nearest=Fraction(1)
count=0
# base AP: {a, a+d, ..., a+11d}; perturb up to 2 elements by +-1 (the n=7 shape: shift ends)
for d in range(1, 8):
    base = [1 + i*d for i in range(12)]
    # perturbation sets: change up to 2 runners by +-1 or +-2 (generalized-AP deficits)
    idxs = list(range(12))
    for npert in range(0, 3):
        for combo in combinations(idxs, npert):
            for deltas in product([-2,-1,1,2], repeat=npert):
                v = base.copy()
                for c,dl in zip(combo, deltas): v[c]+=dl
                if len(set(v)) != 12 or min(v) < 1: continue
                if reduce(gcd, v) != 1: continue
                count += 1
                M = M_exact(v)
                if LO < M < HI:
                    ingap.append((v, M))
                if abs(float(M)-float(LO)) < abs(float(nearest)-float(LO)) and M != LO and M > LO:
                    nearest = M
print(f"  censused {count} generalized-AP families (base d=1..7, up to 2 end perturbations)", flush=True)
print(f"  IN THE GAP (1/13,2/25): {len(ingap)}", flush=True)
for v,M in ingap[:8]:
    print(f"    *** {v} M={M}={float(M):.5f} ***", flush=True)
print(f"  nearest-to-1/13 (from above) M = {nearest} = {float(nearest):.5f}  (jump {float(nearest)-float(LO):.5f} vs window {float(WINDOW):.5f})", flush=True)

# ---- (2) the M-rise vs window: single-perturbation families (the smallest deficit) ----
print("\n=== (2) M-RISE of minimal-deficit AP {1..12}+one perturbation vs the window 1/325 ===", flush=True)
AP = list(range(1,13))
rises=[]
for c in range(12):
    for dl in [-2,-1,1,2]:
        v=AP.copy(); v[c]+=dl
        if len(set(v))!=12 or min(v)<1: continue
        M=M_exact(v)
        if M>LO: rises.append((float(M)-float(LO), v, float(M)))
rises.sort()
print(f"  smallest M-rises of AP+one-perturbation (all must exceed window {float(WINDOW):.5f} to be gap-empty):", flush=True)
for rise,v,M in rises[:8]:
    print(f"    {v}: M={M:.5f}, rise={rise:.5f} {'< window (would survive!)' if rise<float(WINDOW) else '> window (ejected)'}", flush=True)
print(flush=True)
print("VERDICT: if 0 in-gap AND every minimal-deficit M-rise > 1/325 => the n=13 window", flush=True)
print("is too narrow for ANY generalized-AP deficit -- opus's structure x width, at k=12.", flush=True)
