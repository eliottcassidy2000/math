#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
BOUNDING v_max(U) -- the even part (THM-612/Lemma D residual for f=2, m=2 confinement). opus-S62.

mac-mini reduced f=2 confinement to a FINITE per-U check (tighteners bounded by even-part data), leaving
"bound v_max(U)" as the residual, called "a cleaner target". THIS SCRIPT TESTS whether it is cleaner or
whether it IS the AP-rigidity again:

For an even part E=2U (11 runners) we compute the BEST f=2 "tightening gap"
    gap(U) = min over odd pairs (w1,w2) (bounded) of ( M(2U u {w1,w2}) - 1/14 ),
i.e. how close two odd tighteners get S to tight (gap=0 <=> tight q*=28 exists). mac-mini: gap>0 always
(0 hits). QUESTION: as u_max grows, does gap -> 0 (near-feasible, no uniform u_max bound => NOT cleaner)
or stay bounded away (=> u_max effectively bounded)? And is gap SMALLEST for the dilated-AP-like U
(=> bounding u_max = the tight-locus AP rigidity, NOT a cleaner target)?
"""
import sys, itertools
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

import numpy as np
def M_float(S, G=120019):
    t = (np.arange(G) + 0.5) / G
    F = np.full(G, 1.0)
    for v in S:
        fr = (v * t) % 1.0
        F = np.minimum(F, np.minimum(fr, 1.0 - fr))
    return F.max()
def norm(x):
    x = x - int(x)
    if x < 0: x += 1
    return min(x, 1 - x)
def exact_M(S):
    S = sorted(set(S)); cands = set()
    for v in S:
        for k in range(v): cands.add(Fr(2*k+1, 2*v))
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for den in (S[i]+S[j], abs(S[i]-S[j])):
                if den:
                    for s in range(den): cands.add(Fr(s, den))
    best = Fr(0)
    for t in cands:
        val = min(norm(v*t) for v in S)
        if val > best: best = val
    return best
def prim(S): return reduce(gcd, S) == 1
TARGET = 1.0/14
TARGET_F = Fr(1,14)

def best_gap(Uprime, wmax):
    """min over odd pairs w1<w2<=wmax of (M(2U' u {w1,w2}) - 1/14), among M>=1/14. Float search, exact-verify best."""
    E = [2*u for u in Uprime]
    odds = [w for w in range(1, wmax+1, 2) if w not in E]
    best = (1.0, None)
    for w1, w2 in itertools.combinations(odds, 2):
        S = sorted(E + [w1, w2])
        if not prim(S): continue
        M = M_float(S)
        if M >= TARGET - 1e-6:
            g = M - TARGET
            if g < best[0]: best = (g, (w1, w2, M))
    # exact re-verify the winner
    if best[1]:
        w1, w2, _ = best[1]
        Me = exact_M(sorted(E + [w1, w2]))
        best = (float(Me - TARGET_F), (w1, w2, float(Me)))
    return best

print("="*104)
print(" f=2 tightening gap = min_{odd w1,w2} (M(2U u {w1,w2}) - 1/14).  gap=0 <=> tight q*=28 exists.")
print("="*104)
# representative 11-runner even parts U' (E=2U'): dilated-AP-like vs generic, small vs larger u_max
Us = [
 ("AP-like {1..11}",              list(range(1,12))),
 ("AP-like {1..13}\\{6,12}",      [x for x in range(1,14) if x not in (6,12)]),
 ("AP-like {2..12}",              list(range(2,13))),
 ("generic {1,2,3,5,7,8,9,10,11,12,13}", [1,2,3,5,7,8,9,10,11,12,13]),
 ("generic mid {3,4,5,7,8,9,10,11,13,15,17}", [3,4,5,7,8,9,10,11,13,15,17]),
 ("larger umax {1,2,3,4,5,7,9,11,13,17,23}",  [1,2,3,4,5,7,9,11,13,17,23]),
 ("larger umax {1,2,4,5,8,10,13,16,19,22,25}",[1,2,4,5,8,10,13,16,19,22,25]),
]
print(f"  {'even part U (E=2U)':<40} {'u_max':>6} {'gap=min(M-1/14)':>16} {'argmin (w1,w2,M)':>22}")
rows=[]
for name, Up in Us:
    wmax = min(85, 13*max(Up))            # tighteners (odd) up to the compactness region
    g, arg = best_gap(Up, wmax)
    rows.append((name, max(Up), g, arg))
    print(f"  {name:<40} {max(Up):>6} {g:>16.6f} {str(arg):>22}")

print("\n"+"="*104); print(" READING"); print("="*104)
apgaps=[g for n,u,g,a in rows if "AP-like" in n]
gngaps=[g for n,u,g,a in rows if "generic" in n or "larger" in n]
print(f"  AP-like gaps: {[round(x,4) for x in apgaps]}")
print(f"  generic/larger-umax gaps: {[round(x,4) for x in gngaps]}")
print("  * gap=0 for none => 0 tight q*=28 (confirms confinement on this slice).")
print("  * If AP-like U have the SMALLEST gaps (closest to feasible) => bounding u_max = the tight-locus")
print("    AP rigidity (the near-feasible threats ARE the dilated APs, which are imprimitive) => NOT a")
print("    cleaner target than the rigidity. If generic large-umax U have comparably small gaps => u_max")
print("    genuinely needs bounding (a distinct, possibly cleaner, target).")
print("DONE.")
