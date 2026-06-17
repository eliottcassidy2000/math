#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXACT lonely measure for LRC(14) (rational-endpoint arc sweep, not a grid), and an
exploration of inf L>0: scale-invariance L(cS)=L(S), the tight-AP L=0 degeneracy, the
interior-drop extremizer family's exact infimum, and a compactness/no-escape-to-zero probe.
kind-pasteur-2026-06-16-S5.

L(S) = meas{ τ∈[0,1) : ||v τ|| > 1/14 ∀v∈S }  (THM-515 lonely measure = the singular series).
inf L>0 over the relevant primitive 13-sets ⟹ C'(14) ⟹ LRC(14) (THM-398/501). inf≈0.0052 (cores).
DANGER set for speed v: ∪_k [ (14k-1)/(14v), (14k+1)/(14v) ]  (v arcs, each width 1/(7v), total 1/7).
L = 1 - meas(∪_v danger_v), computed EXACTLY with Fractions (sweep).
"""
import sys; sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import itertools, random

def danger_arcs(v):
    """arcs (as (lo,hi) Fractions in [0,1)) where ||v t||<=1/14; wraps handled by splitting at 0/1."""
    w=Fr(1,14*v)
    arcs=[]
    for k in range(v+1):
        c=Fr(k,v); lo=c-w; hi=c+w
        # clip to [0,1) with wraparound
        segs=[(lo,hi)]
        out=[]
        for a,b in segs:
            if b<=0 or a>=1:
                # shift into range
                a%=1; b=a+2*w  # width preserved
            out.append((a,b))
        for a,b in out:
            if a<0:
                arcs.append((Fr(0),b)); arcs.append((1+a,Fr(1)))
            elif b>1:
                arcs.append((a,Fr(1))); arcs.append((Fr(0),b-1))
            else:
                arcs.append((a,b))
    # merge dedupe later
    return arcs

def union_measure(arcs):
    arcs=sorted((a,b) for a,b in arcs if b>a)
    tot=Fr(0); cur_lo=None; cur_hi=None
    for a,b in arcs:
        if cur_hi is None: cur_lo,cur_hi=a,b
        elif a<=cur_hi: cur_hi=max(cur_hi,b)
        else: tot+=cur_hi-cur_lo; cur_lo,cur_hi=a,b
    if cur_hi is not None: tot+=cur_hi-cur_lo
    return tot

def L_exact(S):
    arcs=[]
    for v in S: arcs+=danger_arcs(v)
    return 1-union_measure(arcs)

def prim(S):
    g=reduce(gcd,S); return tuple(x//g for x in S)

# ---- (1) scale-invariance L(cS)=L(S) ----
print("=== (1) scale-invariance L(cS)=L(S) (exact) ===")
S0=(1,2,3,4,5,7,8,9,10,11,13,56)  # a near-tight core (12 elts here for speed of demo)
for c in (1,2,3,5,6):
    print(f"   c={c}: L({c}·S)= {float(L_exact(tuple(c*x for x in S0))):.8f}")
print("   ⟹ L depends only on the PRIMITIVE shape (proof: τ↦cτ is measure-preserving on the circle).")

# ---- (2) the tight AP {1..13}: L=0 (the excluded degeneracy) ----
ap=tuple(range(1,14))
print(f"\n=== (2) tight AP {{1..13}}: L = {L_exact(ap)} (exact) — the n=14 equality case (excluded) ===")

# ---- (3) interior-drop extremizer family ({1..13}\{j})∪{w}: exact infimum ----
print("\n=== (3) interior-drop cores ({1..13}\\{j}) ∪ {w}: EXACT L, sweep j and w ===")
best=(Fr(1),None)
for j in range(1,14):
    core=[x for x in range(1,14) if x!=j]
    row=[]
    for w in list(range(14,200))+[14*m for m in range(1,40)]:
        if w in core: continue
        S=tuple(sorted(core+[w]))
        if reduce(gcd,S)!=1: continue
        Lv=L_exact(S)
        if Lv>0 and Lv<best[0]: best=(Lv,S)
        if w in (14,28,56,84,98): row.append((w,float(Lv)))
    # report the j-family min over the swept w
print(f"   GLOBAL min over interior-drop family (j=1..13, w≤199 & 14m≤546): L={float(best[0]):.8f}")
print(f"      at S={best[1]}  (exact L={best[0]})")
# spotlight the known extremizers
for S in [(1,2,3,4,5,7,8,9,10,11,12,13,56),(1,2,3,4,5,6,7,8,9,10,11,13,84)]:
    print(f"   spotlight S={S}: L={float(L_exact(S)):.8f}")

# ---- (4) compactness probe: does any OTHER primitive 13-set go below the interior-drop inf? ----
print("\n=== (4) compactness probe: random + structured search for L below the interior-drop inf ===")
rng=random.Random(14); lo=(best[0],best[1]); nbelow=0; tested=0
# (a) two-element interior perturbations of the AP
for j1,j2 in itertools.combinations(range(1,14),2):
    core=[x for x in range(1,14) if x not in (j1,j2)]
    for w1 in range(14,90,2):
        for w2 in range(w1+1,140,2):
            S=tuple(sorted(core+[w1,w2]))
            if len(set(S))!=13 or reduce(gcd,S)!=1: continue
            Lv=L_exact(S); tested+=1
            if Lv>0 and Lv<lo[0]: lo=(Lv,S);
            if Lv>0 and Lv<best[0]-Fr(1,100000): nbelow+=1
# (b) random primitive 13-sets in [1,60]
for _ in range(3000):
    S=tuple(sorted(rng.sample(range(1,61),13)))
    if reduce(gcd,S)!=1: continue
    Lv=L_exact(S); tested+=1
    if Lv>0 and Lv<lo[0]: lo=(Lv,S)
print(f"   tested ~{tested} configs; lowest positive L found = {float(lo[0]):.8f} at S={lo[1]}")
print(f"   #(strictly below the interior-drop family inf {float(best[0]):.6f}): {nbelow}")
print("   (if 0 below and the global min stays ~0.005, evidence the inf is attained on a BOUNDED family ⟹ finite check)")
