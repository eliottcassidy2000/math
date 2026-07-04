#!/usr/bin/env python3
"""
THOROUGH tight-locus enumeration (mac-mini-2026-07-03-S32): find ALL primitive tight families (M=1/14)
up to a speed bound, settle whether the locus is {AP,GW} or larger, characterize structure + gaps (g<=3?).
Given confinement (tight => q*=14, THM-612), phases on the 14th-root grid; a tight family's residues mod 14
cover the ±units {1,3,5,9,11,13}. Search:
 (a) AP[k->k+14m]: {1..13} with one speed k lifted to k+14, k+28 (residue-preserving high lift).
 (b) AP[k->j]: {1..13} with one speed replaced by a different small speed (residue change).
 (c) two-move families + broad numpy-prefiltered random. Classify tight ones by residue set + #gaps.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import numpy as np, random

def gcd_all(xs): return reduce(gcd, xs)
_G = np.arange(1, 8400)/8400.0
def approxM(sp):
    v = np.array(sp, dtype=np.float64); ph = np.outer(v, _G) % 1.0
    return np.minimum(ph, 1.0-ph).min(axis=0).max()
def nd(x):
    x = x % 1
    return x if x <= 1-x else 1-x
def M_exact(sp):
    vs = sorted(set(sp)); Q = set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1, len(vs)): Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    best = F(0)
    for q in Q:
        if q < 2: continue
        for a in range(1, q):
            m = min(nd(v*F(a,q)) for v in sp)
            if m > best: best = m
    return best
target = F(1,14); tf = 1.0/14
units14 = {1,3,5,9,11,13}
def gaps_of_residues(R):
    occ = sorted(set(R) | {0})
    gs = sorted(set((occ[(i+1)%len(occ)]-occ[i]) % 14 for i in range(len(occ))))
    return len(gs), gs
def is_tight(S):
    if len(set(S))!=13 or gcd_all(S)!=1: return False
    if abs(approxM(S)-tf) > 2e-4: return False
    return M_exact(S)==target
def classify(S):
    R = sorted(set(v%14 for v in S))
    g,gs = gaps_of_residues([v%14 for v in S])
    # canonical form up to dilation by units mod 14 (multiply residues by unit, sort)
    return tuple(R), g, gs

if __name__ == "__main__":
    import sys
    def out(*a): print(*a); sys.stdout.flush()
    AP = list(range(1,14))
    tights = {}   # canonical residue-set -> example
    def add(S):
        if is_tight(S):
            R,g,gs = classify(S)
            tights.setdefault((R,g), sorted(S))

    out("(a) AP[k -> k+14m]: replace one AP speed by a higher lift of the same residue")
    for k in range(1,14):
        for lift in [k+14, k+28, k+42]:
            add([x for x in AP if x!=k]+[lift])
    out(f"   tight so far: {len(tights)}")

    out("(b) AP[k -> j]: replace one AP speed by a different speed j<=56")
    for k in range(1,14):
        for j in range(14,57):
            if j in AP: continue
            add([x for x in AP if x!=k]+[j])
    out(f"   tight so far: {len(tights)}")

    out("(c) two moves: replace two AP speeds by lifts/replacements <=42")
    for k1,k2 in combinations(range(1,14),2):
        base=[x for x in AP if x not in (k1,k2)]
        for j1 in list(range(14,43)):
            for j2 in list(range(14,43)):
                if j1==j2 or j1 in base or j2 in base: continue
                add(base+[j1,j2])
    out(f"   tight so far: {len(tights)}")

    out("(d) broad numpy-prefiltered random (mixed, up to speed 60)")
    rng = random.Random(32); n=0
    for _ in range(120000):
        S = sorted(set(rng.sample(range(1,61),13)))
        if len(S)!=13: continue
        n+=1; add(S)
    out(f"   tested {n} random; tight so far: {len(tights)}")

    out("\n" + "="*80)
    out(f"DISTINCT tight residue-classes found: {len(tights)}")
    for (R,g), ex in sorted(tights.items()):
        name = "AP" if R==tuple(range(1,14)) else ("GW-type" if len(R)==12 else f"{len(R)}-residue")
        out(f"  residues={list(R)}  g={g}  [{name}]  example={ex}")
    gmax = max(g for (_,g),_ in tights.items()) if tights else 0
    out(f"\n  max #gaps g over all tight families found: {gmax}  (three-gap g(14)<=3 predicts <=3)")
    out(f"  => tight locus (by residue class): {'just AP+GW' if len(tights)<=2 else str(len(tights))+' classes (EXPANDED?)'}")
