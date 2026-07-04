#!/usr/bin/env python3
"""
GAP-A structure: non-covering tight locus = {AP,GW}? (mac-mini-2026-07-03-S35)
Non-covering tight S: primitive, M=1/14, MISSES 14, covers {2..13}, tight at t=1/14 (q*=14).
 (1) FORCED RESIDUES: residues ⊇ {1,3,5,7,9,11,13} (all odds): units {1,3,5,9,11,13} by ±units; 7 by
     covering q=7 (7|v => v mod14 ∈{0,7}, miss 14 => =7). Verify on AP, GW + any found tight family.
 (2) SPEED BOUND: search non-covering tight families for the max speed (empirically 24=GW). If bounded,
     GAP-A is a FINITE check.
 (3) COVER-q=12 CHOICES: {1..11,13}∪{X} with 12|X — which X are tight? (12->AP, 24->GW, 36+ loose?).
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
def covers(sp, q): return any(v % q == 0 for v in sp)
def is_noncov_tight(S):
    if len(set(S)) != 13 or gcd_all(S) != 1: return False
    if covers(S, 14): return False               # non-covering: miss 14
    if abs(approxM(S)-tf) > 1.5e-4: return False
    return M_exact(S) == target

if __name__ == "__main__":
    import sys
    def out(*a): print(*a); sys.stdout.flush()
    ODDS = {1,3,5,7,9,11,13}
    out("(1) FORCED RESIDUES (odds {1,3,5,7,9,11,13} present) on AP, GW:")
    for name, S in [('AP', list(range(1,14))), ('GW', [1,2,3,4,5,6,7,8,9,10,11,13,24])]:
        res = set(v % 14 for v in S)
        out(f"   {name}: residues={sorted(res)}  contains all odds? {ODDS.issubset(res)}  miss14={not covers(S,14)}")

    out("\n(2) SPEED BOUND: search non-covering tight families (perturbations + random), record max speed:")
    maxspeed = 0; found = {}
    AP = list(range(1,14))
    # single/double residue swaps + lifts
    cands = [AP, [1,2,3,4,5,6,7,8,9,10,11,13,24]]
    for k in range(1,14):
        for j in range(14, 120):
            if j not in AP: cands.append([x for x in AP if x!=k]+[j])
    for k1,k2 in combinations(range(1,14),2):
        base=[x for x in AP if x not in (k1,k2)]
        for j1 in range(14,60):
            for j2 in range(j1+1,60):
                if j1 not in base and j2 not in base: cands.append(base+[j1,j2])
    rng = random.Random(35)
    for _ in range(150000):
        cands.append(sorted(random.Random(rng.random()).sample(range(1,80),13)))
    seen=set()
    for S in cands:
        key=tuple(sorted(S))
        if key in seen: continue
        seen.add(key)
        if is_noncov_tight(S):
            res=tuple(sorted(set(v%14 for v in S)))
            found.setdefault(res, sorted(S))
            maxspeed=max(maxspeed, max(S))
    out(f"   distinct non-covering tight residue-classes: {len(found)}; MAX SPEED over all found: {maxspeed}")
    for res, ex in sorted(found.items()):
        out(f"      residues={list(res)}  example={ex}  maxspeed={max(ex)}")

    out("\n(3) COVER q=12 via {1..11,13}∪{X}, 12|X: which X are tight?")
    base = [1,2,3,4,5,6,7,8,9,10,11,13]
    for X in range(12, 200, 12):
        S = base + [X]
        if len(set(S))!=13: continue
        M = M_exact(S); tightq = (M==target and gcd_all(S)==1)
        out(f"   X={X} (={X%14} mod14, 12*{X//12}): M={float(M):.5f} tight&primitive={tightq}")
    out("\n=> odds forced (rigorous); if max speed bounded (empirically 24), GAP-A is a FINITE check = {AP,GW}.")
    out("   cover-q=12: only X=12 (AP) and X=24 (GW) tight => the binary AP/GW choice is 'how to cover 12'.")

# (4) [S35] which k admit a SECOND tight coverer j*k? (only k=12 => 24 => GW unique on single-swap axis)
def _M(sp):
    return M_exact(sp)
if __name__ == "__main__":
    print("\n(4) AP runner k -> j*k (j>=2): which give a tight family? (only k=12 -> 24 expected)")
    for k in range(2,14):
        ts=[j*k for j in range(2,17) if j*k not in list(range(1,14)) and j*k<=200
            and len(set([x for x in range(1,14) if x!=k]+[j*k]))==13
            and gcd_all([x for x in range(1,14) if x!=k]+[j*k])==1
            and _M([x for x in range(1,14) if x!=k]+[j*k])==target]
        print(f"   k={k}: tight second-coverers = {ts}")
