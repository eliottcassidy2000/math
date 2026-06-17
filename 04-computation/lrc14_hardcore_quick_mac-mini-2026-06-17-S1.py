#!/usr/bin/env python3
"""LRC(14) HARD CORE quick: exhaustive single-speed perturbations of {1..13}.
Plus 2-speed perturbations (drop 2 AP speeds, add 2 new in a small window).
Identify EXACTLY which perturbed configs remain tight (M=1/14) or go below.
Output flushed line by line."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools

def P(*a): print(*a, flush=True)
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def g(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at
def gcd_all(S): return reduce(gcd, S, 0)

THRESH = F(1, 14)
base = list(range(1, 14))

P("=== single-speed perturbation of {1..13}, new speed up to 100 ===")
tight = []; below = []; nt = 0
for i in range(13):
    for newv in range(1, 101):
        if newv in base[:i]+base[i+1:]: continue
        S = sorted(base[:i]+[newv]+base[i+1:])
        if len(set(S)) != 13 or gcd_all(S) != 1: continue
        nt += 1
        mm, at = M(S)
        if mm < THRESH: below.append((S, mm, at))
        elif mm == THRESH: tight.append((tuple(S), at, i, newv))
P(f"  tested {nt}; tight(=1/14): {len(tight)}; below: {len(below)}")
for S, at, i, newv in tight:
    P(f"  TIGHT: drop {base[i]} add {newv}: {list(S)} at {at}")
for S, mm, at in below:
    P(f"  BELOW!: {S} M={mm} at {at}")

P("\n=== two-speed perturbation: drop 2 AP speeds, add 2 new in [14..28] ===")
tight2 = 0; below2 = 0; n2 = 0; examples = []
addpool = list(range(14, 29))
for drop in itertools.combinations(range(13), 2):
    keep = [base[k] for k in range(13) if k not in drop]
    for add in itertools.combinations(addpool, 2):
        S = sorted(set(keep) | set(add))
        if len(S) != 13 or gcd_all(S) != 1: continue
        n2 += 1
        mm, at = M(S)
        if mm < THRESH:
            below2 += 1; P(f"  BELOW2!: {S} M={mm} at {at}")
        elif mm == THRESH:
            tight2 += 1
            if len(examples) < 25: examples.append((S, at))
P(f"  tested {n2}; tight: {tight2}; below: {below2}")
for S, at in examples:
    P(f"  TIGHT2: {S} at {at}")

P("\n=== FAMILY (b): which 3-subsets of AP share lonely point tau=5/14 ===")
t = F(5, 14)
tightspeeds = [v for v in range(1, 14) if nrm(v*t) == THRESH]
P(f"  at tau=5/14, ||v*tau||=1/14 for v in {tightspeeds}")
# The full 13-set is tight because EVERY residue mod 14 is hit. Show residues:
res = {v: (5*v) % 14 for v in range(1, 14)}
P(f"  residues 5v mod 14: {res}")
P("DONE", flush=True)
