#!/usr/bin/env python3
"""
LRC(14) PROVE route, part 2: locate the HARD CORE.

SEED FROM DISPROVE: the hard core (M near 1/14) = perturbations of the tight AP {1..13}.
Everything far from the AP is generic (M large) => provably safe.

We quantify:
 (1) How does M behave as we perturb {1..13} by replacing/shifting a few speeds?
 (2) Is the set of "dangerous" 13-sets (M close to 1/14) bounded up to dilation?
 (3) Family (b): does a tight 3-subset force tightness? Test which 3-subsets
     have a common lonely point.
 (4) A SPREAD threshold: if speeds grow fast enough, M >= 1/14 provably.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools, random

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

print("=== (1) SINGLE-SPEED replacement of {1..13} ===")
# Replace one speed by a new value up to some cap; track M.
dangerous = []  # M == 1/14 exactly (still tight)
below = []
nval = 0
for i in range(13):
    for newv in range(1, 60):
        if newv in base[:i] + base[i+1:]: continue
        S = sorted(base[:i] + [newv] + base[i+1:])
        if len(set(S)) != 13: continue
        if gcd_all(S) != 1: continue
        nval += 1
        m, at = M(S)
        if m < THRESH: below.append((S, m, at))
        elif m == THRESH: dangerous.append((tuple(S), at))
print(f"  tested {nval} single-replacements; below 1/14: {len(below)}; "
      f"exactly-tight (M=1/14): {len(dangerous)}")
for S, at in dangerous[:40]:
    print("    TIGHT:", list(S), "at tau=", at)
if below:
    print("  !!! COUNTEREXAMPLES (M<1/14):")
    for S, m, at in below: print("    ", S, m, at)

print()
print("=== (2) DILATION INVARIANCE check: M(cS)=M(S) ===")
S0 = list(range(1, 14))
for c in [1, 2, 3, 5, 7]:
    Sc = [c*v for v in S0]
    print(f"  c={c}: gcd={gcd_all(Sc)}  M={M(Sc)[0]}")
# Note: scaling makes gcd=c (not primitive) but M invariant — confirms primitivity
# is the right normalization and dangerous family is bounded UP TO DILATION.

print()
print("=== (4) SPREAD THRESHOLD: lacunary speeds v_{k+1} >= rho*v_k ===")
# If speeds grow geometrically with ratio rho, when is M >= 1/14 guaranteed?
for rho in [F(2), F(5,2), F(3)]:
    worst = F(10); wS = None
    rng = random.Random(int(rho*100))
    for _ in range(2000):
        S = [1]
        ok = True
        for _ in range(12):
            lo = (rho*S[-1]).__ceil__()
            nxt = lo + rng.randint(0, 5)
            S.append(nxt)
        S = sorted(set(S))
        if len(S) != 13: continue
        if gcd_all(S) != 1: continue
        m, _ = M(S)
        if m < worst: worst = m; wS = S
    print(f"  rho={rho}: worst M over 2000 = {worst}={float(worst):.5f} "
          f"(>=1/14? {worst>=THRESH})  wS={wS}")
