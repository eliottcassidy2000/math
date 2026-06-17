#!/usr/bin/env python3
"""
LRC(14) DISPROVE - TARGETED follow-up.
Focus on regimes most likely to break 1/14:
 (1) Neighbors of the sporadic tight {1..11,13,24} (and {1..13}).
 (2) "Danger-band alignment": choose speeds so that at tau=k/14 many elements
     land near a multiple of 14 (||v tau|| small). The lonely point of the tight
     AP is tau=1/14; we try to make EVERY tau in the candidate set produce some
     element with ||v tau|| < 1/14.
 (3) Random structured search over moderate-range primitive 13-sets, exact confirm.

SEED FROM PROVE: the prove route classifies M=1/14 minimizers (tight locus).
Every known minimizer is an AP-like set whose lonely point is tau=k/14.
To DISPROVE we need a set whose lower envelope never rises to 1/14 anywhere.
Idea borrowed: the minimizers are exactly where some v*tau hits exactly k/14
boundary; a true counterexample must AVOID having any tau where all ||v tau||>=1/14.
So we hunt for sets where the BEST tau is strictly below 1/14.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import random

def nrm(x):
    r = x - int(x)
    r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def g(S, t):
    return min(nrm(v * t) for v in S)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2):
            C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2))
    return C

def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b:
            b = v; at = t
    return b, at

def is_primitive(S):
    return reduce(gcd, S) == 1

THRESH = F(1, 14)
counterexamples = []
best_overall = [None]
seen = set()

def consider(S):
    S = tuple(sorted(set(int(x) for x in S)))
    if len(S) != 13 or any(v <= 0 for v in S): return
    if not is_primitive(S): return
    if S in seen: return
    seen.add(S)
    m, at = M(S)
    if best_overall[0] is None or m < best_overall[0][0]:
        best_overall[0] = (m, S, at)
    if m < THRESH:
        counterexamples.append((m, S, at))

# (1) Neighbors of sporadic tight configs
print("(1) Neighbors of tight {1..13} and sporadic {1..11,13,24}")
sporadics = [tuple(range(1,14)), tuple(list(range(1,12))+[13,24])]
for sp in sporadics:
    m, at = M(sp)
    print(f"   tight {sp}: M={m}={float(m):.6f} at {at}")
    L = list(sp)
    for i in range(13):
        for delta in range(-6, 30):
            v = L[i] + delta
            if v <= 0: continue
            consider(L[:i] + [v] + L[i+1:])
    # double neighbors
    for (i,j) in combinations(range(13),2):
        for di in (-2,-1,1,2):
            for dj in (-2,-1,1,2):
                vi, vj = L[i]+di, L[j]+dj
                if vi<=0 or vj<=0: continue
                NS = L[:]; NS[i]=vi; NS[j]=vj
                consider(NS)

# (2) Danger-band alignment around tau = k/14
print("(2) Danger-band alignment near tau=k/14")
# Pick speeds that are +-1, +-2, ... around multiples of 14 so that at tau=1/14
# v*tau = v/14 is near an integer => ||.|| small.
# Speeds near multiples of 14: {14a +- b}.
pools = []
for a in range(0,6):
    for b in range(-6,7):
        v = 14*a + b
        if v>0: pools.append(v)
pools = sorted(set(pools))
# also include the classic AP fillers
pools_all = sorted(set(pools + list(range(1,30))))
# sample 13-subsets that concentrate near 14k boundaries
random.seed(1234)
for _ in range(200000):
    S = random.sample(pools_all, 13)
    consider(S)

# (3) Random structured search, moderate range
print("(3) Random structured 13-sets in range 1..40")
for _ in range(300000):
    S = random.sample(range(1,41), 13)
    consider(S)

print()
print("="*70)
print(f"Counterexamples (M<1/14): {len(counterexamples)}")
counterexamples.sort()
for m,S,at in counterexamples[:20]:
    print(f"  M={m}={float(m):.6f} S={S} at {at}")
m,S,at = best_overall[0]
print(f"\nSMALLEST M overall: {m} = {float(m):.8f}  (1/14={float(THRESH):.8f})")
print(f"  S={S}, tau={at}")
print(f"  M - 1/14 = {m-THRESH} = {float(m-THRESH):.8f}")
print(f"  M >= 1/14 ? {m >= THRESH}")
