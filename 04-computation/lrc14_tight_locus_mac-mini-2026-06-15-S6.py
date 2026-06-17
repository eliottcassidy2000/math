#!/usr/bin/env python3
"""
LRC(14) DISPROVE-via-known-hard-instances + tight-locus enumeration.
mac-mini-2026-06-15-S6.

For a 13-set S of distinct positive integers, M(S)=max_tau min_v ||v tau||.
LRC(14) <=> M(S) >= 1/14 for all primitive S.  Counterexample := primitive 13-set with M<1/14.

Strategy: compute EXACT M for all classically-tight / near-tight families:
  consecutive AP {1..13}, dilations, {1..11,13,x}, {1..12,x}, {1..10,a,b},
  general APs, 1-element perturbations of {1..13}, plus a randomized screen.
Report any M<1/14 (counterexample) and the COMPLETE tight locus (M==1/14, primitive).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1,2) else 1 - r

def g(S, t):
    return min(nrm(v*t) for v in S)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1,2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1,2): C.add(F(k, d)); k += 1
    C.add(F(1,2)); return C

def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at

def gf(S, t):  # float lower-bound screen of g at t
    m = 2.0
    for v in S:
        r = (v*t) % 1.0
        r = r if r <= 0.5 else 1.0 - r
        if r < m: m = r
    return m

def Mfloat(S):  # float screen of M (max over candidate float taus)
    best = 0.0
    Sl = sorted(set(S)); n = len(Sl)
    for v in Sl:
        k = 0
        while (2*k+1)/(2.0*v) <= 0.5:
            val = gf(S, (2*k+1)/(2.0*v))
            if val > best: best = val
            k += 1
    for i in range(n):
        for j in range(i+1, n):
            for d in (Sl[i]+Sl[j], Sl[j]-Sl[i]):
                if d > 0:
                    k = 1
                    while k/float(d) <= 0.5:
                        val = gf(S, k/float(d))
                        if val > best: best = val
                        k += 1
    return best

def is_prim(S): return reduce(gcd, S) == 1
THR = F(1,14); THRf = 1.0/14.0

below = []; tight = set()
def test(S, lab):
    S = sorted(set(S))
    if len(S) != 13 or min(S) < 1: return
    if Mfloat(S) < THRf + 1e-9:   # only confirm exactly when float says it could be <= 1/14
        m, at = M(S)
        if m < THR: below.append((tuple(S), m, at, lab))
        elif m == THR and is_prim(S): tight.add(tuple(S))

# A) consecutive AP and dilations
test(range(1,14), "consec")
# B) {1..11,13,x}, {1..12,x}
for x in range(14, 500):
    test(list(range(1,12))+[13,x], "A")
    test(list(range(1,13))+[x], "B")
# C) {1..10,a,b}
for a in range(11, 50):
    for b in range(a+1, 100):
        test(list(range(1,11))+[a,b], "C")
# D) general APs
for a in range(1, 40):
    for d in range(1, 40):
        test([a+i*d for i in range(13)], "AP")
# E) 1-element perturbations of {1..13}
base = list(range(1,14))
for i in range(13):
    for nv in range(1, 120):
        test(base[:i]+base[i+1:]+[nv], "perturb")
# F) randomized primitive screen
random.seed(1)
RAND = 120000; cnt = 0
for _ in range(RAND):
    S = random.sample(range(1, 55), 13)
    if not is_prim(S): continue
    cnt += 1
    if Mfloat(S) < THRf - 1e-9:
        m, at = M(S)
        if m < THR: below.append((tuple(S), m, at, "random"))

print("random primitive tested:", cnt)
print("=== M < 1/14 COUNTEREXAMPLES ===")
for b in below: print(b)
if not below: print("  NONE FOUND")
print()
print("=== tight primitive (M == 1/14) ===")
for S in sorted(tight):
    m, at = M(S)
    print(S, " tau =", at)
print("total distinct tight primitive:", len(tight))
