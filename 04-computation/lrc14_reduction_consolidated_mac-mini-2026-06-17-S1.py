#!/usr/bin/env python3
"""
LRC(14) PROVE — CONSOLIDATED FAMILY REDUCTION (lean, flushed, fast).

Goal: shrink the open family of primitive 13-sets by exhibiting large
sub-families with M(S) >= 1/14 PROVABLE, and quantify the residual hard core.

All prints flushed. Trial counts kept modest so it finishes under contention.
"""
from fractions import Fraction as F
from math import gcd, ceil
from functools import reduce
import random, sys

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
def min_ratio(S):
    S = sorted(S); return min(F(S[i+1], S[i]) for i in range(len(S)-1))

THRESH = F(1, 14)

# sanity
m, at = M(list(range(1, 14)))
P("SANITY tight {1..13}: M =", m, "at", at, "(expect 1/14 at 5/14)")
assert m == THRESH

# ---------- FAMILY (a): contains a multiple of 14 ----------
P("\n=== (a) contains a multiple of 14 ===")
rng = random.Random(1)
below = 0; tight = 0; worst = F(10); wS = None; tot = 0
for _ in range(2000):
    mult = 14 * rng.randint(1, 6)
    rest = set()
    while len(rest) < 12:
        x = rng.randint(1, 80)
        if x != mult: rest.add(x)
    S = sorted(rest | {mult})
    if gcd_all(S) != 1: continue
    tot += 1
    mm, _ = M(S)
    if mm < THRESH: below += 1; P("  BELOW:", S, mm)
    elif mm == THRESH: tight += 1
    if mm < worst: worst = mm; wS = S
P(f"  tot={tot} below={below} tight={tight} worst={worst}={float(worst):.5f} wS={wS}")

# Exhaustive: replace ONE of {1..13} by a multiple of 14
P("  -- replace one AP speed by 14j --")
base = list(range(1, 14)); tb = 0
for i in range(13):
    for j in range(1, 12):
        mv = 14*j
        if mv in base[:i]+base[i+1:]: continue
        S = sorted(base[:i]+[mv]+base[i+1:])
        if len(set(S)) != 13 or gcd_all(S) != 1: continue
        mm, _ = M(S)
        if mm <= THRESH:
            tb += 1
            P(f"    i={i} mv={mv}: M={mm} {'TIGHT' if mm==THRESH else 'BELOW'} S={S}")
P(f"    tight-or-below among single-14j replacements: {tb}")

# ---------- FAMILY (c): spread out (ratio condition) ----------
P("\n=== (c) ratio v_{k+1}>=(7/6)v_k  =>  claim M>=1/14 ===")
rng = random.Random(11); viol = 0; tested = 0; worst = F(10); wS = None
for _ in range(3000):
    S = [rng.randint(1, 4)]
    for _ in range(12):
        lo = ceil(F(7, 6) * S[-1])
        S.append(lo + rng.randint(0, 6))
    S = sorted(set(S))
    if len(S) != 13 or gcd_all(S) != 1: continue
    if min_ratio(S) < F(7, 6): continue
    tested += 1
    mm, _ = M(S)
    if mm < THRESH: viol += 1; P("  RATIO-VIOLATION:", S, mm)
    if mm < worst: worst = mm; wS = S
P(f"  tested={tested} (ratio>=7/6) violations={viol} worst={worst}={float(worst):.5f}")

# Threshold scan: worst M per ratio band
P("\n=== (c') worst M per min-ratio band ===")
rng = random.Random(13); bands = {}
for _ in range(15000):
    S = sorted(rng.sample(range(1, 70), 13))
    if gcd_all(S) != 1: continue
    r = min_ratio(S)
    if r < 1: continue
    mm, _ = M(S)
    bi = int(r * 20)
    cur = bands.get(bi)
    if cur is None or mm < cur[0]: bands[bi] = (mm, S, r)
P("  ratio>=   worstM    float    safe?")
for bi in sorted(bands):
    mm, S, r = bands[bi]
    P(f"  {float(F(bi,20)):.2f}     {str(mm):>8} {float(mm):.5f}  {'YES' if mm>=THRESH else 'NO'}")

# ---------- HARD CORE: perturbations of the AP ----------
P("\n=== HARD CORE: single-speed perturbation of {1..13} ===")
base = list(range(1, 14)); tight_perturb = []; below_perturb = []; nt = 0
for i in range(13):
    for newv in range(1, 80):
        if newv in base[:i]+base[i+1:]: continue
        S = sorted(base[:i]+[newv]+base[i+1:])
        if len(set(S)) != 13 or gcd_all(S) != 1: continue
        nt += 1
        mm, at = M(S)
        if mm < THRESH: below_perturb.append((S, mm, at))
        elif mm == THRESH: tight_perturb.append((tuple(S), at))
P(f"  tested {nt} single-replacements; tight(=1/14): {len(tight_perturb)}; below: {len(below_perturb)}")
for S, at in tight_perturb:
    P("    TIGHT:", list(S), "at", at)
for S, mm, at in below_perturb:
    P("    BELOW!:", S, mm, at)

# ---------- lonely-point structure at tau=5/14 ----------
P("\n=== lonely-point structure ===")
for t in (F(5,14), F(1,14), F(3,14)):
    vals = {v: nrm(v*t) for v in range(1, 14)}
    ones = [v for v in range(1, 14) if vals[v] == THRESH]
    P(f"  tau={t}: speeds v with ||v*tau||=1/14: {ones}")
P("DONE", flush=True)
