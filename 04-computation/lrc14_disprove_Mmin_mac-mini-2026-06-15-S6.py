#!/usr/bin/env python3
"""
LRC(14) DISPROVE via broad M-minimization (FIXED: bounded, flushing, float-screened).
Search primitive 13-sets with M(S) < 1/14 across structured families.
Exact rationals for confirmation; float screen for breadth.

Families:
 (a) single/double perturbations of tight AP {1..13}
 (b) dilated/shifted APs d*{1..13}+r, AP variants, unions of two APs
 (c) sets with many coincident danger bands (divisors of HCN; multiples)
 (d) sporadic tight {1..11,13,24} and neighbors
SEED FROM PROVE: tight configs hit M=1/14 with lonely tau in {1/14,5/14}.
Deliberately destroy that lonely point.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import random, sys

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

def Mfloat(S):
    """fast float upper-confirm of best tau value; returns float best M."""
    best = 0.0
    Sf = list(S)
    # build float candidate set quickly
    Ss = sorted(set(S)); maxd = 0
    Cs = set([0.5])
    for v in Ss:
        k = 0
        while (2*k+1) <= v:  # (2k+1)/(2v) <= 1/2
            Cs.add((2*k+1)/(2.0*v)); k += 1
    for i in range(len(Ss)):
        for j in range(i+1, len(Ss)):
            for d in (Ss[i]+Ss[j], Ss[j]-Ss[i]):
                if d > 0:
                    k = 1
                    while 2*k <= d:  # k/d <= 1/2
                        Cs.add(k/float(d)); k += 1
    for tf in Cs:
        m = 1.0
        for v in Sf:
            x = v*tf
            r = x - int(x)
            if r < 0: r += 1
            dd = r if r <= 0.5 else 1-r
            if dd < m: m = dd
            if m <= best: break
        if m > best: best = m
    return best

def is_primitive(S):
    return reduce(gcd, S) == 1

THRESH = F(1, 14)
THRf = 1.0/14.0
counterexamples = []
best_overall = [None]
best_float = [1.0, None]   # smallest float M, set
seen = set()

def consider(S):
    """Float screen first; exact confirm only if float M is suspiciously low."""
    S = tuple(sorted(set(int(x) for x in S)))
    if len(S) != 13 or any(v <= 0 for v in S): return
    if not is_primitive(S): return
    if S in seen: return
    seen.add(S)
    mf = Mfloat(S)
    if mf < best_float[0]:
        best_float[0] = mf; best_float[1] = S
    # confirm exactly when within a margin of 1/14 (float noise safety)
    if mf < THRf + 1e-9:
        m, at = M(S)
        if best_overall[0] is None or m < best_overall[0][0]:
            best_overall[0] = (m, S, at)
        if m < THRESH:
            counterexamples.append((m, S, at))

print("Baseline tight {1..13}:", flush=True)
m, at = M(tuple(range(1,14)))
print(f"  M={m}={float(m):.6f} at {at}; 1/14={THRf:.6f}", flush=True)

# (a) perturbations of {1..13}
print("(a) single+double perturbations of {1..13}", flush=True)
base = list(range(1,14))
for i in range(13):
    for newv in range(1,60):
        S = base[:i]+base[i+1:]
        if newv in S: continue
        consider(S+[newv])
for (i,j) in combinations(range(13),2):
    rest = [base[k] for k in range(13) if k not in (i,j)]
    for (a,b) in combinations(range(1,40),2):
        if a in rest or b in rest: continue
        consider(rest+[a,b])
print(f"  done. counterexamples={len(counterexamples)} bestM={best_overall[0]}", flush=True)

# seed-from-prove: destroy lonely point
print("(seed) destroy lonely tau in {1/14,5/14}", flush=True)
for repl in range(13):
    rest = base[:repl]+base[repl+1:]
    for v in [14,28,42,15,29,27,41,56,70,13,84,98]:
        if v in rest: continue
        consider(rest+[v])

# (b) dilated/shifted APs
print("(b) dilated/shifted APs + AP variants + unions", flush=True)
for d in range(1,9):
    for r in range(0,50):
        S = [d*k+r for k in range(1,14)]
        if all(v>0 for v in S): consider(S)
for start in range(1,40):
    for step in range(1,15):
        consider([start+step*k for k in range(13)])
for len1 in range(2,12):
    len2 = 13-len1
    for s1 in range(1,8):
        for st1 in range(1,6):
            A = [s1+st1*k for k in range(len1)]
            for s2 in range(1,30):
                for st2 in range(1,6):
                    B = [s2+st2*k for k in range(len2)]
                    Sset = set(A)|set(B)
                    if len(Sset)==13: consider(Sset)
print(f"  done. counterexamples={len(counterexamples)} bestM={best_overall[0]}", flush=True)

# (c) small-lcm / many shared divisors (BOUNDED sampling)
print("(c) divisors of HCN + multiples (bounded random sampling)", flush=True)
random.seed(7)
for N in [60,120,180,240,360,720,840,1260,2520]:
    divs = [d for d in range(1,N+1) if N % d == 0]
    if len(divs) >= 13:
        for _ in range(20000):
            consider(random.sample(divs,13))
for base_mult in range(1,6):
    for extra in range(1,40):
        consider([base_mult*k for k in range(1,13)]+[extra])
print(f"  done. counterexamples={len(counterexamples)} bestM={best_overall[0]}", flush=True)

# (d) sporadic tight neighbors
print("(d) sporadic {1..11,13,24} neighbors", flush=True)
sp = list(range(1,12))+[13,24]
m, at = M(tuple(sp))
print(f"  tight {tuple(sp)}: M={m}={float(m):.6f} at {at}", flush=True)
for i in range(13):
    for delta in range(-6,30):
        v = sp[i]+delta
        if v<=0: continue
        consider(sp[:i]+[v]+sp[i+1:])

# danger-band random structured sweep near multiples of 14
print("(e) danger-band random sweep near multiples of 14", flush=True)
pool = sorted(set([14*a+b for a in range(0,6) for b in range(-6,7) if 14*a+b>0] + list(range(1,30))))
random.seed(99)
for _ in range(200000):
    consider(random.sample(pool,13))
# broad random in 1..40
print("(f) broad random 1..40", flush=True)
for _ in range(300000):
    consider(random.sample(range(1,41),13))

print("="*70, flush=True)
print(f"Counterexamples (M<1/14): {len(counterexamples)}", flush=True)
counterexamples.sort()
for m,S,at in counterexamples[:25]:
    print(f"  M={m}={float(m):.6f} S={S} at {at}", flush=True)
if best_overall[0]:
    m,S,at = best_overall[0]
    print(f"\nSMALLEST EXACT M (among float-screened): {m}={float(m):.8f}", flush=True)
    print(f"  S={S} tau={at}", flush=True)
    print(f"  M-1/14={m-THRESH}={float(m-THRESH):.8f}; M>=1/14? {m>=THRESH}", flush=True)
print(f"\nSmallest FLOAT M seen: {best_float[0]:.8f} at S={best_float[1]} (1/14={THRf:.8f})", flush=True)
# exact-confirm the float champion
if best_float[1]:
    mm,att = M(best_float[1])
    print(f"  exact M of float champion: {mm}={float(mm):.8f} at {att}; >=1/14? {mm>=THRESH}", flush=True)
