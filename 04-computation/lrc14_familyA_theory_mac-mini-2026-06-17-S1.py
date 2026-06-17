#!/usr/bin/env python3
"""
LRC(14) PROVE — FAMILY (a) theory test, plus family (b) tight-subset test.

(a) S contains m0 = 14*j. The naive lonely point tau=1/14 puts ||m0*tau||=0 (BAD).
    But the AP tight config {1..13} achieves M=1/14 at tau=5/14 too. Key fact:
    if 14|v then ||v*(k/14)|| = 0 for ALL k. So ANY tau with denominator 14
    fails. The relevant lonely points for these S avoid denominator-14 rationals.

    HYPOTHESIS A: if 14*j in S then M(S) > 1/14 STRICTLY (never tight, never below).
    This would EXCLUDE family (a) entirely from the tight locus and from any
    counterexample. Test exhaustively on small ground sets + heavy random.

(b) A 3-subset {a,b,c} is "co-tight" if it has a common point tau* with
    ||a tau*||=||b tau*||=||c tau*||=1/14 (the AP triples do). Which 3-subsets
    of the AP share the lonely point tau=5/14? Characterize.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random, itertools

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

print("=== (a) HYPOTHESIS A: 14|v in S  =>  M(S) > 1/14 strictly ===")
# Heavy random test across wide ranges and multiple multiples of 14.
random.seed(7)
n_tight = 0; n_below = 0; n_tot = 0; minM = F(10); minS = None
for _ in range(20000):
    nmult = random.choice([1, 1, 1, 2])  # sometimes 2 multiples of 14
    mults = set()
    while len(mults) < nmult:
        mults.add(14 * random.randint(1, 8))
    rest = set()
    cap = random.choice([40, 70, 120, 200])
    while len(rest) < 13 - len(mults):
        x = random.randint(1, cap)
        if x not in mults: rest.add(x)
    S = sorted(rest | mults)
    if len(S) != 13: continue
    if gcd_all(S) != 1: continue
    n_tot += 1
    m, at = M(S)
    if m < THRESH:
        n_below += 1
        print("  !!! BELOW:", S, m, at)
    elif m == THRESH:
        n_tight += 1
        if n_tight <= 10: print("  TIGHT(=1/14):", S, "at", at)
    if m < minM: minM = m; minS = S
print(f"  total {n_tot}; below 1/14: {n_below}; exactly tight: {n_tight}")
print(f"  min M observed = {minM} = {float(minM):.6f}  at S={minS}")

print()
print("=== (a') Can M be tight if 14|v? Targeted: force AP-like + one mult of 14 ===")
# Replace ONE element of {1..13} by a multiple of 14 and see M.
base = list(range(1, 14))
tightcnt = 0
for i in range(13):
    for j in range(1, 16):
        mv = 14*j
        if mv in base[:i]+base[i+1:]: continue
        S = sorted(base[:i] + [mv] + base[i+1:])
        if len(set(S)) != 13: continue
        if gcd_all(S) != 1: continue
        m, at = M(S)
        flag = "TIGHT" if m == THRESH else ("BELOW" if m < THRESH else "")
        if flag:
            tightcnt += 1
            print(f"  i={i} mv={mv}: M={m} at {at}  {flag}  S={S}")
print(f"  tight-or-below count = {tightcnt}")

print()
print("=== (b) AP triples sharing the lonely point tau=5/14 ===")
# At tau=5/14, ||v*5/14|| for v=1..13:
t = F(5, 14)
vals = {v: nrm(v*t) for v in range(1, 14)}
ones = [v for v in range(1, 14) if vals[v] == THRESH]
print("  ||v*5/14|| == 1/14 for v in:", ones)
print("  full table:", {v: str(vals[v]) for v in range(1, 14)})
t2 = F(1, 14)
vals2 = {v: nrm(v*t2) for v in range(1, 14)}
ones2 = [v for v in range(1, 14) if vals2[v] == THRESH]
print("  at tau=1/14, ||v/14||==1/14 for v in:", ones2)
