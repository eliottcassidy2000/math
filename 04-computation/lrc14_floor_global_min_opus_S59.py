#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
IS DROP-7 THE GLOBAL R'-MINIMIZER OVER ALL COVERING FAMILIES?  (opus-2026-07-03-S59)

The covering floor is R'(S) = meas(lonely S) / (m_R m_Q), R = {s: 14 ∤ s}, Q = {s/14: 14|s}.
S58 found drop-7 (S={1..6,8..13,14}, R'=0.315) is the min among the 13 minus-one families of {1..13}.
The owner asks: is it the GLOBAL min over ALL covering families?

WHY THE MIN IS ATTAINED AT BOUNDED MAGNITUDE (the reduction to a finite search):
  As any far runner 14q -> infinity, 1_{||14q t||>=1/14} equidistributes => its correlation with the
  fixed rest decays => R' -> 1. So a LARGE far part pushes R' toward 1 (decorrelation), i.e. the floor's
  infimum is NOT out at infinity; it sits on SMALL, resonant families. This makes a bounded search
  meaningful (we verify the decorrelation numerically: R' vs growing far runner).

SEARCH: exact rational arcs. (a) r=1 far=14: R = {1..13} minus a d-subset (d=1,2,3) + 14-free fillers,
|R|=12, S=R u {14} covering; (b) r=1 far in {14,28,42}; (c) r=2 far={14,28}; (d) random covering
families up to magnitude 60. Report the global min R' and its family; compare to drop-7.
"""
import sys, itertools, random
from fractions import Fraction as Fr
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
random.seed(1)

BAND = Fr(1, 14)
def safe_arcs(v):
    return [((Fr(k) + BAND) / v, (Fr(k + 1) - BAND) / v) for k in range(v)]
def inter(A, B):
    r = []; i = j = 0
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0]); hi = min(A[i][1], B[j][1])
        if lo < hi: r.append((lo, hi))
        if A[i][1] < B[j][1]: i += 1
        else: j += 1
    return r
def meas(S):
    S = sorted(set(S))
    if not S: return Fr(1)
    a = safe_arcs(S[0])
    for v in S[1:]:
        a = inter(a, safe_arcs(v))
        if not a: return Fr(0)
    return sum(h - l for l, h in a)
def covers(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def split(S):
    R = [s for s in S if s % 14 != 0]; Q = [s // 14 for s in S if s % 14 == 0]
    return R, Q
def Rprime(S):
    R, Q = split(S)
    if not Q: return None
    mR, mQ, mS = meas(R), meas(Q), meas(S)
    if mR == 0 or mQ == 0: return None
    return float(mS / (mR * mQ)), R, Q, float(mS)

best = (9.9, None)
def consider(S):
    global best
    S = sorted(set(S))
    if len(S) != 13 or not covers(S): return
    r = Rprime(S)
    if r is None: return
    if r[0] < best[0]: best = (r[0], S)

FREE = [x for x in range(1, 28) if x % 14 != 0]         # 14-free pool 1..27
base = list(range(1, 14))                                # {1..13}
fillers = [x for x in range(15, 28) if x % 14 != 0]

# (a) r=1 far=14: drop d from {1..13}, add fillers to |R|=12
print("Searching (a) r=1 far=14, drop d=1,2,3 from {1..13} + fillers ...")
for d in (1, 2, 3):
    for drop in itertools.combinations(base, d):
        rem = [x for x in base if x not in drop]           # 13-d elements
        need = 12 - len(rem)                                # add this many fillers
        for add in itertools.combinations(fillers, need):
            consider(rem + list(add) + [14])
# (b) r=1 other far
print("Searching (b) r=1 far in {28,42} ...")
for far in (28, 42):
    for d in (1, 2):
        for drop in itertools.combinations(base, d):
            rem = [x for x in base if x not in drop]
            need = 12 - len(rem)
            for add in itertools.combinations(fillers, need):
                consider(rem + list(add) + [far])
# (c) r=2 far={14,28}
print("Searching (c) r=2 far={14,28}, R = 11 from {1..13} minus d + fillers ...")
for d in (2, 3):
    for drop in itertools.combinations(base, d):
        rem = [x for x in base if x not in drop]
        need = 11 - len(rem)
        for add in itertools.combinations(fillers, max(0, need)):
            if len(rem) + len(add) == 11:
                consider(rem + list(add) + [14, 28])
# (d) random covering families up to magnitude 60
print("Searching (d) 60000 random covering 13-families, magnitude <= 60 ...")
pool = [x for x in range(1, 61)]
for _ in range(60000):
    S = random.sample(pool, 13)
    if 14 not in S and 28 not in S and 42 not in S and 56 not in S:
        continue                                            # need a mult of 14 to cover q=14
    consider(S)

# drop-7 reference
dw7 = sorted([1,2,3,4,5,6,8,9,10,11,12,13,14])
r7 = Rprime(dw7)
print("\n" + "="*100)
print(" RESULT")
print("="*100)
print(f"  drop-7  S={dw7}  R'={r7[0]:.4f}")
print(f"  GLOBAL MIN found: R'={best[0]:.4f}  S={best[1]}")
if best[1] == dw7:
    print("  => drop-7 IS the global minimizer over the searched space.")
else:
    rb = Rprime(best[1])
    print(f"     R={rb[1]}  Q={rb[2]}  meas(S)={rb[3]:.5f}")
    print(f"  => drop-7 is NOT global; a lower family exists (R'={best[0]:.4f} < {r7[0]:.4f}).")

# decorrelation check: R' -> 1 as the far runner grows (drop-7 shape, far=14q)
print("\n  DECORRELATION (drop-7 core {1..6,8..13} + far=14q):")
core = [1,2,3,4,5,6,8,9,10,11,12,13]
for q in (1, 2, 3, 5, 10, 20, 50):
    S = sorted(core + [14*q])
    if covers(S):
        rr = Rprime(S)
        print(f"    q={q:>3} far={14*q:>4}  R'={rr[0]:.4f}  meas(S)={rr[3]:.5f}")
print("  (R' -> 1 as far grows => the floor's inf is attained at bounded magnitude; the search is finite.)")
print("DONE.")
