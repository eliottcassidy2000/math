#!/usr/bin/env python3
"""
klein-2026-07-04-S127 (HYP-4082) - THE ONE-SWAP COVERING STRATUM: is the deep well its minimum,
and does every drop-j family close by a kps-style residue formula?

kps-S4 closed the "residue-liar" {1..11,13,12k} = {1..13}\{12} u {12k}: M=k/(12k+5), residue
table, Lean (LRCResidueLiar.lean).  The covering-min extremizer is the deep well {1..12,182} =
{1..13}\{13} u {182} (drop-13).  These are two members of the ONE-SWAP family
    F(j, X) = ({1,...,13} \ {j}) u {X}   (drop runner j, add a large X).
This maps the WHOLE one-swap stratum: for each drop position j, find the covering completions X,
compute M(F(j,X)) exactly, and confirm (a) the covering-min over one-swap families is the deep well
14/183, and (b) each near-min family sits on a ladder M = k/(a k + b) (kps-formula shape),
=> the one-swap stratum is a finite union of formula-closable ladders, deep well the floor.

Exact (Fractions).
"""
from fractions import Fraction as F
from math import gcd
import itertools

def cdist_q(a, q):
    r = a % q
    return min(r, q - r)

def Mval(S, Qcap):
    best = F(0)
    for Q in range(2, Qcap + 1):
        for a in range(1, Q // 2 + 1):
            if gcd(a, Q) != 1: continue
            m = min(F(cdist_q(v * a, Q), Q) for v in S)
            if m > best: best = m
    return best

def is_covering(S, n=14):
    return all(any(v % q == 0 for v in S) for q in range(2, n + 1))

DW = F(14, 183)
print(f"deep well target: 14/183 = {float(DW):.6f}")
print("="*78)
print("ONE-SWAP COVERING FAMILIES F(j,X) = ({1..13}\\{j}) u {X}: min M over covering X")
print("  (X up to 260; report the covering-min per drop-j and the (j,X) achieving global min)")
print("="*78)
glob_min = None; glob_at = None
print(f"{'drop j':>6} {'min covering M':>16} {'>=14/183?':>10} {'argmin X':>9} {'#cov':>5}")
for j in range(1, 14):
    base = [v for v in range(1, 14) if v != j]
    best = None; bestX = None; cov = 0
    for X in range(14, 261):
        if X in base: continue
        S = sorted(base + [X])
        if not is_covering(S): continue
        cov += 1
        M = Mval(S, min(2 * X + 2, 540))
        if best is None or M < best: best, bestX = M, X
    if best is not None:
        ge = best >= DW
        print(f"{j:>6} {str(best):>16} {str(ge):>10} {bestX:>9} {cov:>5}")
        if glob_min is None or best < glob_min: glob_min, glob_at = best, (j, bestX)
print(f"\n  GLOBAL one-swap covering-min = {glob_min} (~{float(glob_min):.6f}) at drop-j,X = {glob_at}")
print(f"  equals deep well 14/183? {glob_min == DW}   (deep well = drop-13, X=182)")

print()
print("="*78)
print("LADDER SHAPE per drop-j: M(F(j, X_j0 * k)) as k grows -- is it k/(a k + b) (kps shape)?")
print("  X_j0 = smallest covering X for drop-j; tabulate M for k=1..5, fit linear/linear.")
print("="*78)
print(f"{'j':>3} {'X0':>5} {'M(k=1..5)':>52}")
for j in range(1, 14):
    base = [v for v in range(1, 14) if v != j]
    # smallest covering X
    X0 = None
    for X in range(14, 400):
        if X in base: continue
        if is_covering(sorted(base + [X])): X0 = X; break
    if X0 is None:
        print(f"{j:>3} {'--':>5} (no covering X<=400)"); continue
    Ms = []
    for k in range(1, 6):
        X = X0 * k
        if X in base: Ms.append("dup"); continue
        S = sorted(base + [X])
        Ms.append(str(Mval(S, min(2 * X + 2, 900))))
    print(f"{j:>3} {X0:>5} {str(Ms):>52}")

print()
print("READING: if the global one-swap covering-min = 14/183 (deep well) and each drop-j ladder is")
print("M=k/(ak+b) >= 14/183, the ONE-SWAP stratum is a finite (13 drop-positions) union of formula-")
print("closable ladders (kps residue-table), with the deep well the unique floor. That closes the")
print("one-runner-spread part of the covering-min; residual = MULTI-swap families (>=2 spread).")
print("DONE")
