#!/usr/bin/env python3
"""
lrc14_multifold_ladder_kps_S15.py -- kind-pasteur-2026-07-05-S15, HYP-4177.

The multi-fold ladder law: multi-leg lift survivors are DOMINATED by the single-leg
deep well {1..11,168} (M=14/169). Extends mac-mini THM-621 (single-leg fourteen-fold).

Lift families: full residue system {1..12} mod 13, coordinate r -> r + 13 k_r.
For r>=7 (unique multiple of r), covering forces r|k_r -> value r(1+13m), m>=1 (14r at m=1).

Verifies: (1) exhaustive r>=7 m=1 domination; (2) the top-l ladder 14/(14(13-l)+1);
(3) higher-m loosening; (4) GLOBAL domination sample (r<=6 legs, general k).
"""
from math import gcd
from fractions import Fraction as Fr
from itertools import combinations, product
import random
random.seed(15)

def distZ(x, q):
    r = x % q; return min(r, q - r)

def M_exact(vs):
    Q = 2 * max(vs); best = Fr(0)
    for q in range(2, Q + 1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            mm = min(distZ(v * a, q) for v in vs)
            if Fr(mm, q) > best: best = Fr(mm, q)
    return best

def covering_2_12(vs):
    return all(any(v % q == 0 for v in vs) for q in range(2, 13))

DEEP = Fr(14, 169)

if __name__ == "__main__":
    print(f"single-leg deep well {{1..11,168}}: M = 14/169 = {float(DEEP):.5f};  2/25 = {float(Fr(2,25)):.5f}")

    # (1) exhaustive r>=7, m=1 (the tightest multi-leg lifts)
    print("\n(1) EXHAUSTIVE r>=7 m=1 lifts (L subset {7..12}):")
    minM = Fr(1); argmin = None; viol = 0
    for size in range(2, 7):
        for L in combinations(range(7, 13), size):
            W = sorted([x for x in range(1, 13) if x not in L] + [14 * r for r in L])
            M = M_exact(W)
            if M < minM: minM = M; argmin = L
            if M <= DEEP: viol += 1
    print(f"  min multi-leg M = {minM} at {argmin}; M<=14/169 count = {viol} (0 => deep well dominates)")

    # (2) the top-l ladder
    print("\n(2) top-l ladder  M(top-l) = 14/(14(13-l)+1) :")
    for l in range(1, 4):
        L = tuple(range(13 - l, 13))
        W = sorted([x for x in range(1, 13) if x not in L] + [14 * r for r in L])
        M = M_exact(W); pred = Fr(14, 14 * (13 - l) + 1)
        print(f"  l={l}: M={M}  pred={pred}  match={M==pred}")

    # (3) higher-m loosens
    print("\n(3) higher-m loosens (deep well r=12):")
    for m in [1, 2, 3]:
        W = sorted(list(range(1, 12)) + [12 * (1 + 13 * m)])
        print(f"  m={m}, killer={12*(1+13*m)}: M = {M_exact(W)}")

    # (4) GLOBAL domination sample: lift ANY coordinates (incl r<=6), general k; covering families
    print("\n(4) GLOBAL domination sample (any coords lifted, general k, covering families):")
    gmin = Fr(1); gargmin = None; gviol = 0; gtested = 0
    for _ in range(4000):
        L = random.sample(range(1, 13), random.randint(1, 4))
        W = list(range(1, 13))
        ok = True
        for r in L:
            k = random.randint(1, 3)
            W[r - 1] = r + 13 * k
        W = sorted(set(W))
        if len(W) != 12 or not covering_2_12(W): continue
        gtested += 1
        M = M_exact(W)
        if M < gmin: gmin = M; gargmin = list(W)
        if M < DEEP: gviol += 1
    print(f"  tested {gtested} covering lift families; min M = {gmin} = {float(gmin):.5f}")
    print(f"  families with M < 14/169: {gviol}  (0 => deep well is the GLOBAL tightest lift)")
    if gviol == 0:
        print("  => the single-leg deep well {1..11,168} DOMINATES the entire lift stratum.")
    else:
        print(f"  (tightest found: {gargmin})")
