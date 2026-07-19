#!/usr/bin/env python3
"""Superseded five-row PG(2,13) scout (boxeph-S110; corrected codex-S67).

This script records a parameter analogy and five deterministic examples.  It does not
construct a Singer set, filter for Covering, run the formerly claimed random census, or
prove a global gap.  THM-1131 gives exact counterexamples to the original conclusion.
"""
from math import gcd
from fractions import Fraction as Fr
from collections import Counter


def Mstar(V, QMAX=1500):
    best = Fr(0)
    for q in range(2, QMAX + 1):
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            m = min(min((a * v) % q, q - ((a * v) % q)) for v in V)
            if Fr(m, q) > best:
                best = Fr(m, q)
    return best


def ddbl(C):
    C = list(C)
    return len(set(x - y for x in C for y in C))


def energy(C):
    C = list(C)
    c = Counter(x + y for x in C for y in C)
    return sum(v * v for v in c.values())


print(f"PG(2,13): points = {13**2 + 13 + 1} = 183 = deep-well modulus; line size = {13 + 1} = 14 = LRC(14)")
print("=> LRC(14) deep well and the Singer (183,14,1) difference set share parameters (183,14).\n")

print("the additive spectrum at (183,14): |C-C|, energy, and M for cores AP -> Sidon (+ far element):")
print(f"{'core':<22} {'|C-C|':>6} {'E(C)':>6} {'M':>14} {'M<1/13?':>8}")
cores = {
    'AP {1..12}': list(range(1, 13)),
    'near-AP {1..11,13}': list(range(1, 12)) + [13],
    'geometric-ish': [1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233],
    'Sidon-like': [1, 2, 5, 11, 22, 33, 45, 60, 78, 97, 120, 144],
    'powers 2^k': [2 ** k for k in range(12)],
}
for name, C in cores.items():
    V = C + [182] if 182 not in C else C + [364]
    M = Mstar(V)
    print(f"{name:<22} {ddbl(C):>6} {energy(C):>6} {str(M):>10}={float(M):.4f} {str(float(M) < 1 / 13):>8}")

print("\nCORRECTION: these five hand-picked rows do not imply a gap or an extremal theorem.")
print("Exact primitive Covering counterexamples: M=1/13, M=3/37 in (1/13,1/12),")
print("and the AP tower M=28/365<1/13; see THM-1131 and its exact pair-sum referee.")
print("No PG/Singer/tournament quotient preserving global M or Covering is supplied.")
