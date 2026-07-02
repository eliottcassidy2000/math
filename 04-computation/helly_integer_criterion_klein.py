#!/usr/bin/env python3
"""
klein-2026-07-01-S93 -- HYP-4000: the 1D-HELLY INTEGER CRITERION (the Lean-ready
restatement of THM-601's verification layer).

CRITERION (for gcd-reduced triple u < v < w at r = 1/14): a non-origin triple-danger
point exists iff there are integers (a,b,c), 0 <= a <= u, 0 <= b <= v, 0 <= c <= w,
(a,b,c) not in {(0,0,0), (u,v,w)}, with the three PAIRWISE inequalities
    14|av - bu| <= u+v-1,   14|bw - cv| <= v+w-1,   14|aw - cu| <= u+w-1.
(1D Helly: three open intervals |x - a/u| < r/u, |x - b/v| < r/v, |x - c/w| < r/w share
a point iff pairwise intersecting; pairwise intersecting <=> the integer inequalities.)

VALIDATION: against the exact interval computation on all 1140 triples <= 20 (including
the 79 known wrapped cases and all gcd cases via reduction), and on all 286 triples of
{1..13} (the THM-601 layer -- expect ZERO with a coincidence).
"""
from fractions import Fraction as F
from math import gcd
import itertools

R = F(1, 14)

def danger(v):
    rv = R / v
    out = []
    for a in range(v + 1):
        lo, hi = F(a, v) - rv, F(a, v) + rv
        lo, hi = max(lo, F(0)), min(hi, F(1))
        if hi > lo: out.append((lo, hi))
    return out

def merge(ivs):
    ivs = sorted(ivs); out = []
    for lo, hi in ivs:
        if out and lo <= out[-1][1]: out[-1] = (out[-1][0], max(out[-1][1], hi))
        else: out.append((lo, hi))
    return out

def inter(A, B):
    out = []
    for l1, h1 in A:
        for l2, h2 in B:
            lo, hi = max(l1, l2), min(h1, h2)
            if hi > lo: out.append((lo, hi))
    return out

DG = {v: merge(danger(v)) for v in range(1, 21)}

def triple_measure(u, v, w):
    return sum(h - l for l, h in inter(inter(DG[u], DG[v]), DG[w]))

def helly_coincidence(u, v, w):
    """Integer criterion: exists non-origin (a,b,c) with the three pairwise conditions."""
    for a in range(u + 1):
        for b in range(v + 1):
            if 14 * abs(a * v - b * u) > u + v - 1: continue
            for c in range(w + 1):
                if (a, b, c) in ((0, 0, 0), (u, v, w)): continue
                if 14 * abs(b * w - c * v) > v + w - 1: continue
                if 14 * abs(a * w - c * u) > u + w - 1: continue
                return (a, b, c)
    return None

print("=" * 90)
print("VALIDATION: Helly integer criterion vs exact interval measure (r = 1/14)")
print("=" * 90)
mismatch = 0
for u, v, w in itertools.combinations(range(1, 21), 3):
    g = gcd(gcd(u, v), w)
    ur, vr, wr = u // g, v // g, w // g
    coin = helly_coincidence(ur, vr, wr)
    clean_measure = (triple_measure(u, v, w) == 2 * R / wr)
    # criterion predicts clean iff NO coincidence
    if (coin is None) != clean_measure:
        mismatch += 1
        if mismatch <= 6:
            print(f"  MISMATCH {(u,v,w)} reduced {(ur,vr,wr)}: coincidence={coin}, "
                  f"measure-clean={clean_measure}, measure={triple_measure(u,v,w)}")
print(f"  triples <= 20: mismatches = {mismatch}/1140")

n13 = 0
for u, v, w in itertools.combinations(range(1, 14), 3):
    g = gcd(gcd(u, v), w)
    if helly_coincidence(u // g, v // g, w // g) is not None:
        n13 += 1
        print(f"  {1,13}-UNIVERSE COINCIDENCE at {(u,v,w)} -- THM-601 would fail!")
print(f"  {{1..13}} universe: coincidences = {n13}/286  (THM-601's Lean layer: expect 0)")
print("\n  => the Lean statement: FORALL 286 gcd-reduced triples of {1..13}: no non-origin")
print("     integer solution -- decidable, no measure theory; monotonicity lifts to all Q.")
print("DONE.")
