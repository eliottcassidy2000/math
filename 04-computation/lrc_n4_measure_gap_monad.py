#!/usr/bin/env python3
"""
lrc_n4_measure_gap_monad.py    monad-researcher-2026-06-01-S1

EXACT (Fraction) study of the n=4 MEASURE GAP, speeds <= 24:
  Conjecture: every primitive triple != {1,2,3} has |SAFE| >= 1/28, eq at {1,6,7}.
  Confirm 1/28 exactly at (1,6,7); list the smallest exact |SAFE| values to expose
  the gap above the unique tight value 0 (= {1,2,3}).
"""
from math import gcd
from fractions import Fraction as F
from itertools import combinations

THR = F(1, 4)

def fnorm(x):
    r = x - (x.numerator // x.denominator)
    return min(r, 1 - r)

def exact_measure(speeds):
    bps = set([F(0), F(1)])
    for s in speeds:
        for num in range(1, 4 * s):
            if num % 4 in (1, 3):
                bps.add(F(num, 4 * s))
    pts = sorted(bps)
    measure = F(0)
    for lo, hi in zip(pts, pts[1:]):
        mid = (lo + hi) / 2
        if min(fnorm(F(s) * mid) for s in speeds) >= THR:
            measure += (hi - lo)
    return measure

def main():
    print("n=4 exact measure gap study (monad-S1), speeds<=24\n")
    vals = []
    for a, b, c in combinations(range(1, 25), 3):
        if gcd(gcd(a, b), c) != 1:
            continue
        m = exact_measure((a, b, c))
        vals.append((m, (a, b, c)))
    vals.sort(key=lambda x: x[0])
    print("smallest |SAFE| values (exact):")
    for m, trip in vals[:14]:
        print(f"   {trip}: |SAFE| = {m} = {float(m):.6f}")
    # the unique zero
    zeros = [t for m, t in vals if m == 0]
    print(f"\n  measure-zero triples (speeds<=24): {zeros}")
    smallest_pos = next(m for m, t in vals if m > 0)
    print(f"  smallest POSITIVE measure: {smallest_pos} = {float(smallest_pos):.6f}")
    print(f"  is it 1/28? {smallest_pos == F(1,28)}")
    # confirm (1,6,7)
    print(f"\n  (1,6,7) exact measure = {exact_measure((1,6,7))} (== 1/28: {exact_measure((1,6,7))==F(1,28)})")

if __name__ == "__main__":
    main()
