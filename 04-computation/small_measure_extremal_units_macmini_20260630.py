#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The small-measure regime, anchored: the consecutive extremal W={1..n-1} has lonely measure EXACTLY 0
at every n, its lonely set = the phi(n) UNITS {k/n: gcd(k,n)=1} (touch-points), and the union bound goes
NEGATIVE so measure-from-below is structurally doomed. mac-mini-2026-06-30-S38 (HYP-3607).
Generalizes klein-S8's n=14 {1..13}->units fact to all n.
"""
from __future__ import annotations
import functools, math
from fractions import Fraction as F
print = functools.partial(print, flush=True)


def phi(n):
    return sum(1 for a in range(1, n + 1) if math.gcd(a, n) == 1)


def norm_frac(num, den):
    f = (num % den)
    return min(f, den - f)


def main():
    print("=" * 78)
    print("THE SMALL-MEASURE REGIME -- the extremal is measure-ZERO, lonely set = the units (S38)")
    print("=" * 78)

    print("\n[1] union bound for the lonely measure goes NEGATIVE (measure can vanish):")
    print(f"    {'n':>3} {'#danger combs':>13} {'each meas':>10} {'lonely bound (2-n)/n':>21}")
    for n in (3, 5, 7, 14):
        print(f"    {n:>3} {n-1:>13} {f'2/{n}':>10} {f'{(2-n)/n:+.3f}':>21}")
    print("    => for n>=3 the n-1 danger combs (each 2/n) can cover all of [0,1); at the extremal they DO")
    print("       (up to measure-0 touch-points). No measure-from-below argument can work at the extremal.")

    print("\n[2] consecutive extremal W={1..n-1}: lonely set = touch-points t=k/n with min_j||jk/n||>=1/n:")
    print(f"    {'n':>3} {'lonely touch-points (k with gcd(k,n)=1)':>42} {'count':>6} {'phi(n)':>7}")
    for n in range(3, 10):
        pts = [k for k in range(1, n)
               if all(norm_frac(w * k, n) * n >= 1 for w in range(1, n))]
        ok = (pts == [k for k in range(1, n) if math.gcd(k, n) == 1])
        print(f"    {n:>3} {str([f'{k}/{n}' for k in pts]):>42} {len(pts):>6} {phi(n):>7}"
              f"{'  [= units mod n]' if ok else '  [MISMATCH]'}")
    print("    PROOF: t=k/n lonely <=> {jk mod n} all nonzero <=> gcd(k,n)=1 (else j=n/gcd puts a runner at")
    print("    0). So the lonely set is EXACTLY the phi(n) units, measure 0, in phi(n)/2 antipodal pairs.")

    print("\n" + "=" * 78)
    print("ANCHOR: the LRC extremal is measure-ZERO with lonely set = the phi(n) units (n=14: (Z/14)* =")
    print("{1,3,5,9,11,13}, klein-S8). The conjecture holds there by EXISTENCE/COUNTING (the units / the odd")
    print("cycle / non-bipartiteness certificate THM-590), NOT by measure. The small-measure regime is the")
    print("heart: measure vanishes, counting carries it.")
    print("=" * 78)


if __name__ == "__main__":
    main()
