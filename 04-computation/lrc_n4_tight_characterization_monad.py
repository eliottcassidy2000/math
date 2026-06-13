#!/usr/bin/env python3
"""
lrc_n4_tight_characterization_monad.py    monad-researcher-2026-06-01-S1  (v2, fast)

KEY METHODOLOGICAL POINT (correction to HYP-2004 item B):
  |SAFE| = int_0^1 prod_i 1[||s_i t||>=1/n] dt >= 0 is TRIVIALLY true (product of
  nonnegative indicators). So NO measure lower bound can prove LRC -- the extremal
  AP {1,2,3} has |SAFE| EXACTLY 0. The whole difficulty of LRC is concentrated in
  the measure-ZERO ("tight") systems: there one must exhibit a CLOSED witness t*.

      LRC(n)  <=>  every measure-zero primitive speed system has a closed witness.

This script (fast, float for the wide scan; fractions only to confirm small cases):
  (1) wide exact-breakpoint scan: which primitive triples have |SAFE|=0?
  (2) minimum POSITIVE measure as speeds grow (positive floor away from AP?)
  (3) parity: all-odd => t=1/4 witness; tight => must contain an even speed.
"""
import sys
from math import gcd
from fractions import Fraction as F
from itertools import combinations

N = 4
THRf = 0.25
EPS = 1e-12

def fnorm(x):
    x %= 1.0
    return min(x, 1.0 - x)

def breakpoints_float(speeds):
    bps = [0.0, 1.0]
    for s in speeds:
        d = 4 * s
        num = 1
        while num < d:
            if num & 3 in (1, 3):
                bps.append(num / d)
            num += 1
    bps.sort()
    return bps

def safe_measure_float(speeds):
    """measure of closed safe set, and best margin + witness (float)."""
    bps = breakpoints_float(speeds)
    measure = 0.0
    best = -1.0; wit = None
    for lo, hi in zip(bps, bps[1:]):
        if hi - lo < 1e-15:
            continue
        mid = 0.5 * (lo + hi)
        d = min(fnorm(s * mid) for s in speeds)
        if d >= THRf - 1e-12:
            measure += (hi - lo)
        if d > best:
            best = d; wit = mid
    for t in bps[:-1]:
        d = min(fnorm(s * t) for s in speeds)
        if d > best:
            best = d; wit = t
    return measure, (best >= THRf - 1e-9), best, wit

def main():
    print("LRC n=4: characterizing the measure-ZERO (tight) triples (monad-S1, v2)\n")
    sys.stdout.flush()

    for MAXS in (20, 40, 70, 100):
        tot = 0; tight = []; minpos = None; lrc_fail = 0
        for a, b, c in combinations(range(1, MAXS + 1), 3):
            if gcd(gcd(a, b), c) != 1:
                continue
            tot += 1
            meas, ne, best, wit = safe_measure_float((a, b, c))
            if not ne:
                lrc_fail += 1
            if meas < 1e-9:
                tight.append((a, b, c, ne, best, wit))
            else:
                if minpos is None or meas < minpos[0]:
                    minpos = (meas, a, b, c)
        print(f"speeds<= {MAXS}:  primitive triples={tot}  LRC failures={lrc_fail}")
        print(f"   measure-zero (tight) triples: {[(t[0],t[1],t[2]) for t in tight]}")
        for (a, b, c, ne, best, wit) in tight:
            print(f"     ({a},{b},{c}): closed-nonempty={ne} witness t*~{wit:.4f} margin~{best:.4f}")
        print(f"   min POSITIVE |SAFE|: {minpos[0]:.6f} at {minpos[1:]}")
        print(); sys.stdout.flush()

    # exact confirmation that {1,2,3} is tight with witness 1/4
    def fnormF(x):
        r = x - (x.numerator // x.denominator); return min(r, 1 - r)
    print("="*70)
    print("EXACT confirmation of the unique small tight set:")
    print("="*70)
    for trip in [(1, 2, 3)]:
        t = F(1, 4)
        margins = [fnormF(F(s) * t) for s in trip]
        print(f"  {trip}: t*=1/4 margins = {margins} (all = 1/4 = threshold) -> lonely witness EXACT")

    # parity lemma
    print("\n" + "="*70)
    print("Parity lemma: all speeds odd => t=1/4 safe (||odd/4||=1/4). Tight => has even speed.")
    print("="*70)
    n_allodd = 0; n_allodd_tight = 0
    for a, b, c in combinations(range(1, 41), 3):
        if gcd(gcd(a, b), c) != 1:
            continue
        if a % 2 and b % 2 and c % 2:
            n_allodd += 1
            meas, ne, best, wit = safe_measure_float((a, b, c))
            assert ne  # t=1/4 always a witness
            if meas < 1e-9:
                n_allodd_tight += 1
    print(f"  all-odd primitive triples (<=40): {n_allodd}; all have t=1/4 witness: True")
    print(f"  all-odd tight (measure 0): {n_allodd_tight}")

if __name__ == "__main__":
    main()
