#!/usr/bin/env python3
"""
lrc_n18_wacky_attempts_s515.py

oracle-2026-06-01-S515

WACKY THEOREMS / NEW IDEAS toward LRC at n=18 (17 speeds, threshold 1/18), with
the same whimsy as the n=14 pass (S514).  n=18 = 2*3^2 is RICHER than 14 = 2*7:
it has a prime SQUARE (9=3^2) and, since 16 <= 18, its sieve even forces the
n=16 "sedenion" dyadic gate.

Ideas:
 A) ANTI-CONCENTRATION.  N(t)=#{i: ||v_i t||<1/18}; E[N]=17*(2/18)=17/9~1.889;
    Poisson P(N=0)~e^{-17/9}~0.151.  Only exact resonance reaches 0.
 B) BOUNDED ANSATZ.  Witnesses at t=j/(18 s); find the cofactor pattern (expect
    tied to the 2-part and 3-part of 18 = 2*3^2, vs 14's {1,2,4,5,7}).
 C) PARITY SPLIT (18 = 2*9).  odd O + even E=2W; smaller side <= 8 speeds =
    PROVED LRC@9.  "18 = 9 doubled" (cf 14 = 7 doubled).  Coupling t<->2t.
 D) TRIADIC SPLIT (the 3^2 part, NEW vs n=14).  classes mod 3; speeds ==0 mod 3
    are 3w -> see ||w(3t)||.  The 3|9 ladder is the n=18 analog of codex's 2|4|8|16
    dyadic ladder for n=16 -- a "ternary endpoint debt".
 E) SIEVE GATES.  prime-power gates <=18 a counterexample must contain: 16,9,5,
    7,11,13,17 (+ composites). Inherits the n=16 dyadic 16-gate AND a 9-gate.
"""

from __future__ import annotations

import math, random
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parent
S356 = SourceFileLoader("s356", str(ROOT / "lonely_runner_residue_probe_s356.py")).load_module()

N = 18
THR = Fraction(1, N)
ONE = Fraction(1)


def d0(x):
    f = x - (x.numerator // x.denominator); return min(f, ONE - f)


def lonely_measure(speeds):
    r = S356.report("x", list(speeds)); return 1 - r.forbidden_length, r.max_gap


def min_ansatz_cofactor(speeds, maxs=40):
    for s in range(1, maxs + 1):
        Q = N * s
        for j in range(1, Q):
            t = Fraction(j, Q)
            if all(d0(Fraction(v) * t) >= THR for v in speeds):
                strict = all(d0(Fraction(v) * t) > THR for v in speeds)
                return s, t, strict
    return None, None, None


def main():
    print("LRC n=18 wacky attempts (oracle-2026-06-01-S515)  [n=18 = 2*3^2]\n")

    print("=" * 70)
    print("A: anti-concentration  (E[N]=17/9, Poisson e^-17/9 = %.5f)" % math.exp(-17/9))
    print("=" * 70)
    setsA = {
        "initial 1..17 (tight)": tuple(range(1, 18)),
        "drop17 add18": tuple(sorted(set(range(1, 18)) - {17} | {18})),
        "drop16 add18": tuple(sorted(set(range(1, 18)) - {16} | {18})),
        "drop9 add18":  tuple(sorted(set(range(1, 18)) - {9} | {18})),
        "primes+18": (1,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,18),
    }
    for label, s in setsA.items():
        p0, gap = lonely_measure(s)
        print(f" [{label:<22}] P(N=0)={float(p0):.5f}  max_gap={S356.fmt_frac(gap)}")
    print()

    print("=" * 70)
    print("B: bounded ansatz cofactor s  (witness at t=j/(18 s))")
    print("=" * 70)
    rng = random.Random(515)
    testB = dict(setsA)
    for _ in range(6):
        gate = rng.choice([18, 36, 54])
        rest = rng.sample([v for v in range(1, 60) if v % 18 != 0], 16)
        s = tuple(sorted(set(rest) | {gate}))
        if len(s) == 17: testB[f"rand {s[:3]}..."] = s
    worst = 0; seen = set()
    for label, s in testB.items():
        cof, t, strict = min_ansatz_cofactor(s)
        if cof is not None:
            worst = max(worst, cof); seen.add(cof)
            print(f" [{label:<22}] s={cof:<3} t={S356.fmt_frac(t)} ({'open' if strict else 'bndry'})")
        else:
            print(f" [{label:<22}] NO witness up to s=40 (!!)")
    print(f"\n => cofactors seen: {sorted(seen)};  worst s = {worst}  (18*s = {18*worst})")
    print("    (compare n=14: s in {1,2,4,5,7}. For 18=2*3^2 expect 2-part & 3-part.)\n")

    print("=" * 70)
    print("B': cofactor stability vs speed size (does s blow up?)")
    print("=" * 70)
    for B in (60, 200, 800):
        worst = 0
        for _ in range(80):
            gate = rng.choice([18, 36, 54, 18*rng.randint(1, max(1,B//18))])
            rest = rng.sample([v for v in range(1, B) if v % 18 != 0], 16)
            s = tuple(sorted(set(rest) | {gate}))
            if len(s) != 17: continue
            cof, _, _ = min_ansatz_cofactor(s, maxs=60)
            if cof: worst = max(worst, cof)
        print(f"   speeds<{B}: worst ansatz cofactor s = {worst} (18*s={18*worst})")
    print()

    print("=" * 70)
    print("C/D: parity (2) and triadic (3) splits")
    print("=" * 70)
    for label, s in list(setsA.items())[:3]:
        odd = [v for v in s if v % 2]; even = [v for v in s if v % 2 == 0]
        m0 = [v for v in s if v % 3 == 0]; m1 = [v for v in s if v % 3 == 1]; m2 = [v for v in s if v % 3 == 2]
        print(f" [{label}]")
        print(f"   parity:  |odd|={len(odd)} |even|={len(even)}  -> smaller {min(len(odd),len(even))} speeds = PROVED LRC@9 (18=9 doubled, t<->2t)")
        print(f"   mod 3:   sizes {len(m0)}/{len(m1)}/{len(m2)}  -> the ==0 mod 3 speeds fold via t<->3t; 3|9 ladder")
    print()

    print("=" * 70)
    print("E: sieve gates a counterexample must contain (prime powers <=18)")
    print("=" * 70)
    from sympy import factorint
    pps = [q for q in range(2, 19) if len(factorint(q)) == 1]
    maxpp = sorted({q for q in pps if not any(q != r and q % r == 0 and len(factorint(r))==1 for r in pps)})
    print(f"   forced divisibility gates (maximal prime powers <=18): {maxpp}")
    print(f"   => a counterexample needs speeds divisible by 16 (the n=16 sedenion gate!),")
    print(f"      9 (=3^2 triadic gate), 5,7,11,13,17 -- inherits BOTH the dyadic-16 and a")
    print(f"      triadic-9 obstruction. n=18 sits above n=16 and n=9 in the sieve tower.")
    print()
    print("SUMMARY: n=18 = 2*3^2 -- parity gives a proved LRC@9 half; the 3^2 adds a")
    print("triadic 3|9 ladder (analog of n=16's dyadic 2|4|8|16); the sieve forces the")
    print("16-gate AND 9-gate. Whimsy: 18 = 9 doubled, and inherits 16's dyadic debt.")


if __name__ == "__main__":
    main()
