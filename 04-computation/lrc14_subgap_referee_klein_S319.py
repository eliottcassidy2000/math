#!/usr/bin/env python3
"""lrc14_subgap_referee_klein_S319.py -- klein-2026-07-19-S319, THM-1290.

Independent exact-M referee for the sub-gap census
(lrc14_subgap_census_klein_S319.c).

Exact M(V) via the complete breakpoint method: the lower envelope
f(t) = min_v ||v t|| is piecewise linear with breakpoints contained in
  { c/(vi+vj) } u { c/(vi-vj) } u { (2c+1)/(2v) } u { c/v },
so max f = max over these candidates.  All arithmetic is exact (Fraction
compare via integers).

Modes:
  --gates          run the six known-family gates (must be EXACT)
  --check FILE     read census output, exact-M every HARD line, classify:
                   BELOW-FLOOR (< 1/14, LRC(14) counterexample!),
                   FLOOR (= 1/14), INGAP (in the open interval -- headline!),
                   EDGE (= HI), ABOVE (> HI)
  --hi N D         the interval's right edge (default 3 41)
"""
import sys
from fractions import Fraction
from math import gcd

def exact_M(V):
    """Exact max_t min_v ||v t|| over t in [0,1), as a Fraction."""
    cands = set()
    n = len(V)
    for i in range(n):
        v = V[i]
        for c in range(1, v):
            cands.add(Fraction(c, v))
        for c in range(2 * v):
            if (2 * c + 1) < 2 * v:
                cands.add(Fraction(2 * c + 1, 2 * v))
        for j in range(i + 1, n):
            s = V[i] + V[j]
            for c in range(1, s):
                cands.add(Fraction(c, s))
            d = abs(V[i] - V[j])
            if d:
                for c in range(1, d):
                    cands.add(Fraction(c, d))
    best = Fraction(0)
    for t in cands:
        m = min(dist_frac(v, t) for v in V)
        if m > best:
            best = m
    return best

def dist_frac(v, t):
    x = v * t
    fl = x - x.numerator // x.denominator  # fractional part
    return min(fl, 1 - fl)

GATES = [
    ([1,2,3,4,5,6,7,8,9,10,11,12,13], Fraction(1,14)),
    ([1,2,3,4,5,6,7,8,9,10,11,12,14], Fraction(1,13)),
    ([1,2,3,4,5,6,7,8,9,10,11,12,26], Fraction(2,27)),
    ([1,2,3,4,5,6,7,8,9,10,11,13,24], Fraction(1,14)),
    ([1,2,3,4,5,6,7,8,9,10,11,13,36], Fraction(3,41)),
    ([1,2,3,5,7,8,9,10,11,12,17,19,104], Fraction(8,105)),
]

def run_gates():
    ok = True
    for V, expect in GATES:
        got = exact_M(V)
        status = "PASS" if got == expect else "FAIL"
        if got != expect:
            ok = False
        print(f"GATE {status}: M({V}) = {got} (expect {expect})")
    print("ALL GATES PASS" if ok else "GATE FAILURE -- DO NOT TRUST THE CENSUS")
    return ok

def classify(m, hi):
    lo = Fraction(1, 14)
    if m < lo:  return "BELOW-FLOOR *** LRC(14) COUNTEREXAMPLE ***"
    if m == lo: return "FLOOR (= 1/14)"
    if m < hi:  return "INGAP *** POPULATES THE INTERVAL ***"
    if m == hi: return "EDGE (= HI)"
    return "ABOVE"

def main():
    args = sys.argv[1:]
    hi = Fraction(3, 41)
    if "--hi" in args:
        i = args.index("--hi")
        hi = Fraction(int(args[i+1]), int(args[i+2]))
        del args[i:i+3]
    if "--gates" in args:
        sys.exit(0 if run_gates() else 1)
    if "--check" in args:
        i = args.index("--check")
        fn = args[i+1]
        n_hard = n_ingap = n_floor = n_below = 0
        for line in open(fn):
            if not line.startswith("HARD"):
                continue
            V = sorted(int(x) for x in line.split()[1:])
            assert len(V) == 13 and len(set(V)) == 13, f"bad line: {line}"
            n_hard += 1
            m = exact_M(V)
            c = classify(m, hi)
            print(f"REFEREE: {V} M = {m} -> {c}")
            if "INGAP" in c: n_ingap += 1
            if "FLOOR" in c and "BELOW" not in c: n_floor += 1
            if "BELOW" in c: n_below += 1
        print(f"REFEREE SUMMARY: hard={n_hard} ingap={n_ingap} floor={n_floor} below-floor={n_below}")
        sys.exit(2 if (n_ingap or n_below) else 0)
    print(__doc__)

if __name__ == "__main__":
    main()
