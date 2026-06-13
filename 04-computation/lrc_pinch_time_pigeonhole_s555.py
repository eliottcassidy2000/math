#!/usr/bin/env python3
"""
lrc_pinch_time_pigeonhole_s555.py    oracle-2026-06-01-S555o

Attempting the user's PINCH-TIME PIGEONHOLE for LRC@14, made precise.

PINCH TIME of a pair (a,b): t = 1/(a+b). Then a sits at a/(a+b), b at b/(a+b),
symmetric about 1/2, both at circular distance min(a,b)/(a+b) >= 1/(a+b) from the
observer. So if a+b <= n, the PINCH CLEARS THE PAIR (>= 1/n). [verified below]

A THIRD runner w at the pinch t=1/(a+b): ||w/(a+b)|| >= 1/(a+b) >= 1/n iff (a+b) does
NOT divide w. So:
   the pinch t=1/(a+b) is LONELY  <=>  no runner is divisible by (a+b).
That is EXACTLY the denominator sieve (THM-369) at q=a+b. The pinch-pigeonhole over
the candidate sums s (pairs summing to s<=n) is: 'some s<=n has no multiple in the
set' = NOT sieve-covered. A counterexample is SIEVE-COVERED (condition A1, S554): it
has a multiple of EVERY q in {2..n}, so EVERY rational pinch t=p/q (q<=n) is spoiled.

So: the pinch pigeonhole = the sieve; it works for the (decorrelated) majority but is
defeated by a sieve-covered counterexample BY CONSTRUCTION. We verify this, and show
the actual lonely time of a sieve-covered set has DENOMINATOR > n (the fine regime),
out of reach of any pinch at a rational t=p/q with q<=n.
"""
from itertools import combinations
from functools import reduce
from math import gcd
from fractions import Fraction

N = 14

def dist0_frac(p):
    """circular distance of a Fraction position to 0."""
    f = p - (p.numerator // p.denominator)
    if f < 0: f += 1
    return min(f, 1 - f)

def pinch_clears_pair(a, b):
    s = a + b; t = Fraction(1, s)
    da = dist0_frac(Fraction(a, s)); db = dist0_frac(Fraction(b, s))
    return da >= Fraction(1, N) and db >= Fraction(1, N), da, db

def lonely_at_rational(v, p, q):
    t = Fraction(p, q)
    return all(dist0_frac(Fraction(s) * t) >= Fraction(1, N) for s in v)

def any_rational_pinch_lonely(v):
    """is ANY t=p/q (q in 2..N, reduced) lonely? (= the rational/sieve regime)"""
    for q in range(2, N+1):
        for p in range(1, q):
            if gcd(p, q) == 1 and lonely_at_rational(v, p, q):
                return True, Fraction(p, q)
    return False, None

def lonely_float(v, G=400000):
    for i in range(G):
        t = (i+0.5)/G
        if all(1.0/N < (s*t) % 1.0 < 1.0 - 1.0/N for s in v): return True, t
    return False, None

def main():
    print("="*74)
    print("(1) PINCH CLEARS THE PAIR: pairs summing to 14, t=1/14, distance min(a,b)/14")
    print("="*74)
    for a in range(1, 7):
        b = 14 - a
        ok, da, db = pinch_clears_pair(a, b)
        print(f"  pair ({a},{b}) sum 14, pinch t=1/14: dist {a}={da}, dist {b}={db}  cleared={ok}")
    print("  small-sum pairs (a+b<=14) cleared with margin; (1,13) is borderline (=1/14).")
    print()
    print("="*74)
    print("(2) THE PINCH t=1/(a+b) IS THE SIEVE: lonely iff no runner divisible by (a+b)")
    print("="*74)
    print("  so the pinch-pigeonhole 'some pinch clears all' = 'some s<=14 has no multiple'")
    print("  = NOT sieve-covered. A counterexample is SIEVE-COVERED (A1) -> defeats it.")
    print()
    print("="*74)
    print("(3) VERIFY: sieve-covered sets have NO lonely rational t=p/q (q<=14), yet ARE")
    print("    lonely at a FINE (denominator > 14) time -- out of reach of the pinch.")
    print("="*74)
    import random
    rnd = random.Random(1); tested = 0; both = 0
    def sieve_covered(v): return all(any(s % q == 0 for s in v) for q in range(2, N+1))
    examples = []
    while tested < 40:
        v = tuple(sorted(rnd.sample(range(1, 70), 13)))
        if reduce(gcd, v) != 1 or not sieve_covered(v): continue
        tested += 1
        rok, rt = any_rational_pinch_lonely(v)
        fok, ft = lonely_float(v)
        if (not rok) and fok: both += 1
        if len(examples) < 5: examples.append((rok, fok, round(ft,5) if ft else None))
    print(f"  of {tested} sieve-covered primitive 13-sets: NO rational pinch lonely AND lonely-fine"
          f" for {both}/{tested}.")
    print(f"  examples (rational-pinch-lonely?, fine-lonely?, fine t): {examples}")
    print()
    print("="*74)
    print("VERDICT")
    print("="*74)
    print("  The pinch-time pigeonhole is the RIGHT SHAPE, but every rational pinch t=p/q")
    print("  (q<=n) is a denominator-sieve witness (THM-369), lonely iff no runner divisible")
    print("  by q. A counterexample is SIEVE-COVERED (A1) -- it has a multiple of every q<=n")
    print("  -- so it defeats EVERY rational pinch BY CONSTRUCTION; a single multiple of n")
    print("  spoils every n-gon vertex t=j/n (that runner sits at 0). The pinch pigeonhole")
    print("  thus = the sieve; it settles the decorrelated majority but NOT the sieve-covered")
    print("  core. The genuine lonely time of a counterexample-candidate lives at a FINE")
    print("  denominator > n (S18 fine regime) -- where the pinch/sieve cannot reach. To")
    print("  beat the core one needs irrational/fine pinches (the measure core S550, or the")
    print("  windows S552), which carry the unresolved multiples-of-n* coupling.")

if __name__ == "__main__":
    main()
