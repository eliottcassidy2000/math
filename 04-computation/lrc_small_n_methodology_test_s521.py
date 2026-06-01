#!/usr/bin/env python3
"""
lrc_small_n_methodology_test_s521.py   claudebox-2026-06-01-S521

Honest test: does the multiplicative-walk methodology (base reduction + finite-Q
residue covering) PROVE small LRC cases on its own?

ANSWER: No -- not without an external speed bound.  A constraint search shows
that for ANY fixed finite Q = {2,...,Qmax}, there is a primitive system blocked
at every q in Q, namely

        (1, lcm(2..Qmax))     [and variants],

whose large speed is ≡ 0 at every q ≤ Qmax (each divides the lcm), so it is unsafe
at all those rational times -- yet it is lonely at q ≈ Qmax+1.  Hence the minimal
lonely denominator is UNBOUNDED over all speeds (~ log of the max speed), and no
fixed Q covers every system.

CONSEQUENCE.  The methodology REDUCES LRC(n) to a residue covering, but closing it
needs Q to scale with the speeds -- i.e. a SPEED BOUND on minimal counterexamples
is essential (it is not a removable nuisance).  With such a bound B(n), one takes
Q ~ log B(n) and the covering becomes a finite check; this is how the methodology
re-proves the known small cases (n<=7), not by itself.

This corrects an earlier overclaim ("residue witnesses => no speed bound needed").
The base reduction (non-fully-covered => lonely at 1/q) remains valid and useful.
"""
from math import gcd, lcm
from functools import reduce
from fractions import Fraction as F

def lonely_at(speeds, n, t):
    thr = F(1, n)
    for v in speeds:
        x = (F(v) * t) % 1
        if min(x, 1 - x) < thr: return False
    return True

def min_lonely_denom(speeds, n, qmax):
    for q in range(2, qmax + 1):
        for a in range(1, q):
            if gcd(a, q) == 1 and lonely_at(speeds, n, F(a, q)): return q
    return None

def main():
    print("Does fixed-Q residue covering prove small LRC? Test the (1, lcm(2..Q)) adversary.\n")
    print("n=3 (2 movers): system (1, lcm(2..Qmax)) is blocked at all q<=Qmax, lonely just above:")
    print(f"  {'Qmax':>5} {'lcm(2..Qmax)':>14} {'min lonely q':>13}")
    for Qmax in [8, 10, 12, 15, 20, 25, 30]:
        N = reduce(lcm, range(2, Qmax + 1))
        d = min_lonely_denom([1, N], 3, Qmax + 8)
        print(f"  {Qmax:>5} {N:>14} {str(d):>13}   (> Qmax: {d is not None and d > Qmax})")
    print("\n=> minimal lonely denominator grows with the speed (no fixed Q covers all systems).")
    print("   The methodology needs a SPEED BOUND to bound Q; then the covering is a finite check.")
    print("   (Base reduction 'non-fully-covered => lonely at 1/q' is unaffected and remains proved.)")

if __name__ == "__main__":
    main()
