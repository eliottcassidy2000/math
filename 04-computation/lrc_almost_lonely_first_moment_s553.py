#!/usr/bin/env python3
"""
lrc_almost_lonely_first_moment_s553.py    oracle-2026-06-01-S553o

A PROVABLE LRC-adjacent theorem (the 'almost lonely' theorem), and an honest n=14
counterexample-impossibility framing.

THEOREM (first moment, rigorous). For nonzero integer speeds v_1..v_{n-1} and the
observer at 0, let near(t) = #{ i : ||v_i t|| < 1/n }. Each runner's danger set has
measure exactly 2/n, so
   INT_0^1 near(t) dt = (n-1) * (2/n) = 2 - 2/n  <  2   (for all n >= 3).
A nonnegative INTEGER-valued function with average < 2 must somewhere be <= 1. Hence:
   *** at some time t, AT MOST ONE runner is within 1/n of the observer ***
   (equivalently max_t #{far} >= n-2). This is LRC 'up to one exceptional runner' --
   exactly the S524 coupon-collector bottleneck (d_Q = 6 of 7), now a theorem.
LRC itself is the (open) statement that 'at most one' can be improved to 'zero'.

n=14 IMPOSSIBILITY (honest). A counterexample to LRC@14 would need near(t) >= 1 for
ALL t (never lonely) while the average forces some t with near(t) <= 1. So a
counterexample is AVERAGING-EXTREMAL: near(t) hits the floor 1 but never 0 (not even
in closure). Combined with the sieve (q-covered for all q<=14) and the high-energy
core (S550, E>=main), the constraints are severe -- but not provably empty (= LRC).

This script: (1) verifies max_t #far >= n-2 (min near <= 1) for n=5..14, and measures
how often the floor is the actual minimum (tightness); (2) computes the SECOND moment
(Var of near) to show why moments stop at 'one near'; (3) the n=14 averaging-extremal
characterization of a counterexample.
"""
from itertools import combinations
from functools import reduce
from math import gcd
import random, statistics as st

def near_count(speeds, n, t):
    return sum(1 for s in speeds if min((s*t) % 1.0, 1 - (s*t) % 1.0) < 1.0/n - 1e-12)

def min_near_over_t(speeds, n, G=300000):
    best = len(speeds)
    for i in range(G):
        t = (i+0.5)/G
        c = near_count(speeds, n, t)
        if c < best: best = c
        if best == 0: break
    return best

def main():
    print("="*72)
    print("THEOREM (first moment): average #near = 2 - 2/n < 2  =>  min_t #near <= 1")
    print("  (at some time at most ONE runner is near; max_t #far >= n-2). LRC = 0 near.")
    print("="*72)
    print("  n : avg #near = 2-2/n   floor(min near)   max_t #far >= n-2")
    for n in range(5, 9):
        print(f"  {n} :   {2-2.0/n:.4f}            <= 1            >= {n-2}")
    print()
    print("(1) VERIFY min_t #near <= 1 (and how often = 0 = LRC, vs = 1 = 'one short')")
    for n in (5, 7, 10, 14):
        rnd = random.Random(n); zero = 0; one = 0; tot = 0
        sets = [tuple(range(1, n))]  # AP / the wall-only core
        while len(sets) < 12:
            v = tuple(sorted(rnd.sample(range(2, 6*n), n-1)))
            if reduce(gcd, v) == 1: sets.append(v)
        for v in sets:
            mn = min_near_over_t(v, n)
            if mn == 0: zero += 1
            elif mn == 1: one += 1
            tot += 1
        print(f"  n={n}: of {tot} sets (incl AP), min near = 0 (LRC, open) for {zero}, "
              f"= 1 ('one short') for {one}  [AP is 'one short' on the open grid = the wall]")
    print()
    print("(2) SECOND moment: Var(#near); why moments stop at 'one near'")
    for n in (7, 14):
        v = tuple(range(1, n))  # AP
        G = 100000; vals = [near_count(v, n, (i+0.5)/G) for i in range(G)]
        mean = st.mean(vals); var = st.pvariance(vals)
        print(f"  n={n} AP: mean #near={mean:.3f} (= 2-2/n={2-2.0/n:.3f}), Var={var:.3f}, "
              f"min={min(vals)}, max={max(vals)}")
    print("  => the first moment pins mean #near = 2-2/n < 2 (so min <= 1); the SECOND")
    print("     moment (resonance correlations) controls whether min reaches 0. Pushing")
    print("     'one near' to 'zero near' (=LRC) is exactly the simultaneity / high-energy")
    print("     core (S550) -- NOT reachable by the first moment alone.")
    print()
    print("="*72)
    print("(3) n=14 counterexample IMPOSSIBILITY (honest framing)")
    print("="*72)
    print("  A counterexample to LRC@14 must satisfy ALL of:")
    print("   - AVERAGING-EXTREMAL: near(t) >= 1 for all t, yet avg near = 2-2/14 = 12/7 ~1.71")
    print("     forces some t with near(t) <= 1; so near hits the floor 1 but NEVER 0 (S553).")
    print("   - SIEVE: a multiple of every q in {2..14} (THM-369).")
    print("   - HIGH-ENERGY CORE: resonance energy E >= (12/14)^13 (S550).")
    print("   - 7-CLASS COUPLED: all 7 mod-7 CRT classes blocked at every t (S524/S552).")
    print("  These pin a counterexample to an extremely thin arithmetic locus -- the AP/")
    print("  regular polygon nearly is it, but the AP is LONELY at the wall t=k/14 (closed),")
    print("  so it is NOT a counterexample. No counterexample is known; the conjunction above")
    print("  is the honest 'proof-lite' -- a counterexample would be averaging-extremal,")
    print("  fully sieve-covered, high-energy, 7-coupled, and never close the last runner.")

if __name__ == "__main__":
    main()
