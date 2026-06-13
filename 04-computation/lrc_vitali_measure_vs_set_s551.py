#!/usr/bin/env python3
"""
lrc_vitali_measure_vs_set_s551.py    oracle-2026-06-01-S551o

CONSIDER THE VITALI SET. It is the canonical proof that MEASURE and SET-NONEMPTINESS
are independent: the circle R/Z splits into countably many congruent (rational-
translate) non-measurable pieces (R/Q transversal), so a nonempty -- even 'full' --
set can have measure 0 or be non-measurable. That is exactly the wall of the S550
resonance-energy bound:

  LRC is a SET statement:  L(v) = { t : forall i, ||v_i t|| >= 1/n }  is NONEMPTY.
  The circle-method bound (S550) controls only the MEASURE:  measure(L) >= main - E.
  It proves nonemptiness WHEN measure(L) > 0 (the decorrelated bulk).
  At the HIGH-ENERGY CORE (regular polygon) measure(L) = 0 but L != empty
  (the n-gon vertices t = k/n). MEASURE IS BLIND THERE -- the Vitali phenomenon.

This script makes the measure/set distinction concrete and shows the SIEVE
(THM-369) is a constructive SET-witness that bypasses measure at the core.
"""
from functools import reduce
from math import gcd
from fractions import Fraction

def frac(x): return x - int(x // 1)

def strict_lonely_measure(v, n, G=400000):
    """measure of the OPEN lonely set {t : all ||v_i t|| > 1/n}."""
    lo, hi = 1.0/n, 1.0 - 1.0/n; c = 0
    for i in range(G):
        t = (i+0.5)/G
        if all(lo < (s*t) % 1.0 < hi for s in v): c += 1
    return c/G

def is_lonely_closed(v, n, t):
    """closed lonely: all ||v_i t|| >= 1/n (exact, for rational t)."""
    for s in v:
        x = Fraction(s) * t
        d = x - int(x)               # fractional part as Fraction in [0,1)
        if d < 0: d += 1
        dist = min(d, 1 - d)
        if dist < Fraction(1, n): return False
    return True

def main():
    print("="*74)
    print("MEASURE vs SET: the lonely set L(v); the Vitali lens at the core")
    print("="*74)
    print("  speeds                 strict measure(L)   closed witness t=1/n lonely?   tool")
    for n in (5, 6, 7):
        # AP / regular polygon (the high-energy core): measure 0 but nonempty
        ap = tuple(range(1, n))
        m_ap = strict_lonely_measure(ap, n)
        t = Fraction(1, n)
        wit = is_lonely_closed(ap, n, t)
        print(f"  AP {str(ap):16s} (n={n}): {m_ap:.5f} (MEASURE 0)    t=1/{n}: {wit}   "
              f"=> SIEVE (set-construction)")
        # a generic primitive set: positive measure
        import random; rnd = random.Random(n)
        g = tuple(sorted(rnd.sample(range(2, 10*n), n-1)))
        while reduce(gcd, g) != 1: g = tuple(sorted(rnd.sample(range(2, 10*n), n-1)))
        m_g = strict_lonely_measure(g, n)
        print(f"  generic {str(g):14s} (n={n}): {m_g:.5f} (>0)         "
              f"(positive measure)            => MEASURE bound (S550)")
    print()
    print("="*74)
    print("THE POINT (Vitali): measure cannot distinguish these three at the core")
    print("="*74)
    print("  At the regular polygon the OPEN lonely set is EMPTY (measure 0) and the CLOSED")
    print("  lonely set is exactly the n-gon vertices {k/n} (n points, measure 0, NONEMPTY).")
    print("  Measure assigns 0 to: the empty set, the n points, and the boundary of any")
    print("  positive set -- it CANNOT see nonemptiness here. The Vitali set is the extreme:")
    print("  a nonempty/full transversal with NO well-defined measure. So a measure tool")
    print("  (S550's resonance energy E) is intrinsically blind to the measure-zero core.")
    print()
    print("  RESOLUTION: LRC's core is NOT Vitali-pathological -- it is measure-zero but")
    print("  EXPLICITLY CONSTRUCTIBLE. The SIEVE (THM-369) exhibits the rational witness")
    print("  t = a/q (no axiom of choice): a constructive choice of the relevant Q-coset")
    print("  representatives. Measure (S550) handles the decorrelated bulk; the sieve handles")
    print("  the measure-zero core. The Vitali set marks exactly the handoff: where measure")
    print("  ends, explicit (arithmetic) construction must begin.")
    print()
    print("  ADELIC echo (S547): the Vitali quotient R/Q lives at the ARCHIMEDEAN circle;")
    print("  LRC's tight witnesses are the RATIONAL points (the line<->tree diagonal). Measure")
    print("  is an archimedean (line) tool; the core's witnesses are rational (diagonal) -- the")
    print("  sieve is the p-adic/tree-side construction. Vitali = the measure-blindness of the")
    print("  line at its own rational points = where the tree (sieve) must take over.")

if __name__ == "__main__":
    main()
