#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S97 -- HYP-3853: (A) THM-597 end-to-end collapse-rate law
verification: the general argmax formula c(S) = (1/n) * sum_{t* in argmax}
(1/v_plus + 1/v_minus) against the direct last-cell slope, on the full tight
locus q = 5, 6, 8, 14 (incl. the q=5 perm-lift beater and the q=8 relift
beater); confirms every witness denominator is divisible by q (THM-597(a))
and that all these sets are clean (argmax = unit fractions).
(B) hpartA LOCAL-COVERING PROBE: can j distinct-spacing arc-combs (far
runners) cover a short window?  Grid probe over comb phases; reports the max
covered fraction and the residual (the windowed-MN deficiency to be certified
by the BBMST-style distortion LP).
"""
from fractions import Fraction as F
from math import gcd
import itertools, random, sys
sys.path.insert(0, '04-computation')
from lonely_profile import profile

def dist(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def argmax_set(S, q):
    """All binding-grid rationals with min_v ||v t|| == 1/q exactly."""
    Sl = sorted(set(S)); dens = set()
    for i, v in enumerate(Sl):
        for w in Sl[i:]:
            dens.add(v + w)
            if w > v: dens.add(w - v)
        dens.add(2 * v)
    A = set()
    for den in dens:
        for m in range(1, den):
            t = F(m, den)
            if min(dist(v * t) for v in Sl) == F(1, q):
                A.add(t)
    return sorted(A)

def binding_sides(S, t, q):
    """v_plus (fastest binding constraint decreasing to the right of t) and
    v_minus (fastest binding on the left): binding v has ||v t|| = 1/q; the
    side is determined by whether v t mod 1 is at +1/q (constraint tightens
    leftwards) or -1/q (tightens rightwards)."""
    vplus = 0; vminus = 0
    for v in sorted(S):
        x = (v * t) - int(v * t)
        if x == F(1, q):          # fractional part +1/q: moving left decreases clearance
            vminus = max(vminus, v)
        elif x == 1 - F(1, q):    # at -1/q: moving right decreases clearance
            vplus = max(vplus, v)
    return vplus, vminus

print("=" * 78)
print("(A) THM-597: c(S) = (1/q) sum_{argmax} (1/v+ + 1/v-)  vs direct slope")
print("=" * 78)
tight_locus = {
    5:  [[1, 2, 3, 4], [1, 3, 4, 7]],
    6:  [[1, 2, 3, 4, 5], [1, 3, 4, 5, 9]],
    8:  [[1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 5, 7, 12], [1, 4, 5, 6, 7, 11, 13]],
    14: [list(range(1, 14)), list(range(1, 12)) + [13, 24]],
}
for q, sets in tight_locus.items():
    units = [a for a in range(1, q) if gcd(a, q) == 1]
    for S in sets:
        A = argmax_set(S, q)
        # (a) every witness denominator divisible by q
        dens_ok = all(t.denominator % q == 0 or q % t.denominator == 0 and False
                      for t in A)
        dens_ok = all((t.denominator % q == 0) if t.denominator >= q else (q % t.denominator == 0)
                      for t in A)
        # actually THM-597(a): q | denominator (in lowest terms den = q*s)
        dens_ok = all(t.denominator % q == 0 for t in A)
        clean = (sorted(A) == sorted(F(a, q) for a in units))
        c_formula = sum(F(1, vp) + F(1, vm) for t in A
                        for vp, vm in [binding_sides(S, t, q)]) / q
        c_direct = profile(S, F(1, q)).critical_slope(q)
        ok = "OK" if c_formula == c_direct else f"MISMATCH direct={c_direct}"
        print(f"  q={q} S={S}")
        print(f"    |argmax|={len(A)} (phi={len(units)}); q|denominators: {dens_ok}; "
              f"clean(argmax=unit fracs): {clean}")
        print(f"    c formula = {c_formula} = {float(c_formula):.6f}  [{ok}]")

print()
print("=" * 78)
print("(B) hpartA probe: max window coverage by j distinct-spacing arc-combs")
print("=" * 78)
# window I = [0, |I|); comb of spacing 1/w with phase ph covers
# {x : dist(w x - ph, Z) < 1/14 } within I.  Probe max covered fraction over phases.
def covered_fraction(ws, phases, width, grid=4000):
    cov = 0
    for i in range(grid):
        x = (i + 0.5) / grid * width
        for w, ph in zip(ws, phases):
            y = w * x - ph
            y -= round(y)
            if abs(y) < 1 / 14:
                cov += 1
                break
    return cov / grid

random.seed(97)
for width in (0.02, 0.01):
    for j, ws in [(6, [101, 137, 173, 211, 251, 293]),
                  (7, [101, 137, 173, 211, 251, 293, 331]),
                  (7, [100, 200, 400, 800, 1600, 3200, 6400])]:
        best = 0.0
        for trial in range(3000):
            phases = [random.random() for _ in ws]
            frac = covered_fraction(ws, phases, width)
            best = max(best, frac)
        A_local = sum(2 / 14 for _ in ws)  # per-comb density 2r = 1/7
        print(f"  |I|={width}, j={j}, spacings 1/{ws[0]}..: mass budget A={A_local:.3f}, "
              f"max covered fraction over 3000 phase draws = {best:.4f} "
              f"{'(FULL COVER reached)' if best >= 0.999 else '(residual > 0)'}")
print("  => the windowed-MN / BBMST-distortion target: certify the residual")
print("     analytically for j <= 13 distinct spacings at density 2r = 1/7 each.")
