#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD 3 -- ANGLE 4 (transference) + ANGLE 7b (ORIGINAL: Beurling-Selberg majorant).

ANGLE 7b: the cover indicator 1[all 6 residues hit] is a Boolean combination of the per-speed
sector indicators. By incl-excl, 1[miss residue r] = product over speeds of (1 - 1[e hits r]).
The Selberg/Beurling majorant replaces each sharp sector-indicator by a band-limited trig
polynomial that DOMINATES it; integrating the majorant gives a PROVABLE upper bound on measS7
(or on the uncovered measure from below). The degree-N Selberg majorant has integral
= true + O(1/N). If we pick N so the O(1/N) error < margin, we get a finite-Fourier UPPER bound
certificate. PROBE: does the simplest band-limited majorant (the FEJER / first-moment majorant)
give an upper bound on measS7 that is below cap with room to spare?

ANGLE 4 (transference): replace the "hard" speed set by a "dense model" -- a random/typical speed
set with the same residue statistics -- where coverage is computable. The transference principle
(Green-Tao style) says the hard config's coverage is close to the dense model's up to a small
error governed by the Gowers/box norm of the difference. PROBE: how close is the TRUE measS7 to
the DENSE-MODEL (iid-residue) prediction = the decorrelated p0_decorr? (this is exactly the
resonance error -- so transference error = resonance error, already bounded). Confirm the framing.
"""
import sys, random
from fractions import Fraction as F
from functools import reduce
from math import gcd
sys.stdout.reconfigure(line_buffering=True)

def measS7(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        s = set()
        for e in E:
            v = F(e) * mid; v = v - F(v.numerator // v.denominator)
            s.add((v.numerator * 7) // v.denominator)
        if len(s & set(range(1, 7))) == 6: tot += hi - lo
    return tot

def first_moment_majorant(E):
    """Union-bound (first-moment) UPPER bound on UNCOVERED measure, hence LOWER bound on cover.
       Actually compute: measS7 >= 1 - sum_r meas(miss r). This is a LOWER bound on coverage =>
       UPPER bound on the lonely (uncovered). For the LRC we want UPPER bound on measS7. Use the
       trivial UPPER bound measS7 <= min_r meas(hit r) = min_r (1 - meas(miss r))."""
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    hit = [F(0)] * 7
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        s = set()
        for e in E:
            v = F(e) * mid; v = v - F(v.numerator // v.denominator)
            s.add((v.numerator * 7) // v.denominator)
        for r in range(1, 7):
            if r in s: hit[r] += hi - lo
    return min(hit[r] for r in range(1, 7))  # measS7 <= meas(hit each sector r) for every r

def rand_full(k, hi, seed, maxtry=400):
    random.seed(seed)
    for _ in range(maxtry):
        E = sorted(set([0] + random.sample(range(1, hi), k - 1)))
        if len(E) == k and reduce(gcd, E) == 1:
            res = set(e % 7 for e in E if e != 0)
            if res >= set(range(1, 7)): return E
    return None

def main():
    print("=== ANGLE 7b: simplest UPPER-bound majorant measS7 <= min_r meas(hit r) ===\n")
    caps = {8: 2243/5880, 9: 1979/4004, 10: 55/91, 11: 66/91, 12: 6/7}
    print("k | consec measS7 | min_r meas(hit r) UPPER bound | cap | bound<cap?")
    for k in (8, 9, 10):
        consec = list(range(k))
        m = float(measS7(consec)); ub = float(first_moment_majorant(consec))
        cap = caps[k]
        print(f"{k} | {m:.4f} | {ub:.4f} | {cap:.4f} | {'YES' if ub < cap else 'NO -- bound too weak'}")
    print("  (the trivial single-sector majorant: is it already below cap? If yes, near-proof for consec.)\n")

    print("Worst-case over configs: max_E min_r meas(hit r) vs cap (does the majorant CAP everything?)")
    for k in (8, 9):
        cap = caps[k]; worst = 0; over = 0; n = 0; we = None
        for s in range(200):
            E = rand_full(k, 14, 40000 + s)
            if E is None: continue
            n += 1; ub = float(first_moment_majorant(E))
            if ub > worst: worst = ub; we = E
            if ub > cap: over += 1
        print(f"  k={k}: max_E (min_r meas hit_r) = {worst:.4f} vs cap={cap:.4f} | {over}/{n} configs exceed cap")
        if we: print(f"     worst E={we}")
    print("  => if max single-sector majorant <= cap, then measS7 <= cap for ALL (1-line proof!).")
    print("     if some config's single-sector hit-meas > cap, the majorant is too weak (expected).\n")

    print("=== ANGLE 4: transference error = resonance error (framing confirmation) ===")
    print("  The dense model (iid residues) = p0_decorr; transference error = |measS7 - dense| = resonance err.")
    print("  Already measured <=0.006 at consec base (poserror_consec). So transference IS the wide route,")
    print("  not a new mechanism -- it RE-NAMES the resonance error. The Gowers-norm bound on the error")
    print("  would be the rigorous version of the 32x-lossy resonance bound.")

if __name__ == "__main__":
    main()
