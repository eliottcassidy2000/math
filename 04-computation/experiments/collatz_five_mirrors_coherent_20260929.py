#!/usr/bin/env python3
"""Mechanism test for the resonant Fourier coefficient of the 3-adic Syracuse law (S21, 2026-09-29).

At the resonant frequency t = 2^s the coefficient mu_hat_h(2^s) = sum_w 2^-A(w) e(2^s Y_h(w)/3^h) is split into
    coh(s) = the sum over the NO-DESCENT words (prefix sums P_j = a_1 + ... + a_j < j log_2 3 for all j <= h)
    rest(s) = the sum over the descending words,
by an exact dynamic programme over (total cost A, level j, prefix sum P) with the phase factors
    e(2^s Y_h / 3^h) = prod_(j=1)^h e((2^(s - A + P_(j-1)) mod 3^j) / 3^j),
and compared with the full coefficient from the closed recursion on the powers of two.  Also: the size of
mu_hat_h at random units (random exponents j_0 mod the order of 2), against the resonant window.
Run: python 04-computation/experiments/collatz_five_mirrors_coherent_20260929.py
"""
from __future__ import annotations

import cmath
import math
import random
import sys
from fractions import Fraction

import numpy as np

LOG23 = math.log2(3)
AMAX = 40


def closed_family(HMAX: int, j_lo: int, j_hi: int):
    """m_n(j) = mu_hat_n(2^j mod 3^n) for j in [j_lo, j_hi] at n = HMAX, by the exact closed recursion."""
    weights = [2.0 ** (-a) for a in range(1, AMAX + 1)]
    prev = {j: 1.0 + 0j for j in range(j_lo - AMAX * HMAX, j_hi + 1)}
    for n in range(1, HMAX + 1):
        mod = 3 ** n
        lo = j_lo - AMAX * (HMAX - n)
        pw = {j: pow(2, j, mod) for j in range(lo - AMAX, j_hi + 1)}
        cur = {}
        for j in range(lo, j_hi + 1):
            s = 0j
            for a in range(1, AMAX + 1):
                s += weights[a - 1] * cmath.exp(2j * math.pi * float(Fraction(pw[j - a], mod))) * prev.get(j - a, 0j)
            cur[j] = s
        prev = cur
    return prev


def coherent_split(h: int, s: int, Awindow: int = 45, c: float = 0.0):
    """Return (coh, coh_abs_sum, mass_ND, rest_ND_cost_mass) where coh = sum over no-descent words of 2^-A e(...),
    coh_abs_sum = sum of |terms| (for the coherent fraction), mass_ND = P_h.  The sum runs over totals A in
    [s - Awindow, s + Awindow] (the no-descent words have A ~ h log_2 3; the window covers all but a negligible mass)."""
    coh = 0j
    coh_abs = 0.0
    mass = 0.0
    lim = [0] + [int(math.floor(j * LOG23 + c - 1e-12)) for j in range(1, h + 1)]  # P_j must be <= lim[j]: excursion above the critical line at most c bits; P_0 = 0
    for A in range(max(h, s - Awindow), s + Awindow + 1):
        # dp over prefix sums with complex weights; states P in [0, lim[j]]
        dp = np.zeros(lim[h] + 2, dtype=complex)
        dp[0] = 1.0
        dpabs = np.zeros(lim[h] + 2)
        dpabs[0] = 1.0
        ok = True
        for j in range(1, h + 1):
            mod = 3 ** j
            # phase for the step into level j depends on P_(j-1): e((2^(s - A + P) mod 3^j)/3^j)
            Pmax_prev = lim[j - 1]
            ph = np.array([cmath.exp(2j * math.pi * float(Fraction(pow(2, s - A + P, mod), mod))) for P in range(Pmax_prev + 1)])
            src = dp[: Pmax_prev + 1] * ph
            srcabs = dpabs[: Pmax_prev + 1]
            new = np.zeros_like(dp)
            newabs = np.zeros_like(dpabs)
            for a in range(1, AMAX + 1):
                w = 2.0 ** (-a)
                hi = min(Pmax_prev + a, lim[j])
                if hi < a:
                    continue
                new[a: hi + 1] += w * src[: hi - a + 1]
                newabs[a: hi + 1] += w * srcabs[: hi - a + 1]
            dp, dpabs = new, newabs
        # only words with total exactly A: P_h == A
        if A <= lim[h]:
            coh += dp[A]
            coh_abs += dpabs[A]
            mass += dpabs[A]
    return coh, coh_abs, mass


if __name__ == "__main__":
    out = []
    for h in (20, 30, 40, 60, 80):
        # resonant s from the closed recursion: scan the window around h log2 3 - 6
        s0 = int(round(h * LOG23)) - 6
        fam = closed_family(h, s0 - 8, s0 + 8)
        s_star = max(fam, key=lambda j: abs(fam[j]))
        full = fam[s_star]
        coh, coh_abs, mass = coherent_split(h, s_star)
        rest = full - coh
        print(f"h={h:3d}: s* = {s_star} (h log2 3 - s* = {h * LOG23 - s_star:.2f}); |full| = {abs(full):.4e}; P_h (no-descent mass) = {mass:.4e}; "
              f"|coh| = {abs(coh):.4e} (|coh|/P_h = {abs(coh)/mass:.4f}, coherent fraction |coh|/sum|terms| = {abs(coh)/coh_abs:.4f}); "
              f"|rest| = {abs(rest):.4e} (|rest|/|full| = {abs(rest)/abs(full):.4f}); arg(coh) - arg(full) = {cmath.phase(coh) - cmath.phase(full):+.3f}")
        sys.stdout.flush()
    # bounded excursions above the critical line: the words with max_j (P_j - j log2 3) < c
    print("== words with bounded excursion above the critical line (P_j < j log2 3 + c for all j) ==")
    for h in (40, 80):
        s0 = int(round(h * LOG23)) - 6
        fam = closed_family(h, s0 - 8, s0 + 8)
        s_star = max(fam, key=lambda j: abs(fam[j]))
        full = fam[s_star]
        for c in (0, 1, 2, 3, 4, 6, 8, 12, 16):
            coh, coh_abs, mass = coherent_split(h, s_star, Awindow=45 + c, c=c)
            print(f"   h={h}, c={c:2d}: mass of the family {mass:.4e} (= {mass / (2 ** (-h)):.3g} x 2^-h); |coh_c| = {abs(coh):.4e} = {abs(coh) / abs(full):.4f} |full|; |full - coh_c|/|full| = {abs(full - coh) / abs(full):.4f}; coherent fraction {abs(coh) / coh_abs:.4f}")
            sys.stdout.flush()
    # random units at h = 30: |mu_hat_30(2^j0)| for random j0 mod L against the resonant window
    h = 30
    L = 2 * 3 ** (h - 1)
    rng = random.Random(3)
    vals = []
    for _ in range(12):
        j0 = rng.randrange(0, L)
        fam = closed_family(h, j0, j0)
        vals.append(abs(fam[j0]))
    s0 = int(round(h * LOG23)) - 6
    fam = closed_family(h, s0 - 6, s0 + 6)
    print(f"h=30: random units |mu_hat| = {['%.2e' % v for v in vals]}; median {np.median(vals):.2e}; square-root scale 0.84 3^(-15) = {0.84 * 3 ** -15:.2e}; resonant window max = {max(abs(v) for v in fam.values()):.2e}")
    print("DONE")
