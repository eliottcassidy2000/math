#!/usr/bin/env python3
"""Excursion profile of the resonant Fourier coefficient (S21).

The coefficient mu_hat_h(2^s) = sum_w 2^-A(w) e(2^s Y_h(w)/3^h) is decomposed by the maximal excursion of the word's
prefix sums above the critical line, m(w) = max_j floor(P_j - j log_2 3) (m = -1 collects the strict no-descent words,
m = 0, 1, 2, ... the words that rise up to m bits above the line).  Exact dynamic programme over
(total cost A, level j, prefix sum P, excursion bin m) with the phase factors e((2^(s - A + P_(j-1)) mod 3^j)/3^j).
Output per level h: the mass and the complex contribution of each bin, the cumulative sum over bins <= m against the
full coefficient (from the closed recursion), and the effective width (the smallest m with the remainder below 10%).
Run: python 04-computation/experiments/collatz_five_mirrors_excursion_profile_20260929.py [h ...]
"""
from __future__ import annotations

import cmath
import math
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, "04-computation/experiments")
from collatz_five_mirrors_coherent_20260929 import LOG23, closed_family  # noqa: E402

AMAX = 40


def profile(h: int, s: int, C: int = 40, Awindow: int = 60):
    """Returns (bins m = -1..C, mass per bin, contribution per bin)."""
    NB = C + 2  # index b = m + 1, m = -1..C
    Pmax = int(math.floor(h * LOG23)) + C + 1
    mass = np.zeros(NB)
    contrib = np.zeros(NB, dtype=complex)
    # excursion bin of a prefix sum P at level j: floor(P - j log2 3), clipped to [-1, C]
    for A in range(max(h, s - Awindow), s + Awindow + 1):
        dp = np.zeros((Pmax + 1, NB), dtype=complex)
        dpm = np.zeros((Pmax + 1, NB))
        dp[0, 0] = 1.0
        dpm[0, 0] = 1.0
        for j in range(1, h + 1):
            mod = 3 ** j
            Pprev = min(Pmax, int(math.floor((j - 1) * LOG23)) + C + 1)
            ph = np.array([cmath.exp(2j * math.pi * float(Fraction(pow(2, s - A + P, mod), mod))) for P in range(Pprev + 1)])
            src = dp[: Pprev + 1] * ph[:, None]
            srcm = dpm[: Pprev + 1]
            new = np.zeros_like(dp)
            newm = np.zeros_like(dpm)
            for a in range(1, AMAX + 1):
                w = 2.0 ** (-a)
                hi = min(Pprev + a, Pmax)
                if hi < a:
                    continue
                # target prefix sums P' = P + a for P = 0..hi-a; their excursion bins at level j
                Pp = np.arange(a, hi + 1)
                e = np.floor(Pp - j * LOG23).astype(int)
                e = np.clip(e, -1, C)
                blk = w * src[: hi - a + 1]        # shape (len, NB)
                blkm = w * srcm[: hi - a + 1]
                # new bin = max(old bin, e): for each old bin b (m = b - 1), the target bin is max(b, e + 1)
                for b in range(NB):
                    tgt = np.maximum(b, e + 1)
                    np.add.at(new, (Pp, tgt), blk[:, b])
                    np.add.at(newm, (Pp, tgt), blkm[:, b])
            dp, dpm = new, newm
        if A <= Pmax:
            contrib += dp[A]
            mass += dpm[A]
    return mass, contrib


if __name__ == "__main__":
    hs = [int(x) for x in sys.argv[1:]] or [40, 80]
    for h in hs:
        s0 = int(round(h * LOG23)) - 6
        fam = closed_family(h, s0 - 8, s0 + 8)
        s_star = max(fam, key=lambda j: abs(fam[j]))
        full = fam[s_star]
        C = 30
        mass, contrib = profile(h, s_star, C=C)
        cum = np.cumsum(contrib)
        print(f"h={h}: s* = {s_star}; |full| = {abs(full):.4e}; total mass covered {mass.sum():.6f}")
        print("   m (max excursion bin): mass, |contribution|, phase, |cumulative|/|full|, |full - cumulative|/|full|")
        for b in range(C + 2):
            m = b - 1
            if mass[b] == 0:
                continue
            print(f"      m={m:3d}: mass {mass[b]:.3e}  |c| {abs(contrib[b]):.3e}  arg {cmath.phase(contrib[b]):+.2f}  cum {abs(cum[b]) / abs(full):.4f}  rem {abs(full - cum[b]) / abs(full):.4f}")
        eff = next((b - 1 for b in range(C + 2) if abs(full - cum[b]) / abs(full) < 0.1), None)
        print(f"   effective width (remainder below 10%): m = {eff}")
        sys.stdout.flush()
    print("DONE")
