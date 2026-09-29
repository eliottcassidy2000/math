#!/usr/bin/env python3
"""The floor part of (T1): coherence of the top part's real angle, split by floor levels (S22).

For a word at the resonant exponent s, the top part (levels j > J with kappa_j = s - T_j >= 0) has the real phase
e(sum_(j>J) 2^(kappa_j)/3^j); the bottom part is 2-adically scrambled and depends on the bottom word and the crossing
exponent kappa_c = kappa_(J+1) only.  Conditional on kappa_c the two parts are independent, so the F >= 1 words
(some level with kappa_j > j log2 3 - K, K = 3) cancel iff their top-part sums W(kappa_c, F = 1) = sum 2^-cost
e(angle_top) are small relative to their mass.  Exact DP over (total cost A, level, prefix sum, F, kappa_c) with the
phase applied only at the top levels.  Reports, per kappa_c and F: mass, |W|, coherence |W|/mass; and the totals.
Run: python 04-computation/experiments/collatz_five_mirrors_floor_angle_20260929.py [h ...]
"""
from __future__ import annotations

import cmath
import math
import sys
import time
from fractions import Fraction

import numpy as np

sys.path.insert(0, "04-computation/experiments")
from collatz_five_mirrors_coherent_20260929 import LOG23, closed_family  # noqa: E402

AMAX = 40
KC = 9  # kappa_c index 0..8 = min(kappa_c, 8); 9 = top part not yet started


def floor_angle(h: int, s: int, K: int = 3, Awindow: int = 60):
    Pmax = 3 * h + Awindow
    W = np.zeros((KC + 1, 2), dtype=complex)
    M = np.zeros((KC + 1, 2))
    for A in range(max(h, s - Awindow), s + Awindow + 1):
        dp = np.zeros((Pmax + 1, 2, KC + 1), dtype=complex)
        dpm = np.zeros((Pmax + 1, 2, KC + 1))
        dp[0, 0, KC] = 1.0
        dpm[0, 0, KC] = 1.0
        for j in range(1, h + 1):
            mod = 3 ** j
            Pprev = min(Pmax, 3 * (j - 1) + Awindow)
            src = dp[: Pprev + 1].copy()
            srcm = dpm[: Pprev + 1].copy()
            for P in range(Pprev + 1):
                kappa = s - A + P
                if kappa < 0:
                    continue
                ph = cmath.exp(2j * math.pi * float(Fraction(pow(2, kappa, mod), mod)))
                floor = kappa > j * LOG23 - K
                blk = src[P] * ph          # shape (2, KC+1)
                blkm = srcm[P]
                new = np.zeros_like(blk)
                newm = np.zeros_like(blkm)
                kc = min(kappa, 8)
                # not-started entries start here with kappa_c = kc
                new[:, kc] += blk[:, KC]
                newm[:, kc] += blkm[:, KC]
                new[:, :KC] += blk[:, :KC]
                newm[:, :KC] += blkm[:, :KC]
                if floor:
                    new[1, :] += new[0, :]
                    new[0, :] = 0
                    newm[1, :] += newm[0, :]
                    newm[0, :] = 0
                src[P] = new
                srcm[P] = newm
            nxt = np.zeros_like(dp)
            nxtm = np.zeros_like(dpm)
            for a in range(1, AMAX + 1):
                w = 2.0 ** (-a)
                hi = min(Pprev + a, Pmax)
                if hi < a:
                    continue
                nxt[a: hi + 1] += w * src[: hi - a + 1]
                nxtm[a: hi + 1] += w * srcm[: hi - a + 1]
            dp, dpm = nxt, nxtm
        if A <= Pmax:
            W += dp[A].T
            M += dpm[A].T
    return W, M


if __name__ == "__main__":
    t0 = time.time()
    hs = [int(x) for x in sys.argv[1:]] or [40, 80, 120]
    for h in hs:
        s0 = int(round(h * LOG23)) - 6
        fam = closed_family(h, s0 - 8, s0 + 8)
        s_star = max(fam, key=lambda j: abs(fam[j]))
        W, M = floor_angle(h, s_star)
        print(f"h={h}: s* = {s_star}; mass covered {M.sum():.4f}; never-started (A > s + all levels negative) mass {M[KC].sum():.2e}   ({time.time() - t0:.0f}s)")
        print("   kappa_c: F=0 mass, |W|, coherence  |  F>=1 mass, |W|, coherence  | ratio of coherences (F>=1 / F=0)")
        for kc in range(KC):
            m0, w0 = M[kc, 0], abs(W[kc, 0])
            m1, w1 = M[kc, 1], abs(W[kc, 1])
            c0 = w0 / m0 if m0 > 0 else float("nan")
            c1 = w1 / m1 if m1 > 0 else float("nan")
            print(f"   {kc}{'+' if kc == 8 else ' '}: {m0:.3e} {w0:.3e} {c0:.4f}  |  {m1:.3e} {w1:.3e} {c1:.4f}  |  {c1 / c0 if c0 > 0 else float('nan'):.3f}")
        m0, w0 = M[:KC, 0].sum(), abs(W[:KC, 0].sum())
        m1, w1 = M[:KC, 1].sum(), abs(W[:KC, 1].sum())
        print(f"   totals: F=0 mass {m0:.3e}, |sum W| {w0:.3e}, coherence {w0 / m0:.4f}; F>=1 mass {m1:.3e}, |sum W| {w1:.3e}, coherence {w1 / m1:.4f}; |W_1|/|W_0| = {w1 / w0:.4f}")
        sys.stdout.flush()
    print("DONE")
