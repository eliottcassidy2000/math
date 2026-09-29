#!/usr/bin/env python3
"""The top part's coherence and floor fraction conditional on J (the number of bottom ceiling levels), S22.

Same DP as collatz_five_mirrors_floor_angle_20260929.py with J added to the state (J <= JMAX): for each J, the mass
of the words with F = 0 / F >= 1 floor levels, the top-part sums W_J(F) = sum 2^-cost e(angle_top) and their coherence,
the F = 0 fraction of the J-class, and mass_J(F=0) / (h^(-1/2) e^(-hI)).  Prediction (section 2d): for the small J that
carry the coefficient, the F = 0 event is a ballot event of the top bridge (fraction ~ 1/h), its coherence stays high,
and the F >= 1 words are spread.
Run: python 04-computation/experiments/collatz_five_mirrors_floor_angle_J_20260929.py [h ...]
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
from collatz_five_mirrors_renewal_bound_20260929 import IRATE  # noqa: E402

AMAX = 40
JMAX = 12


def floor_angle_J(h: int, s: int, K: int = 3, Awindow: int = 60):
    Pmax = 3 * h + Awindow
    # state: (P, F, J) with J = number of levels so far with kappa < 0 (0..JMAX, JMAX collects >= JMAX)
    W = np.zeros((2, JMAX + 1), dtype=complex)
    M = np.zeros((2, JMAX + 1))
    for A in range(max(h, s - Awindow), s + Awindow + 1):
        dp = np.zeros((Pmax + 1, 2, JMAX + 1), dtype=complex)
        dpm = np.zeros((Pmax + 1, 2, JMAX + 1))
        dp[0, 0, 0] = 1.0
        dpm[0, 0, 0] = 1.0
        for j in range(1, h + 1):
            mod = 3 ** j
            Pprev = min(Pmax, 3 * (j - 1) + Awindow)
            src = dp[: Pprev + 1].copy()
            srcm = dpm[: Pprev + 1].copy()
            for P in range(Pprev + 1):
                kappa = s - A + P
                if kappa < 0:
                    # ceiling level: J += 1, no top phase
                    blk = src[P]
                    blkm = srcm[P]
                    new = np.zeros_like(blk)
                    newm = np.zeros_like(blkm)
                    new[:, 1:] = blk[:, :-1]
                    new[:, JMAX] += blk[:, JMAX]
                    newm[:, 1:] = blkm[:, :-1]
                    newm[:, JMAX] += blkm[:, JMAX]
                    src[P] = new
                    srcm[P] = newm
                    continue
                ph = cmath.exp(2j * math.pi * float(Fraction(pow(2, kappa, mod), mod)))
                src[P] = src[P] * ph
                if kappa > j * LOG23 - K:
                    src[P, 1, :] += src[P, 0, :]
                    src[P, 0, :] = 0
                    srcm[P, 1, :] += srcm[P, 0, :]
                    srcm[P, 0, :] = 0
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
            W += dp[A]
            M += dpm[A]
    return W, M


if __name__ == "__main__":
    t0 = time.time()
    hs = [int(x) for x in sys.argv[1:]] or [40, 80, 120]
    table = {}
    for h in hs:
        s0 = int(round(h * LOG23)) - 6
        fam = closed_family(h, s0 - 8, s0 + 8)
        s_star = max(fam, key=lambda j: abs(fam[j]))
        W, M = floor_angle_J(h, s_star)
        base = h ** (-0.5) * math.exp(-h * IRATE)
        print(f"h={h}: s* = {s_star}; mass covered {M.sum():.4f}   ({time.time() - t0:.0f}s)")
        print("   J: mass(F=0), coherence(F=0), F=0 fraction of the J-class, mass(F=0)/(h^-1/2 e^-hI), h*fraction  |  mass(F>=1), coherence(F>=1)")
        for J in range(JMAX + 1):
            m0, w0 = M[0, J], abs(W[0, J])
            m1, w1 = M[1, J], abs(W[1, J])
            frac = m0 / (m0 + m1) if m0 + m1 > 0 else float("nan")
            table[(h, J)] = (frac, w0 / m0 if m0 > 0 else float("nan"))
            print(f"   {J:2d}{'+' if J == JMAX else ' '}: {m0:.3e} {w0 / m0 if m0 > 0 else float('nan'):.3f} {frac:.4f} {m0 / base:.3f} {h * frac:.2f}  |  {m1:.3e} {w1 / m1 if m1 > 0 else float('nan'):.3f}")
        sys.stdout.flush()
    print("F = 0 fraction times h, rows J = 0..8, columns h:")
    for J in range(9):
        print(f"   J={J}: " + "  ".join(f"{h * table[(h, J)][0]:.2f}" for h in hs))
    print("DONE")
