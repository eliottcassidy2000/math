#!/usr/bin/env python3
"""Toward (T1)/(T2): the one-step gap of the geometric Gauss sums and the exponent walk (S22, 2026-09-29).

(1) 2-adic reading and gap lemma.  For a unit t mod 3^n put theta = t/3^n in (0,1) and xi = -t 3^-n in Z_2.  Then
        G_n(t) = sum_(r>=1) 2^-r e(theta/2^r + {xi}_r),   {xi}_r = (xi mod 2^r)/2^r,
    exactly (the finite normalised sum equals the infinite one by periodicity).  Classes by v = v_2(t) and xi mod 4:
      v = 1:                       |G| <= sqrt(5)/4 + 1/4 = 0.8090 for every theta;
      v = 0, xi = 1 mod 4:         |G| <= 0.8090 for every theta;
      v = 0, xi = 3 mod 4:         no two-term gap as theta -> 1 (t -> 3^n, i.e. t = -small);
      v >= 2:                      no two-term gap as theta -> 0 (t = 2^v u with u small).
    Verified below over all units for n <= 9, with the observed suprema per class and per theta-range.
(2) Uniform average of the level phase: (1/L_n) sum_(k mod L_n) e((2^k mod 3^n)/3^n) = 0 for n >= 2 (-1/2 for n = 1).
(3) Contribution of the resonant coefficient by the number S of levels the exponent walk k_n = s - A + P_(n-1) spends
    outside the corridor 0 <= k_n < n log_2 3 - K (exact DP with a scrambled-level counter), at h = 40, 80.
Run: python 04-computation/experiments/collatz_five_mirrors_gap_20260929.py
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


def gauss_exact(n: int, t: int) -> complex:
    mod = 3 ** n
    L = 2 * 3 ** (n - 1)
    R = min(L, 60)
    s = 0j
    for r in range(1, R + 1):
        s += 2.0 ** (-r) * cmath.exp(2j * math.pi * float(Fraction((t * pow(2, -r, mod)) % mod, mod)))
    return s / (1 - 2.0 ** (-L)) if L <= 60 else s


def gauss_2adic(n: int, t: int) -> complex:
    mod = 3 ** n
    s = 0j
    for r in range(1, 61):
        m_r = (-t * pow(3, -n, 2 ** r)) % (2 ** r)
        s += 2.0 ** (-r) * cmath.exp(2j * math.pi * (t / (2 ** r * mod) + m_r / 2 ** r))
    return s


if __name__ == "__main__":
    print("== (1) one-step Gauss sums by class: 2-adic reading and gaps ==")
    for n in range(2, 10):
        mod = 3 ** n
        units = [t for t in range(1, mod) if t % 3]
        G = {}
        for t in units:
            G[t] = gauss_2adic(n, t)
        if n <= 5:
            err = max(abs(G[t] - gauss_exact(n, t)) for t in units)
            print(f"   n={n}: 2-adic reading against the modular sum, max error {err:.1e}")
        cls = {"v=1": [], "v=0,xi=1(4)": [], "v=0,xi=3(4),theta<=1/2": [], "v=0,xi=3(4),theta>1/2": [],
               "v>=2,theta<=1/2": [], "v>=2,theta>1/2": []}
        for t in units:
            v = (t & -t).bit_length() - 1
            xi4 = (-t * pow(3, -n, 4)) % 4
            theta = t / mod
            if v == 1:
                cls["v=1"].append(abs(G[t]))
            elif v == 0 and xi4 == 1:
                cls["v=0,xi=1(4)"].append(abs(G[t]))
            elif v == 0:
                cls["v=0,xi=3(4),theta<=1/2" if theta <= 0.5 else "v=0,xi=3(4),theta>1/2"].append(abs(G[t]))
            else:
                cls["v>=2,theta<=1/2" if theta <= 0.5 else "v>=2,theta>1/2"].append(abs(G[t]))
        print(f"   n={n}: sup |G| by class: " + "; ".join(f"{k}: {max(vv):.4f} (n={len(vv)})" for k, vv in cls.items() if vv))
    print("   proved: v = 1, and v = 0 with xi = 1 mod 4, have |G| <= sqrt(5)/4 + 1/4 = %.4f" % (math.sqrt(5) / 4 + 0.25))
    print("== (2) uniform average of the level phase ==")
    for n in range(1, 8):
        mod = 3 ** n
        L = 2 * 3 ** (n - 1)
        avg = sum(cmath.exp(2j * math.pi * float(Fraction(pow(2, k, mod), mod))) for k in range(L)) / L
        print(f"   n={n}: (1/L) sum_k e(2^k mod 3^n / 3^n) = {avg.real:+.6f}{avg.imag:+.6f}i")
    print("== (3) resonant coefficient by the number S of levels outside the corridor 0 <= k_n < n log2 3 - K ==")
    K = 3
    for h in (40, 80):
        s0 = int(round(h * LOG23)) - 6
        fam = closed_family(h, s0 - 8, s0 + 8)
        s_star = max(fam, key=lambda j: abs(fam[j]))
        full = fam[s_star]
        SMAX = 24
        contrib = np.zeros(SMAX + 1, dtype=complex)
        mass = np.zeros(SMAX + 1)
        Pmax = 3 * h + 60
        for A in range(max(h, s_star - 60), s_star + 61):
            dp = np.zeros((Pmax + 1, SMAX + 1), dtype=complex)
            dpm = np.zeros((Pmax + 1, SMAX + 1))
            dp[0, 0] = 1.0
            dpm[0, 0] = 1.0
            for j in range(1, h + 1):
                mod = 3 ** j
                Pprev = min(Pmax, 3 * (j - 1) + 60)
                # phase for the step into level j depends on P_(j-1): exponent k_j = s - A + P_(j-1)
                ks = np.arange(Pprev + 1) + (s_star - A)
                ph = np.array([cmath.exp(2j * math.pi * float(Fraction(pow(2, int(k), mod), mod))) for k in ks])
                outside = ((ks < 0) | (ks >= j * LOG23 - K)).astype(int)
                src = dp[: Pprev + 1] * ph[:, None]
                srcm = dpm[: Pprev + 1]
                # scrambled-level counter shift
                src_sh = np.zeros_like(src)
                srcm_sh = np.zeros_like(srcm)
                inside = outside == 0
                src_sh[inside] = src[inside]
                srcm_sh[inside] = srcm[inside]
                out = ~inside
                src_sh[out, 1:] = src[out, :-1]
                src_sh[out, SMAX] += src[out, SMAX]
                srcm_sh[out, 1:] = srcm[out, :-1]
                srcm_sh[out, SMAX] += srcm[out, SMAX]
                new = np.zeros_like(dp)
                newm = np.zeros_like(dpm)
                for a in range(1, AMAX + 1):
                    w = 2.0 ** (-a)
                    hi = min(Pprev + a, Pmax)
                    if hi < a:
                        continue
                    new[a: hi + 1] += w * src_sh[: hi - a + 1]
                    newm[a: hi + 1] += w * srcm_sh[: hi - a + 1]
                dp, dpm = new, newm
            if A <= Pmax:
                contrib += dp[A]
                mass += dpm[A]
        cum = np.cumsum(contrib)
        print(f"   h={h}: s* = {s_star}, |full| = {abs(full):.4e}, mass covered {mass.sum():.5f}")
        for S in range(SMAX + 1):
            if mass[S] > 0:
                print(f"      S={S:2d}: mass {mass[S]:.3e}  |contribution| {abs(contrib[S]):.3e}  per mass {abs(contrib[S]) / mass[S]:.3e}  |cum|/|full| {abs(cum[S]) / abs(full):.4f}")
    print("DONE")
