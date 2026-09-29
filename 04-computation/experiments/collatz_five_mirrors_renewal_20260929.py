#!/usr/bin/env python3
"""Renewal structure of the resonant coefficient on the powers of two (S22, 2026-09-29).

Exponent walk: mu_hat_h(2^s) = E prod_(j=1)^h omega_j(kappa_j), omega_j(k) = e((2^k mod 3^j)/3^j), kappa_j = s - T_j,
T_j = a_j + ... + a_h (suffix sums; from the top the exponent decreases by i.i.d. geometric steps).  Levels with
kappa_j < 0 (ceiling-scrambled: modular inverse powers) and levels with kappa_j > j log2 3 - K (floor-scrambled:
wrapped powers) carry phases far from 1; the corridor 0 <= kappa_j < j log2 3 - K carries phases within 2 pi 2^-K of 1.
(C1) profile |m_n(k)| = |mu_hat_n(2^k mod 3^n)| for k in [-60, floor(n log2 3) + 60] at n = 20, 40, 60, 80;
(C2) N_n = sup_(1<=m<=60) |m_n(-m)| (the negative-power family) for n <= 80, against 3^(-n/2) and M(n);
(C3) exact decomposition of the resonant coefficient by J = #{j : kappa_j < 0} and by F = #{j : kappa_j >= 0,
     kappa_j > j log2 3 - K}, K = 3, at h = 40, 80 (mass, |contribution|, per-mass ratio, cumulative);
(C4) mean and geometric mean of the one-step Gauss-sum modulus |G_n(2^k)| over the scrambled exponent windows
     k in [-40, -1] and k in [floor(n log2 3) + 1, floor(n log2 3) + 40], against the random-unit mean 0.543;
(C5) the real-phase identity e(2^s Y_h/3^h) = e(2^(s-A) C_w(3,2)/3^h) for A <= s, checked at h = 6.
Run: python 04-computation/experiments/collatz_five_mirrors_renewal_20260929.py
"""
from __future__ import annotations

import cmath
import itertools
import math
import sys
import time
from fractions import Fraction

import numpy as np

sys.path.insert(0, "04-computation/experiments")
from collatz_five_mirrors_coherent_20260929 import LOG23, closed_family  # noqa: E402

AMAX = 40


def closed_family_all(HMAX: int, j_lo: int, j_hi: int):
    """{n: {j: m_n(j)}} for j in [j_lo, j_hi], every level n <= HMAX (window valid at every level)."""
    weights = [2.0 ** (-a) for a in range(1, AMAX + 1)]
    prev = {j: 1.0 + 0j for j in range(j_lo - AMAX * HMAX, j_hi + 1)}
    out = {}
    for n in range(1, HMAX + 1):
        mod = 3 ** n
        lo = j_lo - AMAX * (HMAX - n)
        pw = {j: pow(2, j, mod) for j in range(lo - AMAX, j_hi + 1)}
        ph = {j: cmath.exp(2j * math.pi * float(Fraction(pw[j], mod))) for j in pw}
        cur = {}
        for j in range(lo, j_hi + 1):
            s = 0j
            for a in range(1, AMAX + 1):
                s += weights[a - 1] * ph[j - a] * prev.get(j - a, 0j)
            cur[j] = s
        prev = cur
        out[n] = {j: cur[j] for j in range(j_lo, j_hi + 1)}
    return out


def gauss(n: int, k: int) -> complex:
    mod = 3 ** n
    return sum(2.0 ** (-a) * cmath.exp(2j * math.pi * float(Fraction(pow(2, k - a, mod), mod))) for a in range(1, 61))


def decompose(h: int, s: int, kind: str, K: int = 3, SMAX: int = 30, Awindow: int = 60):
    """Exact DP over (A, level, prefix sum, counter): contribution and mass of the words by the counter value.
    kind J: count levels with kappa_j < 0; kind F: levels with kappa_j >= 0 and kappa_j > j log2 3 - K."""
    contrib = np.zeros(SMAX + 1, dtype=complex)
    mass = np.zeros(SMAX + 1)
    Pmax = 3 * h + Awindow
    for A in range(max(h, s - Awindow), s + Awindow + 1):
        dp = np.zeros((Pmax + 1, SMAX + 1), dtype=complex)
        dpm = np.zeros((Pmax + 1, SMAX + 1))
        dp[0, 0] = 1.0
        dpm[0, 0] = 1.0
        for j in range(1, h + 1):
            mod = 3 ** j
            Pprev = min(Pmax, 3 * (j - 1) + Awindow)
            ks = np.arange(Pprev + 1) + (s - A)
            ph = np.array([cmath.exp(2j * math.pi * float(Fraction(pow(2, int(k), mod), mod))) for k in ks])
            if kind == "J":
                counted = ks < 0
            else:
                counted = (ks >= 0) & (ks > j * LOG23 - K)
            src = dp[: Pprev + 1] * ph[:, None]
            srcm = dpm[: Pprev + 1]
            src_sh = np.zeros_like(src)
            srcm_sh = np.zeros_like(srcm)
            inside = ~counted
            src_sh[inside] = src[inside]
            srcm_sh[inside] = srcm[inside]
            src_sh[counted, 1:] = src[counted, :-1]
            src_sh[counted, SMAX] += src[counted, SMAX]
            srcm_sh[counted, 1:] = srcm[counted, :-1]
            srcm_sh[counted, SMAX] += srcm[counted, SMAX]
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
    return mass, contrib


if __name__ == "__main__":
    t0 = time.time()
    print("== (C5) real-phase identity at h = 6: e(2^s Y/3^h) = e(2^(s-A) C_w(3,2)/3^h) for A <= s ==")
    h = 6
    worst = 0.0
    for w in itertools.product(range(1, 5), repeat=h):
        A = sum(w)
        P = [0]
        for a in w:
            P.append(P[-1] + a)
        C = sum(3 ** (h - j) * 2 ** P[j - 1] for j in range(1, h + 1))
        Y = (pow(2, -A, 3 ** h) * C) % 3 ** h
        for s in range(A, A + 6):
            lhs = cmath.exp(2j * math.pi * float(Fraction((pow(2, s, 3 ** h) * Y) % 3 ** h, 3 ** h)))
            rhs = cmath.exp(2j * math.pi * float(Fraction(2 ** (s - A) * C, 3 ** h)))
            worst = max(worst, abs(lhs - rhs))
    print(f"   max |lhs - rhs| over all words with letters <= 4 and s in [A, A+5]: {worst:.1e}")
    print("== (C1) profile |m_n(k)| on the powers of two, k from -60 to floor(n log2 3) + 60 ==")
    HM = 80
    fam = closed_family_all(HM, -60, int(math.floor(HM * LOG23)) + 60)
    print(f"   (recursion to level {HM} done, {time.time() - t0:.0f}s)")
    for n in (20, 40, 60, 80):
        row = fam[n]
        kc = int(math.floor(n * LOG23))
        ks = [k for k in range(-60, kc + 61)]
        vals = [abs(row[k]) for k in ks]
        kmax = ks[int(np.argmax(vals))]
        print(f"   n={n}: argmax k = {kmax} (= n log2 3 - {n * LOG23 - kmax:.2f}), max {max(vals):.3e}; 3^(-n/2) = {3.0 ** (-n / 2):.3e}")
        for lo, hi, tag in ((-60, -41, "k in [-60,-41]"), (-40, -21, "k in [-40,-21]"), (-20, -1, "k in [-20,-1]"),
                            (0, kc // 2, "k in [0, kc/2]"), (kc // 2 + 1, kc - 12, "k in (kc/2, kc-12]"),
                            (kc - 11, kc, "k in [kc-11, kc]"), (kc + 1, kc + 20, "k in [kc+1, kc+20]"),
                            (kc + 21, kc + 60, "k in [kc+21, kc+60]")):
            seg = [abs(row[k]) for k in range(lo, hi + 1)]
            print(f"      {tag:22s}: rms {math.sqrt(sum(v * v for v in seg) / len(seg)):.3e}  max {max(seg):.3e}  min {min(seg):.3e}")
        print("      k = kc-12..kc+8: " + " ".join(f"{abs(row[k]):.2e}" for k in range(kc - 12, kc + 9)))
    print("== (C2) negative-power family: N_n = sup_(1<=m<=60) |m_n(-m)| against 3^(-n/2) and M(n) ==")
    for n in range(2, HM + 1):
        row = fam[n]
        kc = int(math.floor(n * LOG23))
        N = max(abs(row[-m]) for m in range(1, 61))
        margm = max(range(1, 61), key=lambda m: abs(row[-m]))
        Mn = max(abs(row[k]) for k in range(-60, kc + 61))
        if n <= 12 or n % 4 == 0:
            print(f"   n={n:2d}: N_n = {N:.3e} (at m = {margm:2d}), N_n 3^(n/2) = {N * 3 ** (n / 2):.3f}, M(n) = {Mn:.3e}, |m_n(-1)| = {abs(row[-1]):.3e}, |m_n(-2)| = {abs(row[-2]):.3e}")
    print("== (C4) one-step Gauss-sum modulus on the scrambled exponent windows ==")
    for n in (10, 20, 30, 40):
        kc = int(math.floor(n * LOG23))
        for lo, hi, tag in ((-40, -1, "ceiling k in [-40,-1]"), (kc + 1, kc + 40, "floor k in [kc+1,kc+40]"), (0, kc - 6, "corridor k in [0,kc-6]")):
            g = [abs(gauss(n, k)) for k in range(lo, hi + 1)]
            print(f"   n={n}: {tag:26s}: mean |G| {np.mean(g):.4f}, geometric mean {math.exp(np.mean(np.log(g))):.4f}, min {min(g):.4f}, max {max(g):.4f}")
    print("== (C3) decomposition of the resonant coefficient by J (ceiling-scrambled levels) and F (floor-scrambled levels, K = 3) ==")
    for h in (40, 80):
        s0 = int(round(h * LOG23)) - 6
        famh = closed_family(h, s0 - 8, s0 + 8)
        s_star = max(famh, key=lambda j: abs(famh[j]))
        full = famh[s_star]
        for kind in ("J", "F"):
            mass, contrib = decompose(h, s_star, kind)
            cum = np.cumsum(contrib)
            print(f"   h={h}, kind={kind}: s* = {s_star}, |full| = {abs(full):.4e}, mass covered {mass.sum():.5f}  ({time.time() - t0:.0f}s)")
            for S in range(31):
                if mass[S] > 0:
                    print(f"      {kind}={S:2d}: mass {mass[S]:.3e}  |contribution| {abs(contrib[S]):.3e}  arg {cmath.phase(contrib[S]):+.2f}  per mass {abs(contrib[S]) / mass[S]:.3e}  |cum|/|full| {abs(cum[S]) / abs(full):.4f}")
            sys.stdout.flush()
    print("DONE")
