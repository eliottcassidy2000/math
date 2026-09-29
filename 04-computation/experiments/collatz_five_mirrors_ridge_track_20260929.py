#!/usr/bin/env python3
"""Tracking the dominant ridge of the negative-power family and testing its mechanism (S22).

Along the ridge path (n, m*_n) found adaptively from (11, 414): the amplitude v = 3^(n/2) |m_n(-m*)|, the one-step
Gauss-sum modulus |G_n(2^-m*)|, the 2-adic class of the frequency t = 2^-m* mod 3^n (v_2(t), xi mod 4, theta = t/3^n)
and the longest run of equal 2-adic digits of 3^-n inside positions [m* - 4, m* + 16] (Lemma G: one-step coherence at
exponent -m is a run of equal digits straddling position m).  The same diagnostics at the control position m* + 40.
Then the ridge followed to level NMAX (window [-WIDTH, 0]) and Ntilde_n = sum_(m<=60) 2^-m |m_n(-m)| to NMAX, to see
the ridge arrive at the small exponents.
Run: python 04-computation/experiments/collatz_five_mirrors_ridge_track_20260929.py [NMAX] [WIDTH]
"""
from __future__ import annotations

import cmath
import math
import sys
import time
from fractions import Fraction

import numpy as np

sys.path.insert(0, "04-computation/experiments")
from collatz_five_mirrors_coherent_20260929 import LOG23  # noqa: E402
from collatz_five_mirrors_renewal_20260929 import closed_family_all, gauss  # noqa: E402


def digit_run(n: int, m: int, lo: int, hi: int) -> tuple[int, int]:
    """Longest run of equal binary digits of the 2-adic 3^-n within positions [lo, hi] (position 1 = 2^0); returns (length, start)."""
    P = hi + 2
    x = pow(3, -n, 2 ** P)
    bits = [(x >> (i - 1)) & 1 for i in range(1, P + 1)]
    best, start = 0, lo
    i = lo
    while i <= hi:
        j = i
        while j + 1 <= hi and bits[j] == bits[i - 1] if False else (j + 1 <= hi and bits[j] == bits[i - 1]):
            j += 1
        L = j - i + 1
        if L > best:
            best, start = L, i
        i = j + 1
    return best, start


if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 320
    WIDTH = int(sys.argv[2]) if len(sys.argv) > 2 else 450
    t0 = time.time()
    fam = closed_family_all(NMAX, -WIDTH, 0)
    print(f"(recursion to level {NMAX}, window [-{WIDTH}, 0], {time.time() - t0:.0f}s)")
    print("ridge path from (11, 414): n, m*, v = 3^(n/2)|m_n(-m*)|, |G_n(2^-m*)|, class (v_2 t, xi mod 4, theta), digit run (len@start) in [m*-4, m*+16]; control at m*+40: v, |G|, run")
    m_prev = 414
    path = {}
    for n in range(11, NMAX + 1):
        row = fam[n]
        lo, hi = max(1, m_prev - 6), min(WIDTH, m_prev + 2)
        mstar = max(range(lo, hi + 1), key=lambda m: abs(row[-m]))
        path[n] = mstar
        m_prev = mstar
        if n <= 40 or n % 5 == 0:
            mod = 3 ** n
            t = pow(2, -mstar, mod)
            v2 = (t & -t).bit_length() - 1
            xi4 = (-t * pow(3, -n, 4)) % 4
            theta = t / mod
            g = abs(gauss(n, -mstar))
            run = digit_run(n, mstar, max(1, mstar - 4), mstar + 16)
            mc = min(WIDTH, mstar + 40)
            gc = abs(gauss(n, -mc))
            runc = digit_run(n, mc, max(1, mc - 4), mc + 16)
            print(f"   n={n:3d}: m*={mstar:3d} v={abs(row[-mstar]) * 3 ** (n / 2):6.2f} |G|={g:.3f} class=({v2},{xi4},{theta:.3f}) run={run[0]:2d}@{run[1]:3d} | control m={mc}: v={abs(row[-mc]) * 3 ** (n / 2):5.2f} |G|={gc:.3f} run={runc[0]:2d}")
        if mstar <= 2:
            print(f"   ridge reached m = {mstar} at n = {n}")
            break
    print("ridge slope: " + ", ".join(f"{(path[b] - path[a]) / (b - a):+.2f} over {a}..{b}" for a, b in ((11, 50), (50, 100), (100, 150), (150, 200), (200, 250), (250, 300)) if a in path and b in path))
    print("Ntilde_n = sum_(m<=60) 2^-m |m_n(-m)| to NMAX, with 3^(n/2) and (log2 3 - 1)^-n normalisations, and N_n = max_(m<=60)")
    for n in range(100, NMAX + 1, 4):
        row = fam[n]
        Nt = sum(2.0 ** (-m) * abs(row[-m]) for m in range(1, 61))
        Nm = max(range(1, 61), key=lambda m: abs(row[-m]))
        print(f"   n={n:3d}: Ntilde 3^(n/2) = {Nt * 3 ** (n / 2):7.3f}   (log2 3 - 1)^-n Ntilde = {Nt / (LOG23 - 1) ** n:.3e}   N_n 3^(n/2) = {abs(row[-Nm]) * 3 ** (n / 2):6.2f} at m = {Nm}   ridge m* = {path.get(n, '-')}")
    Nts = {n: sum(2.0 ** (-m) * abs(fam[n][-m]) for m in range(1, 61)) for n in range(20, NMAX + 1)}
    for lo, hi in ((20, NMAX), (160, NMAX), (200, NMAX)):
        ns = np.arange(lo, hi + 1)
        slope, icpt = np.polyfit(ns, [math.log(Nts[n]) for n in ns], 1)
        print(f"   rate of Ntilde over {lo}..{hi}: {math.exp(slope):.4f}")
    print(f"   max (log2 3 - 1)^-n Ntilde_n over 20..{NMAX}: {max(Nts[n] / (LOG23 - 1) ** n for n in Nts):.3e} at n = {max(Nts, key=lambda n: Nts[n] / (LOG23 - 1) ** n)}")
    print("DONE")
