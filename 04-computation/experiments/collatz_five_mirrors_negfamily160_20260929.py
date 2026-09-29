#!/usr/bin/env python3
"""The negative-power family of the 3-adic Syracuse law to level 160 (S22 obligation 3).

Ntilde_n = sum_(m=1)^60 2^-m |mu_hat_n(2^-m mod 3^n)| and N_n = max_(m<=60) |mu_hat_n(2^-m)| by the closed recursion
(window of exponents [-60, 0], valuation truncation a <= 40; certified against a <= 60 at n <= 80 in _renewal_cert_).
Reports Ntilde_n 3^(n/2), the least-squares rate over sliding windows, and the margin against log2 3 - 1 = 0.58496.
Run: python 04-computation/experiments/collatz_five_mirrors_negfamily160_20260929.py [NMAX]
"""
from __future__ import annotations

import math
import sys

import numpy as np

sys.path.insert(0, "04-computation/experiments")
from collatz_five_mirrors_coherent_20260929 import LOG23  # noqa: E402
from collatz_five_mirrors_renewal_20260929 import closed_family_all  # noqa: E402

if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 160
    fam = closed_family_all(NMAX, -60, 0)
    Nt = {}
    Nm = {}
    for n in range(1, NMAX + 1):
        Nt[n] = sum(2.0 ** (-m) * abs(fam[n][-m]) for m in range(1, 61))
        Nm[n] = max(abs(fam[n][-m]) for m in range(1, 61))
    print("n, Ntilde_n, Ntilde_n 3^(n/2), N_n 3^(n/2), (log2 3 - 1)^-n Ntilde_n")
    for n in range(4, NMAX + 1, 4):
        print(f"   n={n:3d}: {Nt[n]:.3e}  {Nt[n] * 3 ** (n / 2):.3f}  {Nm[n] * 3 ** (n / 2):.3f}  {Nt[n] / (LOG23 - 1) ** n:.3e}")
    for lo, hi in ((20, 80), (40, 120), (80, 160), (20, 160), (100, 160)):
        if hi > NMAX:
            continue
        ns = np.arange(lo, hi + 1)
        ys = np.array([math.log(Nt[n]) for n in ns])
        slope, icpt = np.polyfit(ns, ys, 1)
        print(f"   rate of Ntilde_n over {lo}..{hi}: {math.exp(slope):.4f} (3^-1/2 = 0.5774, critical 0.5850); prefactor {math.exp(icpt):.3f}")
    print("diagnostic: n, v_2(n), argmax m of |m_n(-m)|, and 3^(n/2) |m_n(-m)| at m = 1, 2, 4, 8, 16, 32, 60")
    for n in list(range(60, 70)) + list(range(84, NMAX + 1, 2)):
        row = fam[n]
        am = max(range(1, 61), key=lambda m: abs(row[-m]))
        v2 = (n & -n).bit_length() - 1
        print(f"   n={n:3d} v2={v2}: argmax m={am:2d} ({abs(row[-am]) * 3 ** (n / 2):.2f}); " + " ".join(f"{abs(row[-m]) * 3 ** (n / 2):.2f}" for m in (1, 2, 4, 8, 16, 32, 60)))
    print(f"   max Ntilde_n 3^(n/2) over 20..{NMAX}: {max(Nt[n] * 3 ** (n / 2) for n in range(20, NMAX + 1)):.3f}; max (log2 3 - 1)^-n Ntilde_n over 20..{NMAX}: {max(Nt[n] / (LOG23 - 1) ** n for n in range(20, NMAX + 1)):.3e}")
    print("DONE")
