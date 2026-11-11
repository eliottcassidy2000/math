#!/usr/bin/env python3
"""Verifier of the corrections to HYP-9660: the exceptions to the worst-prime-power rule above 101 that the corrected
text names (171, 205, 219, 381, 395, 445) are recomputed here from the FULL 2m x 2m operators of THM-5128 (B), with
dense eigenvalues at every k (no estimates), and compared with rates_exact/rates_by_modulus.txt and with the rates
of their proper divisors."""
import cmath
import math
import os
import sys
import time

import numpy as np

EXP = "C:/tmp/wt-s91/04-computation/experiments/s91_f1_local_densities_20261110"


def T_full(m, k):
    inv2 = pow(2, -1, m)
    kk = (k * inv2) % m
    M = np.zeros((2 * m, 2 * m), dtype=complex)
    eta = np.arange(m)
    for beta in (0, 1):
        rows = 2 * eta + beta
        M[rows, 2 * ((2 * eta) % m)] += 1
        M[rows, 2 * ((2 * eta + 1) % m) + 1] += np.exp(2j * math.pi * ((kk * (eta + beta)) % m) / m)
    return M


def rate_of(m):
    best, arg = 0.0, 0
    for k in range(1, (m - 1) // 2 + 1):
        r = float(np.abs(np.linalg.eigvals(T_full(m, k))).max()) / 2
        if r > best:
            best, arg = r, k
    return best, arg


table = {}
with open(os.path.join(EXP, "rates_exact/rates_by_modulus.txt")) as fh:
    for line in fh:
        if not line.startswith("#") and line.strip():
            m, r, k, sg = line.split()
            table[int(m)] = float(r)

ok = True
t0 = time.time()
for m in (171, 205, 219, 381, 395, 445):
    r, k = rate_of(m)
    divs = [d for d in range(3, m, 2) if m % d == 0]
    top = max(table[d] for d in divs)
    good = abs(r - table[m]) < 6e-7 and r > top + 1e-3
    ok &= good
    print(f"m = {m}: max_k rho(T_k)/2 = {r:.6f} at k = +-{k} (table {table[m]:.6f}); proper divisors "
          f"{[(d, round(table[d], 6)) for d in divs]}: above all of them by {r - top:.4f}   "
          f"{'PASS' if good else 'FAIL'}   [{time.time() - t0:.0f} s]", flush=True)
print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
sys.exit(0 if ok else 1)
