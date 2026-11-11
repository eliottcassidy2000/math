#!/usr/bin/env python3
"""Audit proof-5128-AB: sign and exact extremes of S(h) = (2h+1) phi(h) - 4 q(h) = 2 (2h+1) (rho - 1 + phi/2) on B_n.
Prints, for every n <= NB, the exact minimum and maximum of 2^n S(h) (integers) and where they are attained."""
import sys

import numpy as np

nb = int(sys.argv[1]) if len(sys.argv) > 1 else 24
size = 1 << nb
q = np.zeros(size, dtype=np.int64)
for k in range(nb):
    lo, hi = 1 << k, 1 << (k + 1)
    idx = np.arange(lo, hi, dtype=np.int64)
    half = idx >> 1
    q[lo:hi] = q[half] + (idx & 1) * ((half + 1) >> 1)
allpos = True
for n in range(1, nb + 1):
    lo, hi = 1 << (n - 1), 1 << n
    h = np.arange(lo, hi, dtype=np.int64)
    rev = np.zeros_like(h)
    t = h.copy()
    for _ in range(n):
        rev = (rev << 1) | (t & 1)
        t >>= 1
    Sn = (2 * h + 1) * rev - 4 * q[lo:hi] * (1 << n)
    i0, i1 = int(np.argmin(Sn)), int(np.argmax(Sn))
    allpos &= bool(Sn.min() > 0)
    if n <= 10 or n % 4 == 0:
        print(f"n = {n:2d}: min 2^n S = {int(Sn[i0])} at h = {bin(int(h[i0]))[2:]}  (S = {Sn[i0] / 2 ** n:.6f});  "
              f"max S = {Sn[i1] / 2 ** n:.6f} at h = {bin(int(h[i1]))[2:]}  (S/n = {Sn[i1] / 2 ** n / n:.4f})")
print("S(h) > 0 for every 1 <= h < 2^%d: %s" % (nb, allpos))
