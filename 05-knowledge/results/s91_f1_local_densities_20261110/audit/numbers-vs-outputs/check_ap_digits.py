#!/usr/bin/env python3
"""Auditor numbers-vs-outputs: partial sums A_p^0, A_p^1 of THM-5128 (D4) at 2^26 with 7 decimals (own numba loop,
compensated by summing per dyadic block in float64), to check the 4- and 5-decimal roundings quoted."""
import math

import numpy as np
from numba import njit

PR = np.array([3, 5, 7, 11, 13, 17, 31, 43], dtype=np.int64)


@njit(cache=False)
def run(nb, primes):
    size = 1 << nb
    q = np.zeros(size, np.int32)
    out = np.zeros((primes.shape[0], 2, nb + 1))
    acc = np.zeros((primes.shape[0], 2))
    for h in range(size):
        if h > 0:
            half = h >> 1
            q[h] = q[half] + (((half + 1) >> 1) if (h & 1) else 0)
        y = 2 * h + 1
        x = y - 2 * q[h]
        w0 = math.log(2.0 * (x + y) / (x + 2.0 * y)) / x
        w1 = math.log1p(x / (2.0 * y)) / x
        for j in range(primes.shape[0]):
            p = primes[j]
            chi = -1.0 / (p * p)
            if y % p == 0 and x % p == 0:
                chi += 1.0
            acc[j, 0] += chi * w0
            acc[j, 1] += chi * w1
        # h + 1 is a power of two: record the partial sum over h < 2^n
        if ((h + 1) & h) == 0:
            n = 0
            t = h + 1
            while t > 1:
                t >>= 1
                n += 1
            for j in range(primes.shape[0]):
                out[j, 0, n] = acc[j, 0]
                out[j, 1, n] = acc[j, 1]
    return out


o = run(26, PR)
for j, p in enumerate(PR.tolist()):
    print(f"p = {p:2d}: A^0 at 2^20, 22, 24, 26 = " + " ".join(f"{o[j, 0, n]:+.7f}" for n in (20, 22, 24, 26))
          + " | A^1 = " + " ".join(f"{o[j, 1, n]:+.7f}" for n in (20, 22, 24, 26)))
