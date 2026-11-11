#!/usr/bin/env python3
"""Auditor 'citations' (S91): the fast F1 enumeration used in check_singular_series_prime_powers.py is compared
with a plain triple loop (no early exits) at X = 2^11, and its totals at X = 2^22 with THM-5128 (D1)."""
import math
import sys

import numpy as np
from numba import njit


@njit
def psi_table(n):
    psi = np.zeros(n, dtype=np.int64)
    for h in range(1, n):
        k = h >> 1
        if h & 1:
            psi[h] = psi[k] + (k >> 1) + 1
        else:
            psi[h] = psi[k] + k
    return psi


@njit
def census(X, x, y, hmax):
    counts = np.zeros((2, 9), dtype=np.int64)
    s = 0
    while (1 << s) + 1 <= X + 1:
        pw = 1 << s
        lim = X + 1 if s == 0 else X
        sheet = 0 if s == 0 else 1
        for h1 in range(hmax):
            x1 = x[h1]
            y1 = pw * y[h1]
            if (y[h1] // 2) + y1 > lim:
                break
            if x1 + y1 > lim:
                continue
            for h2 in range(hmax):
                M = x1 * x[h2] + y1 * y[h2]
                if M > lim:
                    if (x1 // 2 + y1) * y[h2] > lim:
                        break
                    continue
                counts[sheet, M % 9] += 1
        s += 1
    return counts


def plain(X):
    hmax = X // 2 + 2
    psi = [0] * hmax
    for h in range(1, hmax):
        k = h >> 1
        psi[h] = psi[k] + ((k >> 1) + 1 if h & 1 else k)
    xs = [2 * p + 1 for p in psi]
    ys = [2 * h + 1 for h in range(hmax)]
    counts = [[0] * 9, [0] * 9]
    for s in range(0, X.bit_length() + 1):
        lim = X + 1 if s == 0 else X
        for h1 in range(hmax):
            for h2 in range(hmax):
                M = xs[h1] * xs[h2] + (ys[h1] * ys[h2] << s)
                if M <= lim:
                    counts[0 if s == 0 else 1][M % 9] += 1
    return counts


ok = True
for X in (1 << 9, 1 << 11):
    hmax = X // 2 + 2
    psi = psi_table(hmax)
    x = 2 * psi + 1
    y = 2 * np.arange(hmax, dtype=np.int64) + 1
    fast = census(X, x, y, hmax)
    slow = plain(X)
    same = [list(map(int, fast[0])), list(map(int, fast[1]))] == slow
    ok = ok and same
    print("X = %d: fast enumeration equals the plain triple loop: %s  (totals %d, %d)" % (X, same, int(fast[0].sum()), int(fast[1].sum())))
X = 1 << 22
k0, k1 = 0.1615389173, 0.2125139092
print("X = 2^22: leading terms kappa_s X ln X = %.4g, %.4g (THM-5128 (D1); the census gave 10532579 and 13218857, difference O(X) = %d)"
      % (k0 * X * math.log(X), k1 * X * math.log(X), X))
sys.exit(0 if ok else 1)
