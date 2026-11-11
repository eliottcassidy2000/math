#!/usr/bin/env python3
"""Auditor proof-5128-CF, check A11: THM-5128 (C1) uniformly in the smaller state, and the summation (C2).

For every smaller word h1 with p not dividing g(h1), s in {0, 1, 3}:
    e(h1) = max_c |#{h2 > h1 : M <= B, M = c mod p} - G/p|,      bound claimed: C_p (Y/y1)^(1 - theta_p).
Printed: max over h1 of e(h1)/(Y/y1)^(1 - theta_p) (and the h1 where it is attained), and
    sum over h1 of e(h1), against Y^(1 - theta_p/2)   (the summation step of (C2)).
theta_p from the rates quoted in (B) (the supremum of the admissible exponents).

usage: python3 a11_c1_uniform.py K     (default 24)
"""
import sys
import numpy as np
from numba import njit
from a_lib import psi_table

K = int(sys.argv[1]) if len(sys.argv) > 1 else 24
H = (1 << (K - 1)) + 2
psi = psi_table(H)
x = 2 * psi + 1
hh = np.arange(H, dtype=np.int64)
g = np.gcd(x, 2 * hh + 1)
rates = {3: 0.8660, 5: 0.7418, 7: 0.9320, 13: 0.8272}


@njit(cache=True)
def scan(k, s, p, th, x, g):
    B = (np.int64(1) << k) + (1 if s == 0 else 0)
    Y = float(np.int64(1) << k) / 2.0 ** s
    worst = 0.0
    arg = -1
    tot = 0.0
    cnt = np.zeros(p, dtype=np.int64)
    h1 = 0
    while ((2 * h1 + 1) * (2 * h1 + 1)) << s < B:
        y1 = 2 * h1 + 1
        if g[h1] % p != 0:
            x1 = x[h1]
            by = y1 << s
            h2max = ((B - 1) // by - 1) // 2
            cnt[:] = 0
            G = 0
            for h2 in range(h1 + 1, h2max + 1):
                M = x1 * x[h2] + by * (2 * h2 + 1)
                if M <= B:
                    cnt[M % p] += 1
                    G += 1
            e = 0.0
            for c in range(p):
                d = abs(cnt[c] - G / p)
                if d > e:
                    e = d
            tot += e
            r = e / (Y / y1) ** (1 - th)
            if r > worst:
                worst = r
                arg = h1
        h1 += 1
    return worst, arg, tot, Y


for p in (3, 5, 7, 13):
    r = rates[p]
    th = np.log(1 / r) / np.log(2 / r ** 2)
    print("p = %d, theta_p = %.4f" % (p, th))
    for s in (0, 1, 3):
        row = []
        for k in range(16, K + 1, 2):
            worst, arg, tot, Y = scan(k, s, p, th, x, g)
            row.append("k=%d: max %.3f (h1=%d), sum/Y^(1-th/2) = %.4f" % (k, worst, arg, tot / Y ** (1 - th / 2)))
        print("  s = %d: " % s + " | ".join(row))
