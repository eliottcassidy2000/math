#!/usr/bin/env python3
"""S91 (THM-5128 (D2)-(D4)): the constants A_p of the 1/ln X term, and the expansion against the census.

A_p^sigma = sum_{h >= 0} chi(h) w^sigma(h),   chi(h) = [p | gcd(x(h), y(h))] - 1/p^2,
w^0(h) = (1/x) ln(2 (x + y)/(x + 2y)),   w^1(h) = (1/x) ln(1 + x/(2y)),
summed by dyadic blocks (the block sums decay geometrically; the series is only conditionally convergent).
THM-5128 (D): D(X) = #T(X)/p^2 + 2 A_p X + O(X^(1 - theta_p/2)), i.e. the share F_small(X) of triples whose smaller
state is p-divisible satisfies  (F_small(X) - 1/p^2) * #T(X) / X  ->  2 A_p.

Printed: partial sums of A_p^0, A_p^1 at 2^20, 2^22, 2^24, 2^26 (NUMERICAL); and, from first_order_exact.out and
first_order_exact_28.out if present, the measured (F_small - 1/p^2) #T / (2X) at X = 2^18 ... 2^28 next to A_p.
Checks: for p = 5, 11, 17 (fast primes with an early divisible state) the last two partial sums agree to 3 %, and
the measured value at the largest X is within 12 % of A_p.
"""
import math
import os
import re
import sys

import numpy as np

PRIMES = [3, 5, 7, 11, 13, 17, 31, 43]


def main():
    nb = 26
    size = 1 << nb
    q = np.zeros(size, dtype=np.int32)
    for k in range(nb):
        lo, hi = 1 << k, 1 << (k + 1)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + np.where((idx & 1) == 1, (half + 1) >> 1, 0).astype(np.int32)
    h = np.arange(size, dtype=np.int64)
    y = 2 * h + 1
    x = y - 2 * q.astype(np.int64)
    xf, yf = x.astype(np.float64), y.astype(np.float64)
    w0 = np.log(2 * (xf + yf) / (xf + 2 * yf)) / xf
    w1 = np.log1p(xf / (2 * yf)) / xf
    del xf, yf
    ok = True
    A = {}
    print(f"{'p':>3} | A_p^0 at 2^20, 2^22, 2^24, 2^26 | A_p^1 at 2^20, 2^22, 2^24, 2^26")
    for p in PRIMES:
        div = ((y % p == 0) & (x % p == 0))
        chi = div.astype(np.float64) - 1.0 / (p * p)
        c0 = np.cumsum(chi * w0)
        c1 = np.cumsum(chi * w1)
        vals0 = [float(c0[(1 << n) - 1]) for n in (20, 22, 24, 26)]
        vals1 = [float(c1[(1 << n) - 1]) for n in (20, 22, 24, 26)]
        A[p] = (vals0[-1], vals1[-1])
        print(f"{p:3d} | " + " ".join(f"{v:+.5f}" for v in vals0) + " | " + " ".join(f"{v:+.5f}" for v in vals1))
        if p in (5, 11, 17):
            ok &= abs(vals0[-1] - vals0[-2]) <= 0.03 * abs(vals0[-1]) and abs(vals1[-1] - vals1[-2]) <= 0.03 * abs(vals1[-1])
    # compare with the census
    rows = {}
    for name in ("first_order_exact.out", "first_order_exact_28.out"):
        if not os.path.exists(name):
            continue
        K = None
        tot = None
        for line in open(name):
            m = re.match(r"X = 2\^(\d+): triples sheet 0: (\d+), sheet 1: (\d+)", line)
            if m:
                K = int(m.group(1))
                tot = (int(m.group(2)), int(m.group(3)))
                continue
            f = line.split()
            if K is not None and len(f) == 7 and f[0] in ("0", "1") and f[1].isdigit():
                sheet, p = int(f[0]), int(f[1])
                fs = float(f[2]) / (p * p)
                rows[(K, sheet, p)] = (fs - 1.0 / (p * p)) * tot[sheet] / (2.0 * (1 << K))
    if rows:
        ks = sorted({k for (k, _, _) in rows})
        print("measured (F_small - 1/p^2) #T / (2X), which tends to A_p:")
        print(f"{'p':>3} {'sheet':>5} | " + " ".join(f"2^{k:<6d}" for k in ks) + " |   A_p")
        for p in PRIMES:
            for sheet in (0, 1):
                if (ks[-1], sheet, p) not in rows:
                    continue
                print(f"{p:3d} {sheet:5d} | " + " ".join(f"{rows[(k, sheet, p)]:+8.5f}" for k in ks) + f" | {A[p][sheet]:+.5f}")
                if p in (5, 11, 17):
                    ok &= abs(rows[(ks[-1], sheet, p)] - A[p][sheet]) <= 0.12 * abs(A[p][sheet])
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
