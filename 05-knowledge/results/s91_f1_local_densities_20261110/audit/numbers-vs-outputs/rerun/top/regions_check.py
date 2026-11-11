#!/usr/bin/env python3
"""S91 (THM-5128 (B4)(ii)): numerical check of the regions lemma with its power-saving rate.

For an odd prime p, a block B_n (n binary digits), b >= a > 0, X and an interval I, the lemma says that for every
residue vector w modulo p
        #{h in B_n, h in I, a x(h) + b y(h) <= X, (x, y) = w}  =  p^-2 #{h in B_n, h in I, a x + b y <= X}  +  O_p(2^(n (1 - theta_p))),
theta_p = ln(1/rho_p)/ln(2/rho_p^2), uniformly in a, b, X, I.  For p = 3, 5, 7, 11 and n = 12, ..., 22 this script
takes 60 random regions per (p, n) (a = x(h1), b = 2^s y(h1) for random h1, s, as in the application; X random so
that the region cuts the block; I a random interval) and prints
        max over regions and w of |deviation| / 2^(n (1 - theta_p)).
The ratio must stay bounded as n grows (checked: the value at n = 22 is at most 3 times the largest value at
n <= 16, and at most 8 in absolute terms).  Rates rho_p are taken from twisted_spectrum.out (0.8660, 0.7418,
0.9320, 0.8587) with a margin of 0.005.
"""
import math
import sys

import numpy as np

RATES = {3: 0.8660254, 5: 0.741760, 7: 0.932039, 11: 0.858710}


def main():
    rng = np.random.default_rng(20261110)
    nb = 22
    size = 1 << nb
    q = np.zeros(size, dtype=np.int64)
    for k in range(nb):
        lo, hi = 1 << k, 1 << (k + 1)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + np.where((idx & 1) == 1, (half + 1) >> 1, 0)
    h = np.arange(size, dtype=np.int64)
    y = 2 * h + 1
    x = y - 2 * q
    ok = True
    print(f"{'p':>3} {'theta_p':>8} | " + " ".join(f"n={n:<2d}   " for n in range(12, 23, 2)))
    for p, rate in RATES.items():
        rho = rate + 0.005
        theta = math.log(1 / rho) / math.log(2 / rho ** 2)
        code = (x % p) * p + (y % p)
        row = []
        for n in range(12, 23, 2):
            lo, hi = 1 << (n - 1), 1 << n
            xs, ys, cs, hs = x[lo:hi], y[lo:hi], code[lo:hi], h[lo:hi]
            worst = 0.0
            for _ in range(60):
                h1 = int(rng.integers(0, 1 << int(rng.integers(1, 12))))
                s = int(rng.integers(0, 6))
                a = int(x[h1])
                b = int(y[h1]) << s
                # X such that a fraction between 0.1 and 0.9 of the block lies in Omega
                vals = a * xs + b * ys
                X = int(np.quantile(vals, rng.uniform(0.1, 0.9)))
                i0 = int(rng.integers(lo, lo + (hi - lo) // 3))
                i1 = int(rng.integers(lo + 2 * (hi - lo) // 3, hi))
                sel = (vals <= X) & (hs >= i0) & (hs < i1)
                cnt = np.bincount(cs[sel], minlength=p * p).astype(np.float64)
                dev = float(np.max(np.abs(cnt - cnt.sum() / (p * p))))
                worst = max(worst, dev / 2 ** (n * (1 - theta)))
            row.append(worst)
        print(f"{p:3d} {theta:8.4f} | " + " ".join(f"{v:7.3f}" for v in row))
        ok &= row[-1] <= 3 * max(row[:3]) and row[-1] <= 8
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
