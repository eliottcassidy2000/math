#!/usr/bin/env python3
"""S91 (THM-5128 (C)): direct test of the first-order by-value law, in its two equivalent forms.

Among the triples (h1, h2, s) of a sheet with N <= X, with c_s the distinguished class:
  symmetric form   share(c) = 1/p + ([c = c_s] - 1/p) (2 f(X) - 1/p^2)  + R_sym(c),   f = share with h1 p-divisible
  one-sided form   share(c) = 1/p + ([c = c_s] - 1/p) F_small(X)        + R_small(c), F_small = share whose SMALLER
                                                                                     state is p-divisible
(2f = F_small + F_large; the two forms differ by ([c = c_s] - 1/p)(F_large - 1/p^2), and F_large -> 1/p^2 at the
mixing rate).  THM-5128 (C) says both remainders are O(X^-delta).  This script enumerates all triples for X = 2^K
and prints, per sheet and prime: F_small p^2, F_large p^2 and  p max_c |R_sym(c)|,  p max_c |R_small(c)|  with the
sampling scale sqrt(p / #triples).  (The squarefree product law of THM-5128 (C) is not checked by this script;
the independent audit checked it at several squarefree moduli: audit/recompute-5128/a11_squarefree.out.)

Checked (exit code): every triple with a p-divisible state lies in the distinguished class (exact);
p max|R_sym| < 0.005 at p = 5 for the largest X; and at the largest X the symmetric remainder is not larger than the
one-sided one at the slow primes 3 and 7.  Everything else is printed as data (the remainders oscillate: the
dominant eigenvalues of the twisted operators are complex).

Usage: python3 first_order_exact.py K1 [K2 ...]        (default 18 20 22 24 26)
"""
import math
import sys

import numpy as np
from numba import njit

PRIMES = np.array([3, 5, 7, 11, 13, 17, 31], dtype=np.int64)


@njit(cache=False)
def build(size, primes):
    q = np.zeros(size, np.int32)
    mask = np.zeros(size, np.int64)
    for h in range(size):
        if h > 0:
            half = h >> 1
            q[h] = q[half] + (((half + 1) >> 1) if (h & 1) else 0)
        y = 2 * h + 1
        x = y - 2 * q[h]
        m = 0
        for j in range(primes.shape[0]):
            if y % primes[j] == 0 and x % primes[j] == 0:
                m |= (1 << j)
        mask[h] = m
    return q, mask


@njit(cache=False)
def enumerate_all(q, mask, X, K, primes):
    np_ = primes.shape[0]
    cnt = np.zeros((2, np_, 2, 2, 32), np.int64)      # sheet, prime, smaller divisible, larger divisible, class
    H = q.shape[0]
    for s in range(K + 1):
        sheet = 0 if s == 0 else 1
        d0 = 1 if s == 0 else 0
        for h1 in range(H):
            y1 = 2 * h1 + 1
            Y = y1 << s
            if Y > X:
                break
            x1 = y1 - 2 * q[h1]
            y2max = (2 * (X + d0)) // (2 * Y + x1)
            if y2max < 1:
                continue
            h2max = (y2max - 1) // 2 + 1
            if h2max >= H:
                h2max = H - 1
            for h2 in range(0, h2max + 1):
                y2 = 2 * h2 + 1
                x2 = y2 - 2 * q[h2]
                N = x1 * x2 + Y * y2 - d0
                if N <= X:
                    if h1 <= h2:
                        ms, ml = mask[h1], mask[h2]
                    else:
                        ms, ml = mask[h2], mask[h1]
                    for j in range(np_):
                        cnt[sheet, j, (ms >> j) & 1, (ml >> j) & 1, N % primes[j]] += 1
    return cnt


def main():
    ks = [int(a) for a in sys.argv[1:]] or [18, 20, 22, 24, 26]
    kmax = max(ks)
    q, mask = build((1 << kmax) // 2 + 2, PRIMES)
    ok = True
    last = {}
    for K in ks:
        X = 1 << K
        cnt = enumerate_all(q, mask, X, K, PRIMES)
        print(f"X = 2^{K}: triples sheet 0: {int(cnt[0, 0].sum())}, sheet 1: {int(cnt[1, 0].sum())}")
        print(f"  {'sheet':>5} {'p':>3} {'F_small p^2':>11} {'F_large p^2':>11} {'p max|R_sym|':>12} {'p max|R_small|':>14} {'noise':>7}")
        for sheet in (0, 1):
            for j, p in enumerate(PRIMES.tolist()):
                c = cnt[sheet, j, :, :, :p].astype(np.float64)
                tot = c.sum()
                cs = (p - 1) if sheet == 0 else 0
                div = c[1, :, :].sum(axis=0) + c[0, 1, :]
                ok &= abs(div.sum() - div[cs]) < 0.5                  # divisible triples sit in the distinguished class
                f_small = c[1].sum() / tot
                f_large = c[:, 1].sum() / tot
                share = c.sum(axis=(0, 1)) / tot
                ind = np.array([1.0 if k == cs else 0.0 for k in range(p)])
                r_sym = share - 1 / p - (ind - 1 / p) * (f_small + f_large - 1 / p ** 2)
                r_small = share - 1 / p - (ind - 1 / p) * f_small
                a, b = p * np.abs(r_sym).max(), p * np.abs(r_small).max()
                last[(sheet, p)] = (a, b)
                print(f"  {sheet:5d} {p:3d} {f_small * p * p:11.3f} {f_large * p * p:11.3f} {a:12.5f} {b:14.5f} {math.sqrt(p / tot):7.5f}")
    for sheet in (0, 1):
        ok &= last[(sheet, 5)][0] < 0.005
        for p in (3, 7):
            ok &= last[(sheet, p)][0] <= last[(sheet, p)][1]
    print("exactness of the divisible triples, p max|R_sym| < 0.005 at p = 5, and R_sym <= R_small at p = 3, 7 (largest X):", ok)
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
