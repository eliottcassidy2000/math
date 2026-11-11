#!/usr/bin/env python3
"""Auditor proof-5128-CF, check A14: the two forms of THM-5128 (C3) (own census), X = 2^k, k <= K (default 26).

  R_small(c) = share(c) - 1/p - ([c = c_sigma] - 1/p) F_small
  R_sym(c)   = share(c) - 1/p - ([c = c_sigma] - 1/p) (2 f - 1/p^2),   f = share of triples whose FIRST state is divisible
Printed: p max|R_sym|, p max|R_small|, F_small p^2, F_large p^2, and 2f - F_small - F_large (the diagonal term).
The file quotes, at X = 2^26 (sheets 0 / 1): p max|R_sym| = 0.0003/0.0002 (5), 0.002 (17), 0.004 (13), 0.009/0.007 (3),
0.011/0.008 (31), 0.059/0.044 (7); one-sided 0.117 at p = 7; F_large p^2 = 1.47 at p = 7.
"""
import sys
import numpy as np
from a_lib import psi_table, census_mod

K = int(sys.argv[1]) if len(sys.argv) > 1 else 26
KMIN = 18
HN = (1 << (K - 1)) + 2
psi = psi_table(HN)
x64 = 2 * psi + 1
x = x64.astype(np.int32)
h = np.arange(HN, dtype=np.int64)
g = np.gcd(x64, 2 * h + 1)
for p in (3, 5, 7, 11, 13, 17, 31):
    cls = (g % p == 0).astype(np.uint8)
    hist, diag = census_mod(K, KMIN, x, p, cls, 2)
    hist = np.cumsum(hist, axis=1)
    diag = np.cumsum(diag, axis=1)
    print("p = %d" % p)
    for sheet in (0, 1):
        csig = (p - 1) if sheet == 0 else 0
        ind = np.zeros(p)
        ind[csig] = 1.0
        row = []
        for k in range(KMIN, K + 1, 2):
            hk = hist[sheet, k - KMIN]
            dk = diag[sheet, k - KMIN]
            T = 2 * hk.sum() + dk.sum()
            C = (2 * hk.sum(axis=(0, 1)) + dk.sum(axis=0)) / T
            Fs = 2 * hk[1].sum() / T
            Fl = 2 * hk[:, 1].sum() / T
            f = (hk[1].sum() + hk[:, 1].sum() + dk[1].sum()) / T      # first state divisible, ordered triples
            Rs = C - 1 / p - (ind - 1 / p) * Fs
            Ry = C - 1 / p - (ind - 1 / p) * (2 * f - 1 / p ** 2)
            row.append("k=%d: sym %.5f small %.5f Fs %.3f Fl %.3f diag %.1e"
                       % (k, p * np.abs(Ry).max(), p * np.abs(Rs).max(), Fs * p * p, Fl * p * p, 2 * f - Fs - Fl))
        print("  sheet %d: " % sheet + " | ".join(row))
