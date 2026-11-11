#!/usr/bin/env python3
"""Auditor proof-5128-CF, check A6: the prime 2 BY VALUE (triples with N <= X), which THM-5128 (E) asserts
("N is asymptotically equidistributed on the odd classes modulo every power of 2, on each sheet") but proves
only through "x and y are independent uniform odd residues in the limit" (a statement about h < 2^n).

usage: python3 a6_two_adic_by_value.py K      (default 24)

Measured, per sheet, for X = 2^k:
  (a) the law of N mod 2^j, j = 2..5, on the odd classes (max relative deviation from uniform, and the noise scale);
  (b) the by-value share of triples whose LARGER / SMALLER word is odd (y = 3 mod 4), and whose x is 3 mod 4;
  (c) the prediction of (b) from the uniform law of rho = x/y on [1/2, 1] (bottom digit of h decides rho <= 3/4):
        larger state :  share(h2 odd) = int_{1/2}^{1} (1/r) ln((3r/4 + 2^s)/(r/2 + 2^s)) dr / int_{1/2}^{1} F_s
        smaller state:  share(h1 odd) = int_{1/2}^{3/4} F_s / int_{1/2}^{1} F_s        (F_s as in (D1)).
"""
import sys
import time
import numpy as np
import mpmath as mp
from a_lib import psi_table, census_mod

K = int(sys.argv[1]) if len(sys.argv) > 1 else 24
KMIN = 16
t0 = time.time()
HN = (1 << (K - 1)) + 2
psi = psi_table(HN)
x64 = 2 * psi + 1
x = x64.astype(np.int32)
h = np.arange(HN, dtype=np.int64)
y64 = 2 * h + 1
cls = (2 * (x64 % 4 == 3) + (y64 % 4 == 3)).astype(np.uint8)
M = 32
hist, diag = census_mod(K, KMIN, x, M, cls, 4)
hist = np.cumsum(hist, axis=1)
diag = np.cumsum(diag, axis=1)

mp.mp.dps = 25


def Fs(s, r):
    return mp.log((r + mp.mpf(2) ** s) / (r / 2 + mp.mpf(2) ** s)) / r


def Fodd(s, r):
    return mp.log((3 * r / 4 + mp.mpf(2) ** s) / (r / 2 + mp.mpf(2) ** s)) / r


half, one, tq = mp.mpf(1) / 2, mp.mpf(1), mp.mpf(3) / 4
den = [mp.quad(lambda r: Fs(s, r), [half, one]) for s in range(0, 60)]
num_large = [mp.quad(lambda r: Fodd(s, r), [half, one]) for s in range(0, 60)]
num_small = [mp.quad(lambda r: Fs(s, r), [half, tq]) for s in range(0, 60)]
pred_large = (float(num_large[0] / den[0]), float(sum(num_large[1:]) / sum(den[1:])))
pred_small = (float(num_small[0] / den[0]), float(sum(num_small[1:]) / sum(den[1:])))
print("limit predictions from the uniform law of rho:  share(larger word odd) = %.5f (sheet 0), %.5f (sheet 1);"
      "  share(smaller word odd) = %.5f, %.5f" % (pred_large[0], pred_large[1], pred_small[0], pred_small[1]))

for sheet in (0, 1):
    print("\nsheet %d" % sheet)
    print("   k |          T | max rel. dev. of N mod 4, 8, 16, 32 on odd classes | noise sqrt(16/T) |"
          " larger word odd | smaller word odd | larger x=3 mod 4 | smaller x=3 mod 4")
    for k in range(KMIN, K + 1, 2):
        hk = hist[sheet, k - KMIN]
        dk = diag[sheet, k - KMIN]
        T = 2 * hk.sum() + dk.sum()
        C = 2 * hk.sum(axis=(0, 1)) + dk.sum(axis=0)
        assert C[0::2].sum() == 0
        devs = []
        for j in (2, 3, 4, 5):
            m = 1 << j
            Cm = np.zeros(m)
            for c in range(M):
                Cm[c % m] += C[c]
            odd = Cm[1::2]
            devs.append(np.abs(odd / (T / (m / 2)) - 1).max())
        pairs = 2 * hk.sum()          # ordered triples with h1 != h2
        ylarge = 2 * hk[:, [1, 3]].sum() / pairs
        ysmall = 2 * hk[[1, 3]].sum() / pairs
        xlarge = 2 * hk[:, [2, 3]].sum() / pairs
        xsmall = 2 * hk[[2, 3]].sum() / pairs
        print("  %2d | %10d | %.5f %.5f %.5f %.5f | %.5f | %.5f | %.5f | %.5f | %.5f"
              % (k, T, devs[0], devs[1], devs[2], devs[3], np.sqrt(16.0 / T), ylarge, ysmall, xlarge, xsmall))
print("\ntotal time %.1f s" % (time.time() - t0))
