#!/usr/bin/env python3
"""Auditor proof-5128-CF, check A3: THM-5128 (C2), (C3), (D2), (D3) against a direct census (own code).

usage: python3 a3_class_law.py K [nA]      (X = 2^k for k <= K, default 24; A_p from states below 2^nA, default 24)

For each odd prime p in the list:
  C(c)   = #{ordered triples of the sheet with N <= X, N = c mod p}
  D      = #{ordered triples with h1 != h2 whose smaller word has p | gcd(x, y)}
  L      = the same for the larger word
  R(c)   = C(c) - T/p - ([c = c_sigma] - 1/p) D          (claim (C2): O(X^(1 - theta_p/2)))
  exact  : a triple whose smaller state is p-divisible has N = c_sigma mod p        (must hold exactly)
  (C3)   : L/T - 1/p^2 small
  (D2)   : (D - T/p^2)/(2X) against A_p^sigma = sum_h chi(h) w^sigma(h) (own partial sums by blocks)
  (D3)   : share(c_sigma) against mu + (1 - 1/p) (2 A/kappa)/ln X.
"""
import sys
import time
import numpy as np
from numba import njit
from a_lib import psi_table, census_mod

K = int(sys.argv[1]) if len(sys.argv) > 1 else 24
NA = int(sys.argv[2]) if len(sys.argv) > 2 else 24
KMIN = 14
t0 = time.time()
HN = 1 << max(K - 1, NA)
psi = psi_table(HN + 2)
x64 = 2 * psi + 1
x = x64.astype(np.int32)
h = np.arange(HN + 2, dtype=np.int64)
y64 = 2 * h + 1
g = np.gcd(x64, y64)
kappa = (0.16153891726308, 0.21251390923738)
rates = {3: 0.8660, 5: 0.7418, 7: 0.9320, 11: 0.8587, 13: 0.8272, 17: 0.8289, 31: 0.9002}   # from the file, (B)


def theta_half(p):
    r = rates[p]
    return 0.5 * np.log(1 / r) / np.log(2 / r ** 2)


def A_partial(p):
    """sum over h < 2^n of chi(h) w^sigma(h), for n = 0..NA, block by block (float64)."""
    xf = x64[: 1 << NA].astype(np.float64)
    yf = y64[: 1 << NA].astype(np.float64)
    chi = (g[: 1 << NA] % p == 0).astype(np.float64) - 1.0 / p ** 2
    w0 = np.log(2 * (xf + yf) / (xf + 2 * yf)) / xf
    w1 = np.log1p(xf / (2 * yf)) / xf
    c0 = np.cumsum(chi * w0)
    c1 = np.cumsum(chi * w1)
    out0 = [c0[(1 << n) - 1] for n in range(0, NA + 1)]
    out1 = [c1[(1 << n) - 1] for n in range(0, NA + 1)]
    # worst partial sum inside the last blocks (natural order), to see whether the order of summation matters
    lo = 1 << (NA - 4)
    osc0 = float(np.max(np.abs(c0[lo:] - c0[-1])))
    osc1 = float(np.max(np.abs(c1[lo:] - c1[-1])))
    return out0, out1, osc0, osc1


print("X = 2^k, k <= %d; A_p from h < 2^%d" % (K, NA))
for p in (3, 5, 7, 11, 13, 17, 31):
    cls = (g % p == 0).astype(np.uint8)
    hist, diag = census_mod(K, KMIN, x, p, cls, 2)
    hist = np.cumsum(hist, axis=1)
    diag = np.cumsum(diag, axis=1)
    A0, A1, o0, o1 = A_partial(p)
    A = (A0[-1], A1[-1])
    th = theta_half(p)
    print("\np = %d   theta_p/2 = %.4f   A_p^0, A_p^1 at h < 2^%d: %+.5f %+.5f" % (p, th, NA, A[0], A[1]))
    print("   block partial sums A_p^0 at n = 12, 14, ..: " + " ".join("%+.5f" % A0[n] for n in range(12, NA + 1, 2)))
    print("   block partial sums A_p^1 at n = 12, 14, ..: " + " ".join("%+.5f" % A1[n] for n in range(12, NA + 1, 2)))
    print("   max |partial sum - final| over 2^%d <= h < 2^%d (natural order): %.5f %.5f" % (NA - 4, NA, o0, o1))
    for sheet in (0, 1):
        csig = (p - 1) if sheet == 0 else 0
        print("   sheet %d: k |        T   p max|R|/T   max|R|/X^(1-th/2) | forced-class violations | L/T p^2 | (D - T/p^2)/(2X)  A_p |"
              " share(c_sig): measured, (D3) formula, limit" % sheet)
        for k in range(KMIN, K + 1, 2):
            X = float(2 ** k)
            hk = hist[sheet, k - KMIN]
            dk = diag[sheet, k - KMIN]
            T = 2 * hk.sum() + dk.sum()
            C = 2 * hk.sum(axis=(0, 1)) + dk.sum(axis=0)
            D = 2 * hk[1].sum()
            L = 2 * hk[:, 1].sum()
            viol = int(hk[1].sum() - hk[1][:, csig].sum())
            ind = np.zeros(p)
            ind[csig] = 1.0
            R = C - T / p - (ind - 1.0 / p) * D
            mu = 1.0 / p + (1.0 / p ** 2) * (1 - 1.0 / p)
            form = mu + (1 - 1.0 / p) * (2 * A[sheet] / kappa[sheet]) / np.log(X)
            print("            %2d | %10d  %8.5f   %8.5f | %d | %.4f | %+.5f  %+.5f | %.5f %.5f %.5f"
                  % (k, T, p * np.abs(R).max() / T, np.abs(R).max() / X ** (1 - th), viol, L / T * p * p,
                     (D - T / p ** 2) / (2 * X), A[sheet], C[csig] / T, form, mu))
print("\ntotal time %.1f s" % (time.time() - t0))
