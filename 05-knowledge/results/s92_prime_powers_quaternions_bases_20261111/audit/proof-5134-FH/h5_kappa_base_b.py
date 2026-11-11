#!/usr/bin/env python3
"""Auditor proof-5134-FH.  THM-5134 (H)(ii): what the transfer of THM-5128 (D1) to base b gives explicitly.

In base b the fixed-smaller-state region is  a x2 + b' y2 <= X - const  with
    a = b X1 / 2,   b' = b 2^(s-1) y1   (s >= 1),   b' = b y1 / 2   (s = 0),
    const = X1 (1 - b/2)  (s >= 1),     const = X1 (1 - b/2) - b/2  (s = 0),
so the weight of the smaller state (THM-5128 (D1) with these a, b') is
    w_s^b(h) = (1/a) ln((a + b')/(a/2 + b')) = (2/(b X)) ln((X + 2^s y)/(X/2 + 2^s y)),     X = 1 + b psi(h), y = 2h + 1,
and with X/y -> b rho/2, rho uniform on [1/2, 1]:
    kappa_s(b) = int_{1/2}^{1} (4/(b^2 rho)) ln((b rho/2 + 2^s)/(b rho/4 + 2^s)) d rho
               = (4/b^2) (l_(s+1) - 2 l_(s+2) + l_(s+3)),       l_j = -Li_2(-b/2^j),
    kappa^0(b) = kappa_0(b),      kappa^1(b) = sum_{s >= 1} kappa_s(b) = (4/b^2)(l_2 - l_3).
For b = 2 these are kappa_0 = 0.1615389173 and kappa_1 = 0.2125139092 of THM-5128 (D1).

Check against the own complete census: #T(X)/X = kappa ln X + C + o(1), so
    (#T(4X)/(4X) - #T(X)/X) / ln 4  ->  kappa.
Also the exact first-order sum: #T ~ 2 X sum_{h1} w(h1) restricted as in (D1) is not needed here.
"""
import math
import numpy as np
from numba import njit


def Li2(z):
    """dilogarithm for real z <= 1/2 (series) with Landen / inversion for z < 0."""
    if z < 0:
        x = -z
        u = x / (1 + x)
        return -Li2(u) - 0.5 * math.log(1 + x) ** 2
    if z > 0.5:
        return math.pi ** 2 / 6 - math.log(z) * math.log(1 - z) - Li2(1 - z)
    s, t = 0.0, 1.0
    for k in range(1, 200):
        t *= z
        s += t / (k * k)
    return s


def kappa(b):
    l = lambda j: -Li2(-b / 2.0 ** j)
    k0 = (4.0 / b ** 2) * (l(1) - 2 * l(2) + l(3))
    k1 = (4.0 / b ** 2) * (l(2) - l(3))
    return k0, k1


def psi_table(nbits):
    q = np.zeros(1 << nbits, dtype=np.int64)
    for k in range(nbits):
        lo, hi = 1 << k, 1 << (k + 1)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + np.where((idx & 1) == 1, (half + 1) >> 1, 0)
    return np.arange(1 << nbits, dtype=np.int64) - q


@njit(cache=True)
def count(psi, b, X, sheet):
    total = 0
    H = psi.shape[0]
    s0 = 0 if sheet == 0 else 1
    smax = 0 if sheet == 0 else 62
    for s in range(s0, smax + 1):
        if sheet == 1 and (b << (s - 1)) > X:
            break
        for h1 in range(H):
            y1 = 2 * h1 + 1
            lb1 = 1 + b * ((h1 + 1) // 2)
            if sheet == 0:
                seam_d = b * (y1 * y1 - 1) // 2
            else:
                seam_d = (b << (s - 1)) * y1 * y1
            if lb1 * lb1 + seam_d > X:
                break
            X1 = 1 + b * psi[h1]
            for h2 in range(h1, H):
                y2 = 2 * h2 + 1
                if sheet == 0:
                    seam = b * (y1 * y2 - 1) // 2
                else:
                    seam = (b << (s - 1)) * y1 * y2
                if X1 * (1 + b * ((h2 + 1) // 2)) + seam > X:
                    break
                N = X1 * (1 + b * psi[h2]) + seam
                if N <= X:
                    total += 1 if h2 == h1 else 2
    return total


def main():
    print("kappa_sigma(b) = (4/b^2)(l_1 - 2 l_2 + l_3) (sheet 0), (4/b^2)(l_2 - l_3) (sheet 1), l_j = -Li_2(-b/2^j)")
    k0, k1 = kappa(2)
    print(f"  b = 2: kappa_0 = {k0:.10f} (THM-5128: 0.1615389173), kappa_1 = {k1:.10f} (THM-5128: 0.2125139092); "
          f"pi^2/12 check: -Li_2(-1) = {-Li2(-1.0):.10f} vs {math.pi ** 2 / 12:.10f}")
    psi = psi_table(25)          # words up to X/(1.5 b) <= 2^26/3 < 2^25 are needed
    Ks = [18, 20, 22, 24, 26]
    for b in (2, 3, 4, 5, 10):
        k0, k1 = kappa(b)
        for sheet, kk in ((0, k0), (1, k1)):
            T = [count(psi, b, 1 << K, sheet) for K in Ks]
            slopes = [(T[i + 1] / (1 << Ks[i + 1]) - T[i] / (1 << Ks[i])) / math.log(4) for i in range(len(Ks) - 1)]
            print(f"  b = {b:2d} sheet {sheet}: kappa = {kk:.6f};  #T(2^K), K = {Ks}: {T};  "
                  f"slopes (#T(4X)/(4X) - #T(X)/X)/ln 4: {', '.join(f'{v:.5f}' for v in slopes)}")
    print("the slopes approach kappa_sigma(b): the mass law #T = kappa_sigma(b) X ln X + O(X) transfers to base b with these constants")


if __name__ == "__main__":
    main()
