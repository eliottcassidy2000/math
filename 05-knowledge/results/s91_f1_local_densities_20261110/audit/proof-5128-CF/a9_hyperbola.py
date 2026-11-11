#!/usr/bin/env python3
"""Auditor proof-5128-CF, check A9: the count of second states used in THM-5128 (D1), state by state.

For a smaller word h1 and a padding s:  G = #{h2 > h1 : x1 x2 + 2^s y1 y2 <= B}  (B = X + [s = 0]),
  Gamma_X(h1) = int_0^1 (Z(v) - max(h1, 1))^+ dv,   Z(v) = X / (2 y1 (rho_1 (1 - v/2) + 2^s)),   int_0^1 Z = X w_s(h1).
Claim: G = Gamma_X + O((Y/y1)^(1/2) + 1), Y = X/2^s.   Printed: (G - Gamma_X) / sqrt(Y/y1) for a sample, and its
maximum over ALL smaller states for s = 0, 1, 2, 5 at X = 2^K.  Also (C1) at p = 5, 7 for h1 = 0, 1, 3:
(#{M = c} - G/p) / (Y/y1)^(1 - theta_p) over the classes c.

usage: python3 a9_hyperbola.py K     (default 24)
"""
import sys
import numpy as np
from numba import njit
from a_lib import psi_table

K = int(sys.argv[1]) if len(sys.argv) > 1 else 24
X = 1 << K
H = (1 << (K - 1)) + 2
psi = psi_table(H)
x = 2 * psi + 1


@njit(cache=True)
def G_count(h1, s, B, x):
    y1 = 2 * h1 + 1
    x1 = x[h1]
    by = y1 << s
    h2max = ((B - 1) // by - 1) // 2
    cnt = 0
    for h2 in range(h1 + 1, h2max + 1):
        if x1 * x[h2] + by * (2 * h2 + 1) <= B:
            cnt += 1
    return cnt


@njit(cache=True)
def G_classes(h1, s, B, x, p):
    y1 = 2 * h1 + 1
    x1 = x[h1]
    by = y1 << s
    h2max = ((B - 1) // by - 1) // 2
    cnt = np.zeros(p, dtype=np.int64)
    for h2 in range(h1 + 1, h2max + 1):
        M = x1 * x[h2] + by * (2 * h2 + 1)
        if M <= B:
            cnt[M % p] += 1
    return cnt


def Gamma(h1, s, Xv):
    """closed form of int_0^1 (Z(v) - m)^+ dv, Z(v) = c0 / (alpha - beta v) increasing in v."""
    y1 = 2 * h1 + 1.0
    rho = x[h1] / y1
    m = max(h1, 1)
    c0 = Xv / (2 * y1)
    alpha = rho + 2.0 ** s
    beta = rho / 2
    # Z(v) >= m  <=>  alpha - beta v <= c0/m  <=>  v >= (alpha - c0/m)/beta
    v0 = (alpha - c0 / m) / beta
    v0 = min(max(v0, 0.0), 1.0)
    if v0 >= 1.0:
        return 0.0
    integral = (c0 / beta) * np.log((alpha - beta * v0) / (alpha - beta))
    return integral - m * (1 - v0)


print("X = 2^%d" % K)
print(" s    h1 |          G      Gamma_X    (G - Gamma)/sqrt(Y/y1) |  X w_s(h1)")
for s in (0, 1, 3):
    B = X + 1 if s == 0 else X
    for h1 in (0, 1, 2, 3, 5, 8, 21, 100, 1000):
        y1 = 2 * h1 + 1
        if (y1 * y1) << s >= B:
            continue
        G = G_count(h1, s, B, x)
        Gam = Gamma(h1, s, float(B))
        Y = X / 2.0 ** s
        xw = X / x[h1] * np.log((x[h1] + 2.0 ** s * y1) / (x[h1] / 2 + 2.0 ** s * y1))
        print(" %d %5d | %10d %12.1f   %+8.4f | %12.1f" % (s, h1, G, Gam, (G - Gam) / np.sqrt(Y / y1), xw))

print("\nmaximum of |G - Gamma_X| / (sqrt(Y/y1) + 1) over all smaller states:")
for s in (0, 1, 2, 5):
    B = X + 1 if s == 0 else X
    worst = 0.0
    arg = -1
    tot = 0
    totg = 0.0
    h1 = 0
    while ((2 * h1 + 1) ** 2) << s < B:
        G = G_count(h1, s, B, x)
        Gam = Gamma(h1, s, float(B))
        r = abs(G - Gam) / (np.sqrt(X / 2.0 ** s / (2 * h1 + 1)) + 1)
        tot += G
        totg += Gam
        if r > worst:
            worst, arg = r, h1
        h1 += 1
    print("  s = %d: max = %.3f at h1 = %d (over %d smaller states);  sum G = %d, sum Gamma_X = %.0f, difference/Y^(3/4) = %+.3f"
          % (s, worst, arg, h1, tot, totg, (tot - totg) / (X / 2.0 ** s) ** 0.75))

rates = {5: 0.7418, 7: 0.9320}
print("\n(C1): max over classes of |#{M = c} - G/p| / (Y/y1)^(1 - theta_p), s = 0")
for p in (5, 7):
    th = np.log(1 / rates[p]) / np.log(2 / rates[p] ** 2)
    for h1 in (0, 1, 3):
        row = []
        for k in range(14, K + 1, 2):
            Bk = (1 << k) + 1
            c = G_classes(h1, 0, Bk, x, p)
            G = c.sum()
            row.append("%.4f" % (np.abs(c - G / p).max() / (2.0 ** k / (2 * h1 + 1)) ** (1 - th)))
        print("  p = %d, h1 = %d (g = %d): k = 14, 16, ..: %s" % (p, h1, np.gcd(x[h1], 2 * h1 + 1), " ".join(row)))
