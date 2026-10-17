#!/usr/bin/env python3
"""S39 audit: THM-4664 statement 2, sizes of the two remainders E and S_Gamma (alpha != 0), recomputed directly
(complex sums from integer exponent histograms; independent of the session's loops), p = 5, 7, 11, 13.
Also the alpha = 0 cases of statement 1: S(0, 0, gamma) = p^2 and S(0, beta, gamma) = -E + S_Gamma (beta != 0).
Usage: python -B a8c_thm4664_remainder_sizes.py
"""
import numpy as np


def F(x, y, z, p):
    u = (1 + x * y) % p
    return ((u ** 3 * z + y * y * u * (4 + 3 * x * y)) % p, (y + 3 * x * u * u * z + 3 * x * y * y * (4 + 3 * x * y)) % p,
            (2 * x - 3 * x * x * y - x ** 3 * z) % p)


for p in (5, 7, 11, 13):
    inv = [0] + [pow(t, -1, p) for t in range(1, p)]
    zeta = np.exp(2j * np.pi * np.arange(p) / p)
    pts = [F(x, y, z, p) for x in range(p) for y in range(p) for z in range(p)]
    G = [F(x, y, (2 - 3 * x * y) * inv[x] * inv[x] % p, p) for x in range(1, p) for y in range(p)]
    mxE = mxG = 0.0
    worst0 = 0.0
    for al in range(p):
        for be in range(p):
            for ga in range(p):
                h = np.zeros(p)
                for c in range(1, p):
                    for w in range(p):
                        qw = (-w ** 3 + w * w) % p; pw = (-3 * w * w + 4 * w) % p
                        h[(al * qw * inv[c] * inv[c] + be * pw * inv[c] + ga * c) % p] += 1
                E = h @ zeta
                hg = np.zeros(p)
                for (f1, f2, f3) in G:
                    hg[(al * f1 + be * f2) % p] += 1
                SG = hg @ zeta
                if al:
                    mxE = max(mxE, abs(E)); mxG = max(mxG, abs(SG))
                else:
                    hs = np.zeros(p)
                    for (f1, f2, f3) in pts:
                        hs[(be * f2 + ga * f3) % p] += 1
                    S = hs @ zeta
                    if be == 0 and ga == 0:
                        continue
                    pred = p * p if be == 0 else (-E + SG)
                    worst0 = max(worst0, abs(S - pred))
    print('p = %2d: max|E|/p = %.3f, max|S_Gamma|/p = %.3f (alpha != 0); alpha = 0 cases of statement 1: max deviation %.1e'
          % (p, mxE / p, mxG / p, worst0))
