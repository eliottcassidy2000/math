#!/usr/bin/env python3
"""Auditor prior-art-wording, check a2.

README section 5 (2): "the Euler product has minimum 8/pi^2 = 0.811 of the mean, and the S82 audit measured 0.809 on
the N coprime to 255255".  What does THM-5128 (D3) itself predict for the set {N coprime to 255255 = 3*5*7*11*13*17}?

Own enumeration of F1 (ordered triples; word-pair weight s + 1) on odd N in [2^k, 2^(k+1)); mean of the
representation count over N coprime to 255255 divided by the overall mean, by sheet and in total, against the
limit of THM-5128 (D3):
    sheet 1: prod (1 - 1/p^2);   sheet 0: prod (1 - 1/p^2 + 1/(p (p - 1)))   over p = 3, 5, 7, 11, 13, 17.
"""
import sys
import math
import numpy as np
from numba import njit


@njit(cache=False)
def table(nmax):
    psi = np.zeros(nmax, dtype=np.int64)
    for h in range(1, nmax):
        if h % 2 == 0:
            psi[h] = psi[h // 2] + h // 2
        else:
            k = h // 2
            psi[h] = psi[k] + k // 2 + 1
    return psi


@njit(cache=False)
def enumerate_window(x, y, A, B, r0, r1, w1):
    n = len(x)
    for h1 in range(n):
        y1 = np.int64(y[h1])
        x1 = np.int64(x[h1])
        if y1 >= B + 1:
            break
        for h2 in range(n):
            y2 = np.int64(y[h2])
            if y1 * y2 >= B + 1:
                break
            xx = x1 * np.int64(x[h2])
            N = xx + y1 * y2 - 1
            if A <= N < B:
                r0[(N - A) // 2] += 1
            m = 2 * y1 * y2
            s = 1
            while m < B:
                N = xx + m
                if A <= N < B:
                    r1[(N - A) // 2] += 1
                    w1[(N - A) // 2] += s + 1
                m *= 2
                s += 1


def main():
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    A, B = (1 << k) + 1, (1 << (k + 1))
    nmax = (B // 2) + 2
    psi = table(nmax)
    x = (2 * psi + 1).astype(np.int32)
    del psi
    y = (2 * np.arange(nmax, dtype=np.int64) + 1).astype(np.int32)
    nodd = (B - A + 1) // 2
    r0 = np.zeros(nodd, dtype=np.int32)
    r1 = np.zeros(nodd, dtype=np.int32)
    w1 = np.zeros(nodd, dtype=np.int32)
    enumerate_window(x, y, A, B, r0, r1, w1)
    Ns = A + 2 * np.arange(nodd, dtype=np.int64)
    primes = (3, 5, 7, 11, 13, 17)
    cop = np.ones(nodd, dtype=bool)
    cop2 = np.ones(nodd, dtype=bool)
    for p in primes:
        cop &= (Ns % p != 0)
        cop2 &= (Ns % p != 0) & ((Ns + 1) % p != 0)
    tri = r0 + r1                 # ordered triples
    wp = r0.astype(np.int64) + w1  # ordered word pairs (weight s + 1)
    print(f"odd N in [2^{k}, 2^{k + 1}): {nodd} integers; mean triples {tri.mean():.4f} (sheet 0 {r0.mean():.4f}, "
          f"sheet 1 {r1.mean():.4f}); mean word pairs {wp.mean():.4f} = 2 * {wp.mean() / 2:.4f}")
    print(f"share of N coprime to 255255: {cop.mean():.4f}")
    s1 = math.prod(1 - 1 / p**2 for p in primes)
    s0 = math.prod(1 - 1 / p**2 + 1 / (p * (p - 1)) for p in primes)
    w0_tri, w0_wp = 0.43186124, 0.19640777
    print("mean over N coprime to 255255 / overall mean:")
    print(f"    sheet 1 (triples): measured {r1[cop].mean() / r1.mean():.4f}   limit of THM-5128 (D3) {s1:.4f}")
    print(f"    sheet 0 (triples): measured {r0[cop].mean() / r0.mean():.4f}   limit of THM-5128 (D3) {s0:.4f}")
    print(f"    all triples:       measured {tri[cop].mean() / tri.mean():.4f}   limit {w0_tri * s0 + (1 - w0_tri) * s1:.4f}")
    print(f"    word pairs (r1 of the S82 audit up to the factor 2): measured {wp[cop].mean() / wp.mean():.4f}   "
          f"limit {w0_wp * s0 + (1 - w0_wp) * s1:.4f}   (S82 audit at k = 26: 0.809; 8/pi^2 = {8 / math.pi**2:.4f})")
    print(f"mean over N with N (N + 1) coprime to 255255 / overall mean: triples {tri[cop2].mean() / tri.mean():.4f}, "
          f"word pairs {wp[cop2].mean() / wp.mean():.4f}   (limit: prod (1 - 1/p^2) = {s1:.4f} on both sheets)")
    print("done")


if __name__ == "__main__":
    main()
