#!/usr/bin/env python3
"""S93 coordinator: the degrees Deg_K(P_u) of THM-5141 (half split, s1 = s2 = s), exactly, for small s.

For every P_u = (a, b) of the difference set of the reversal graph of length s with b odd, and every dyadic K:
Deg_K(P_u) = #{(a', b') : |a'|, |b'| < 2^s, b' odd, height of (a' 2^s + b, a 2^s + b') <= K}.
Printed: the maximum, the r^2-weighted mean (r = 2^zeros), against the bounds (D1) 4K^2 and (D2) 3 * 2^s (3 + 2 ln K).
Usage: python3 degree_statistics.py S
"""
import math
import sys

import numpy as np
from numba import njit, prange


@njit(cache=True)
def gcd_nb(a, b):
    while b:
        a, b = b, a % b
    return a


@njit(parallel=True, cache=True)
def degrees(s, A, B, out):
    size = 1 << s
    for idx in prange(A.shape[0]):
        a = A[idx]
        b = B[idx]
        for ap in range(-size + 1, size):
            X = ap * size + b
            ax = abs(X)
            for bp in range(-size + 1, size, 2):       # b' odd (size is even, so -size + 1 is odd)
                Xi = a * size + bp
                ai = abs(Xi)
                g = gcd_nb(ax, ai)
                m = (ax if ax > ai else ai) // g
                lg = 0
                while (1 << lg) < m:
                    lg += 1                             # smallest lg with m <= 2^lg
                if lg <= 2 * s + 1:
                    out[idx, lg] += 1


def main():
    s = int(sys.argv[1])
    n = 3 ** s
    idx = np.arange(n, dtype=np.int64)
    D = np.zeros(n, np.int64)
    R = np.zeros(n, np.int64)
    z = np.zeros(n, np.int64)
    t = idx.copy()
    for i in range(s):
        e = (t % 3) - 1
        t //= 3
        D += e << i
        R += e << (s - 1 - i)
        z += (e == 0)
    keep = (R & 1) == 1
    A, B, Z = D[keep], R[keep], z[keep]
    out = np.zeros((A.shape[0], 2 * s + 3), np.int64)
    degrees(s, A, B, out)
    cum = np.cumsum(out, axis=1)
    w4 = 4.0 ** Z
    print(f"s = {s}: {A.shape[0]} difference vectors with b odd; weights sum {w4.sum():.0f} (6^s / 2 = {6 ** s / 2:.0f})")
    for lg in range(0, s + 1):
        K = 1 << lg
        col = cum[:, lg]
        mx = int(col.max())
        mean = float((w4 * col).sum() / w4.sum())
        d1 = 4 * K * K
        d2 = 3 * (1 << s) * (3 + 2 * math.log(K)) if K <= (1 << s) else float('inf')
        print(f"  K = 2^{lg:2d}: max Deg = {mx:7d}, r^2-weighted mean = {mean:9.2f} = {mean / K:.3f} K; "
              f"(D1) 4K^2 = {d1}, (D2) = {d2:.0f}; max/K = {mx / K:.2f}", flush=True)


if __name__ == "__main__":
    main()
