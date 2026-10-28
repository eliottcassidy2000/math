"""Third, independent existence check (numpy, single integer key) for sums of four Eisenstein fourth powers.

Key K(x, y) = y * W + x with W = 2^30 (valid while |x| < 2^28 for pair sums, i.e. B <= 30).  n = P + Q with P, Q pair
sums iff K(P) + K(Q) = n (the encoding is injective and additive in this range).  Prints the set of n in [lo, hi]
WITHOUT a representation for bases |a|, |b| <= B, plus small-box minimal B for each n (B = 1..Bmax).
"""
import sys

import numpy as np


def fourth_powers(B):
    a, b = np.meshgrid(np.arange(-B, B + 1, dtype=np.int64), np.arange(-B, B + 1, dtype=np.int64), indexing='ij')
    a, b = a.ravel(), b.ravel()
    s, t = a * a - b * b, 2 * a * b - b * b
    x, y = s * s - t * t, 2 * s * t - t * t
    keys = np.unique(y * (1 << 30) + x)          # dedupe +-alpha (same value)
    return keys


def misses(B, lo, hi):
    f = fourth_powers(B)
    P = np.unique((f[:, None] + f[None, :]).ravel())
    out = []
    for n in range(lo, hi + 1):
        q = n - P
        idx = np.searchsorted(P, q)
        idx[idx == len(P)] = 0
        if not np.any(P[idx] == q):
            out.append(n)
    return out


if __name__ == '__main__':
    lo, hi = 0, 200
    for B in (20, 30):
        m = misses(B, lo, hi)
        print(f'B = {B}: n in [{lo}, {hi}] without a representation ({len(m)}): {m}')
    # minimal box for each n <= 50 (n not = 5 mod 8, n != 6)
    first = {}
    for B in range(1, 21):
        m = set(misses(B, 0, 50))
        for n in range(0, 51):
            if n not in m and n not in first:
                first[n] = B
    print('least B with a representation, n = 0..50:', first)
