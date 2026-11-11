import sys
import numpy as np
from numba import njit

@njit(cache=True)
def build_q(size):
    q = np.zeros(size, np.int64)
    for h in range(1, size):
        hh = h >> 1
        q[h] = q[hh] + (((hh + 1) >> 1) if (h & 1) else 0)
    return q

def pairs_with_prefix(x, n, t):
    lo, hi = 1 << (n - 1), 1 << n
    h = np.arange(lo, hi, dtype=np.int64)
    key = ((h >> (n - t)) << 40) | x[lo:hi]
    key.sort()
    d = np.diff(key) != 0
    idx = np.flatnonzero(np.concatenate(([True], d, [True])))
    k = np.diff(idx).astype(np.int64)
    return int((k * (k - 1) // 2).sum()), int(len(k))

for n in (16, 18, 20, 22, 24, 26):
    q = build_q(1 << n)
    x = 2 * np.arange(1 << n, dtype=np.int64) + 1 - 2 * q
    N = 1 << (n - 1)
    row = []
    for frac in (0.0, 0.125, 0.25, 0.3333, 0.5):
        t = max(1, int(np.ceil(frac * n)))
        p, classes = pairs_with_prefix(x, n, t)
        row.append(f"t={t}: {p / N:.4f}")
    # distinct values inside the best prefix family of length n/4: average number of distinct values per prefix
    t = int(np.ceil(n / 4))
    p, classes = pairs_with_prefix(x, n, t)
    print(f"n = {n}: colliding pairs per word with common prefix >= t: " + "; ".join(row)
          + f"; distinct (prefix, value) classes at t = n/4: {classes} of {N} words ({classes / N:.4f})")
