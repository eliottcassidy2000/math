#!/usr/bin/env python3
"""Auditor bridges-note, script A9: which triples does every_base.py (4) miss?

Counts the ordered triples (h1, h2, s >= 1) with N_b = X1 X2 + b 2^(s-1) y1 y2 <= 2^24 in three ways:
  complete            : table h < 2^23, exact lower-bound pruning
  table 2^21          : same pruning, table cut at h < 2^21 (every_base.py line 220)
  table 2^21 + skip   : in addition, skip the smaller state h1 whenever X1^2 + coef y1^2 > X
                        (every_base.py lines 135-138; unsound because psi is not monotone)
The last line should reproduce the triple counts printed in every_base.out.
"""
import numpy as np
from numba import njit


@njit(cache=False)
def count(psi, b, X, H, skip):
    tot = 0
    s = 1
    while b * (1 << (s - 1)) + 1 <= X:
        coef = b * (1 << (s - 1))
        for h1 in range(H):
            y1 = 2 * h1 + 1
            X1 = 1 + b * psi[h1]
            if X1 + coef * y1 * y1 > X:
                if coef * y1 * y1 > X:
                    break
            if skip and X1 * X1 + coef * y1 * y1 > X:
                continue
            for h2 in range(h1, H):
                y2 = 2 * h2 + 1
                if X1 + coef * y1 * y2 > X:
                    break
                N = X1 * (1 + b * psi[h2]) + coef * y1 * y2
                if N > X:
                    continue
                tot += 1 if h2 == h1 else 2
        s += 1
    return tot


def main():
    X = 1 << 24
    size = 1 << 23
    psi = np.zeros(size, np.int64)
    for h in range(1, size):
        k = h >> 1
        psi[h] = psi[k] + k if h % 2 == 0 else psi[k] + k // 2 + 1
    printed = {2: 55892593, 3: 35516685, 4: 24945520, 5: 18836356, 6: 14916908, 10: 7632149}
    for b in (2, 3, 4, 5, 6, 10):
        c_full = count(psi, b, X, size, False)
        c_21 = count(psi, b, X, 1 << 21, False)
        c_21s = count(psi, b, X, 1 << 21, True)
        print("b = %2d: complete %d ; table 2^21 %d ; table 2^21 + skip %d ; every_base.out %d ; reproduced: %s ; missing %d (%.3f %%)"
              % (b, c_full, c_21, c_21s, printed[b], c_21s == printed[b], c_full - printed[b], 100.0 * (c_full - printed[b]) / c_full))


if __name__ == "__main__":
    main()
