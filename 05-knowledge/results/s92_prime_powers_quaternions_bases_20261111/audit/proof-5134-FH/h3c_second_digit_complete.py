#!/usr/bin/env python3
"""Auditor proof-5134-FH.  Side check (the numbers belong to HYP-9665 (C), not to THM-5134 (H)): the statistic
"second digit (N_b - 1)/b modulo b is flat up to ..." printed by every_base.py (4) comes from the truncated census
(words < 2^21, diagonal skip).  Here it is recomputed on the complete census, sheet s >= 1, N_b <= 2^24.
"""
import numpy as np
from numba import njit


def psi_table(nbits):
    q = np.zeros(1 << nbits, dtype=np.int64)
    for k in range(nbits):
        lo, hi = 1 << k, 1 << (k + 1)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + np.where((idx & 1) == 1, (half + 1) >> 1, 0)
    return np.arange(1 << nbits, dtype=np.int64) - q


@njit(cache=True)
def digit2(psi, b, X, out):
    H = psi.shape[0]
    total = 0
    for s in range(1, 63):
        if (b << (s - 1)) > X:
            break
        for h1 in range(H):
            y1 = 2 * h1 + 1
            lb1 = 1 + b * ((h1 + 1) // 2)
            if lb1 * lb1 + (b << (s - 1)) * y1 * y1 > X:
                break
            X1 = 1 + b * psi[h1]
            for h2 in range(h1, H):
                y2 = 2 * h2 + 1
                seam = (b << (s - 1)) * y1 * y2
                if X1 * (1 + b * ((h2 + 1) // 2)) + seam > X:
                    break
                N = X1 * (1 + b * psi[h2]) + seam
                if N <= X:
                    w = 1 if h2 == h1 else 2
                    total += w
                    out[((N - 1) // b) % b] += w
    return total


def main():
    psi = psi_table(24)
    printed = {2: 0.000, 3: 0.040, 4: 0.001, 5: 0.002, 6: 0.052, 10: 0.003}
    print("complete census, sheet s >= 1, N_b <= 2^24: max_c |b * share((N_b - 1)/b = c mod b) - 1|   (every_base.out printed value)")
    for b in (2, 3, 4, 5, 6, 10):
        out = np.zeros(b, dtype=np.int64)
        T = digit2(psi, b, 1 << 24, out)
        d2 = b * out / T
        print(f"  b = {b:2d}: #T = {T}; flat to {abs(d2 - 1).max():.4f}   (printed {printed[b]:.3f});  factors {np.round(d2, 4).tolist()}")


if __name__ == "__main__":
    main()
