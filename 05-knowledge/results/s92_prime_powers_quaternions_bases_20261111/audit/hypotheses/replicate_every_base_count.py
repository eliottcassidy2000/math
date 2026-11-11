#!/usr/bin/env python3
# Auditor "hypotheses": why do the triple counts of every_base.out (4) differ from the complete census?
# Replica of the LOOP STRUCTURE of the session's census_b (sheet s >= 1, pairs h1 <= h2 counted twice), with two
# switches: the cap on the words (the session passes q[: 2^21]) and the early "continue" on
#   X1*X1 + coef*y1*y1 > X          (the diagonal value exceeds X)
# which is not a valid exit, because X2 = 1 + b psi(h2) is not monotone in h2: some h2 > h1 can have X2 < X1.
# Counts only (no statistics).
import sys
import numpy as np
from numba import njit
from f1core import build_psi


@njit(cache=True)
def count(psi, b, X, K, use_skip):
    H = psi.shape[0]
    tot = 0
    lost_example_h1 = -1
    for s in range(1, K + 1):
        coef = b << (s - 1)
        for h1 in range(H):
            y1 = 2 * h1 + 1
            X1 = 1 + b * psi[h1]
            if (1 + b * (h1 // 2)) * (1 + b * (h1 // 2)) + coef * y1 * y1 > X:
                break
            if use_skip and X1 * X1 + coef * y1 * y1 > X:
                continue
            for h2 in range(h1, H):
                y2 = 2 * h2 + 1
                if (1 + b * (h2 // 2)) * X1 + coef * y1 * y2 > X:
                    break
                N = X1 * (1 + b * psi[h2]) + coef * y1 * y2
                if N > X:
                    continue
                tot += 1 if h2 == h1 else 2
    return tot


def main():
    K = 24
    X = 1 << K
    print("# replicate_every_base_count.py (auditor hypotheses)")
    print("# sheet s >= 1, N_b <= 2^24: session's counts (every_base.out): b=2 55892593, b=3 35516685, b=4 24945520, b=5 18836356, b=6 14916908, b=10 7632149")
    full = build_psi(1 << 24)
    for b in (2, 3, 4, 5, 6, 10):
        capped = full[: 1 << 21]
        a = count(full, b, X, K, False)
        c = count(capped, b, X, K, False)
        d = count(capped, b, X, K, True)
        e = count(full, b, X, K, True)
        print("b = %2d: complete %d | words < 2^21, valid exits %d | words < 2^21 with the session's 'continue' %d | all words with the 'continue' %d" % (b, a, c, d, e))


if __name__ == "__main__":
    main()
