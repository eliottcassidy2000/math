#!/usr/bin/env python3
"""S94 coordinator: first miss N_L of the two-apex family G(2,0) (pre-registered Q1, Q2).

Words of length exactly L over the letters (a_i, b_i) (H(0 w) = H(w), so these have the value set of all words of
length at most L).  36-state transfer of two_apex_transfer.py (the S93 coordinator's eval_w1_w2_two_apex.py, checked
there against brute force).  Every thread writes the byte 1 into a shared map (idempotent writes).
Usage: python3 census_twoapex_bitmap.py L LOG2BOUND
"""
import sys
import time

import numpy as np
from numba import njit, prange

from two_apex_transfer import step, final, brute, count_word  # noqa: F401


@njit(cache=True)
def dfs_mark(start, n, d0, seen, bound):
    stack = np.zeros((n + 1, 36), np.int64)
    for k in range(36):
        stack[d0, k] = start[k]
    letters = np.zeros(n + 1, np.int64)
    out = np.zeros(36, np.int64)
    depth = d0
    letters[depth] = -1
    while depth >= d0:
        letters[depth] += 1
        if letters[depth] == 4:
            depth -= 1
            continue
        l = letters[depth]
        step(stack[depth], l & 1, (l >> 1) & 1, out)
        if depth + 1 == n:
            h = final(out)
            if h < bound:
                seen[h >> 1] = 1
        else:
            for k in range(36):
                stack[depth + 1, k] = out[k]
            depth += 1
            letters[depth] = -1


@njit(parallel=True, cache=True)
def census(n, bound, seen):
    for pre in prange(256):                      # four-letter prefixes
        vec = np.zeros(36, np.int64)
        out = np.zeros(36, np.int64)
        vec[0] = 1
        vec[18] = 1
        for pos in range(4):
            l = (pre >> (2 * pos)) & 3
            step(vec, l & 1, (l >> 1) & 1, out)
            for k in range(36):
                vec[k] = out[k]
        dfs_mark(vec, n, 4, seen, bound)


def main():
    L = int(sys.argv[1])
    K = int(sys.argv[2])
    bound = 1 << K
    seen = np.zeros(bound >> 1, np.uint8)
    t0 = time.time()
    census(L, bound, seen)
    miss = np.flatnonzero(seen == 0)
    odd = [int(2 * i + 1) for i in miss[:12]]
    print(f"L = {L}: words of length {L} (value set of all words of length <= {L}); bound 2^{K}; "
          f"{time.time() - t0:.0f} s", flush=True)
    print(f"  odd N < 2^{K} not attained: {miss.shape[0]}; smallest: {odd}")
    above = [x for x in odd if x > 21]
    if above:
        print(f"  N_{L} = {above[0]}; N_{L} / 3^{L + 1} = {above[0] / 3 ** (L + 1):.4f}")
    else:
        print(f"  N_{L} >= 2^{K} (no miss above 21 below the bound)")


if __name__ == "__main__":
    main()
