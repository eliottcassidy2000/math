#!/usr/bin/env python3
"""Adversarial verifier, lane twoapex: maximum of H over the words of length 16 and 17 (the lane's maxima at these
lengths lie above the histogram bound 2^26 and rest on its single run).  Primitive words only (a word starting with
0 or ending with 3 has the value of a shorter word); own transfer; no histogram.

usage: python v7_max.py THREADS
"""
import sys
import time

import numpy as np
from numba import njit, prange, set_num_threads

from v2_census import prefixes, suffixes


@njit(parallel=True, cache=True)
def scan(al, be, T, bound):
    n1 = al.shape[0]
    n2 = be.shape[0]
    mx = np.zeros(T, np.int64)
    ai = np.zeros(T, np.int64)
    aj = np.zeros(T, np.int64)
    ties = np.zeros(T, np.int64)
    big = np.zeros(T, np.int64)
    for c in prange(T):
        lo = (n1 * c) // T
        hi = (n1 * (c + 1)) // T
        m = np.int64(0)
        bi = 0
        bj = 0
        tie = 0
        nb = 0
        for i in range(lo, hi):
            for j in range(n2):
                h = np.int64(0)
                for s in range(24):
                    h += np.int64(al[i, s]) * np.int64(be[j, s])
                if h >= bound:
                    nb += 1
                    if h > m:
                        m = h
                        bi = i
                        bj = j
                        tie = 1
                    elif h == m:
                        tie += 1
        mx[c] = m
        ai[c] = bi
        aj[c] = bj
        ties[c] = tie
        big[c] = nb
    return mx, ai, aj, ties, big


def main():
    T = int(sys.argv[1])
    set_num_threads(T)
    for n, claimed in ((16, 73529109), (17, 226501279)):
        t0 = time.time()
        n2 = 6
        al = prefixes(n - n2)
        be = suffixes(n2)
        mx, ai, aj, ties, big = scan(al, be, T, 1 << 26)
        k = int(np.argmax(mx))
        m = int(mx.max())
        nt = int(ties[mx == m].sum())
        print(f"n = {n}: maximum over the primitive words {m} (lane: {claimed}; equal: {m == claimed}), attained by {nt} "
              f"primitive word(s); primitive words with H >= 2^26: {int(big.sum())}; ratio to 2 * 3^n: {m / (2 * 3 ** n):.5f}; "
              f"{time.time() - t0:.1f} s", flush=True)


if __name__ == "__main__":
    main()
