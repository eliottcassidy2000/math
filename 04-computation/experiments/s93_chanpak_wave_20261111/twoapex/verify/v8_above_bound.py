#!/usr/bin/env python3
"""Adversarial verifier, lane twoapex: the words of length 16 and 17 with H >= 2^26 (not in the lane's value
histogram).  Gives the EXACT collision energy E_17, the exact largest multiplicity at length 17 over ALL N, and the
exact number of distinct values at length 17.  For N >= 2^26:  c_16(N) = P_16(N),  c_17(N) = 2 P_16(N) + P_17(N)
(P_n = primitive words of exact length n; no word of length <= 15 reaches 2^26).

usage: python v8_above_bound.py THREADS
"""
import sys
import time

import numpy as np
from numba import njit, prange, set_num_threads

from v2_census import prefixes, suffixes

BOUND = 1 << 26


@njit(parallel=True, cache=True)
def scan(al, be, hist, bound):
    n1 = al.shape[0]
    n2 = be.shape[0]
    T = hist.shape[0]
    big = np.zeros(T, np.int64)
    for c in prange(T):
        lo = (n1 * c) // T
        hi = (n1 * (c + 1)) // T
        nb = 0
        for i in range(lo, hi):
            for j in range(n2):
                h = np.int64(0)
                for s in range(24):
                    h += np.int64(al[i, s]) * np.int64(be[j, s])
                if h >= bound:
                    nb += 1
                    hist[c, (h - bound) >> 1] += 1
        big[c] = nb
    return big


def main():
    T = int(sys.argv[1])
    set_num_threads(T)
    size17 = ((2 * 3 ** 17 - BOUND) >> 1) + 1
    P = {}
    for n in (16, 17):
        t0 = time.time()
        size = ((2 * 3 ** n - BOUND) >> 1) + 1
        hist = np.zeros((T, size), np.uint16)
        al = prefixes(n - 6)
        be = suffixes(6)
        big = scan(al, be, hist, BOUND)
        p = hist.sum(axis=0, dtype=np.int64)
        assert int(p.sum()) == int(big.sum()), "uint16 overflow"
        P[n] = p
        print(f"n = {n}: primitive words with H >= 2^26: {int(big.sum())}; distinct values: {int(np.count_nonzero(p))}; "
              f"largest multiplicity of a value: {int(p.max())}; {time.time() - t0:.1f} s", flush=True)
        del hist, al, be
    c17 = P[17].copy()
    c17[: len(P[16])] += 2 * P[16]
    e_above = int((c17 * c17).sum())
    e_low = 19376960947461            # verifier's v5_assemble17.out (= lane's z2: 1.937696e13)
    E = e_low + e_above
    print(f"all words of length 17 with H >= 2^26: {int(c17.sum())} (lane: 301885533)")
    print(f"sum over N >= 2^26 of c_17(N)^2 = {e_above}; largest c_17(N) for N >= 2^26: {int(c17.max())} at N = "
          f"{BOUND + 2 * int(np.argmax(c17)) + 1}; distinct values >= 2^26: {int(np.count_nonzero(c17))}")
    print(f"EXACT E_17 = {E}; E_17 (3/16)^17 = {E * (3 / 16) ** 17:.5f}; 4^34 / E_17 = {4 ** 34 / E:.1f} = 3^17 / {3 ** 17 * E / 4 ** 34:.3f}")
    print(f"distinct values of the words of length exactly 17: {33554430 + int(np.count_nonzero(c17))} = "
          f"{(33554430 + int(np.count_nonzero(c17))) / 3 ** 17:.4f} * 3^17")
    # mean multiplicity above the bound by octave-like slices
    for lo, hi in ((1 << 26, 1 << 27), (1 << 27, 2 * 3 ** 17 + 1)):
        seg = c17[(lo - BOUND) >> 1: (hi - BOUND) >> 1]
        print(f"   N in [{lo}, {hi}): odd N {len(seg)}, zeros {int((seg == 0).sum())}, mean multiplicity {seg.mean():.2f}, max {int(seg.max())}")


if __name__ == "__main__":
    main()
