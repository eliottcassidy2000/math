#!/usr/bin/env python3
"""Verifier: general (non-affine) arc cubes at order 7 (claim M8), own code.  usage: python v3_general_cubes.py k
Reports longest run of consecutive odd values in one k-cube, number of perfect cubes, and for EVERY X in 3..31 the
number of cubes whose value set contains all odd numbers of [X, 2X] (the lane only looked at X = 24..31)."""
import sys, time, itertools
import numpy as np
from numba import njit, prange
from vlib import *


@njit(cache=True)
def one_set(Hb, arcs, k, N, RX, cover):
    K = 1 << k
    offs = np.zeros(K, np.int64)
    allb = 0
    for t in range(k):
        allb |= 1 << arcs[t]
    for s in range(K):
        o = 0
        for t in range(k):
            if (s >> t) & 1:
                o |= 1 << arcs[t]
        offs[s] = o
    free = (N - 1) ^ allb
    nX = RX.shape[0]
    m = 0
    br = 0; blo = np.uint64(0); bhi = np.uint64(0); perfect = 0
    one = np.uint64(1)
    done = False
    while not done:
        lo = np.uint64(0); hi = np.uint64(0)
        for s in range(K):
            b = np.int64(Hb[m | offs[s]])
            if b < 64:
                lo |= one << np.uint64(b)
            else:
                hi |= one << np.uint64(b - 64)
        for xi in range(nX):
            if (lo & RX[xi]) == RX[xi]:
                cover[xi] += 1
        a = lo; c = hi; run = 0
        while a != np.uint64(0) or c != np.uint64(0):
            run += 1
            a2 = a & ((a >> one) | (c << np.uint64(63)))
            c2 = c & (c >> one)
            a = a2; c = c2
        if run > br:
            br = run; blo = lo; bhi = hi
        if run == K:
            perfect += 1
        m = ((m | ~free) + 1) & free
        if m == 0:
            done = True
    return br, blo, bhi, perfect


@njit(parallel=True, cache=True)
def scan(Hb, subsets, k, N, RX):
    ns = subsets.shape[0]
    nX = RX.shape[0]
    best = np.zeros(ns, np.int64); BLO = np.zeros(ns, np.uint64); BHI = np.zeros(ns, np.uint64)
    cover = np.zeros((ns, nX), np.int64); perfect = np.zeros(ns, np.int64)
    for si in prange(ns):
        br, blo, bhi, pf = one_set(Hb, subsets[si], k, N, RX, cover[si])
        best[si] = br; BLO[si] = blo; BHI[si] = bhi; perfect[si] = pf
    return best, BLO, BHI, cover, perfect


def main():
    k = int(sys.argv[1]); n = 7
    pi, pj = pair_arrays(n)
    H, S = build_table(n, pi, pj)
    Hb = ((H - 1) // 2).astype(np.uint8)
    prs = pairs_of(n); na = len(prs)
    subsets = np.array(list(itertools.combinations(range(na), k)), np.int64)
    Xs = list(range(3, 32))
    RX = np.zeros(len(Xs), np.uint64)
    for i, X in enumerate(Xs):
        r = 0
        for v in range(X, 2 * X + 1):
            if v % 2 == 1:
                r |= 1 << ((v - 1) // 2)
        RX[i] = np.uint64(r)
    t0 = time.time()
    best, BLO, BHI, cover, perfect = scan(Hb, subsets, k, 1 << na, RX)
    j = int(np.argmax(best))
    bits = int(BLO[j]) | (int(BHI[j]) << 64)
    vals = [2 * b + 1 for b in range(128) if (bits >> b) & 1]
    print(f"order 7, k = {k}: {len(subsets)} arc sets x {1 << (na - k)} bases ({time.time()-t0:.1f}s)")
    print(f"   longest run of consecutive odd values in one {k}-cube: {int(best.max())}; example arcs "
          f"{[prs[a] for a in subsets[j]]}, value set {vals}")
    print(f"   perfect cubes (2^k consecutive odd values): {int(perfect.sum())}")
    for i, X in enumerate(Xs):
        print(f"   cubes containing all odd numbers of [{X}, {2*X}]: {int(cover[:, i].sum())}")


if __name__ == '__main__':
    main()
