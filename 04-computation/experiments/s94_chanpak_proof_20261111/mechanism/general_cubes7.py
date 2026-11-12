#!/usr/bin/env python3
"""General (not necessarily affine) arc cubes at order 7, exhaustive and exact.
A k-cube = k arcs (any positions) and a base tournament; its 2^k members are obtained by reversing subsets of the arcs.
Question ("positional system with carries"): does some k-cube contain ALL odd numbers of a dyadic interval [X, 2X]?
At order 7 the strong values are >= 25 and 63 is not a value, so only 24 <= X <= 31 is possible, and k >= 4 is needed.
Also reported: the longest run of consecutive odd numbers inside one k-cube (k = 3, 4), and the number of 'perfect'
k-cubes (2^k distinct consecutive odd values).   usage: python general_cubes7.py k"""
import sys, time, itertools
import numpy as np
from numba import njit, prange
from arc_cube import table

@njit(parallel=True, cache=True)
def scan(H, subsets, k, na, RX):
    ns = subsets.shape[0]
    N = 1 << na
    best_run = np.zeros(ns, np.int64); best_start = np.zeros(ns, np.int64)
    cover = np.zeros((ns, RX.shape[0]), np.int64); perfect = np.zeros(ns, np.int64)
    for si in prange(ns):
        bits = np.zeros(k, np.int64)
        allb = 0
        for t in range(k):
            bits[t] = 1 << subsets[si, t]
            allb |= bits[t]
        br = 0; bs = 0
        for m in range(N):
            if m & allb:
                continue
            lo = np.uint64(0); hi = np.uint64(0)
            for s in range(1 << k):
                mm = m
                for t in range(k):
                    if (s >> t) & 1:
                        mm |= bits[t]
                v = H[mm]
                b = (v - 1) // 2
                if b < 64:
                    lo |= np.uint64(1) << np.uint64(b)
                else:
                    hi |= np.uint64(1) << np.uint64(b - 64)
            for xi in range(RX.shape[0]):
                if (lo & RX[xi]) == RX[xi]:
                    cover[si, xi] += 1
            # longest run of consecutive set bits over lo:hi
            run = 0; start = 0; cur = 0; cs = 0
            for b in range(128):
                if b < 64:
                    bit = (lo >> np.uint64(b)) & np.uint64(1)
                else:
                    bit = (hi >> np.uint64(b - 64)) & np.uint64(1)
                if bit:
                    if cur == 0:
                        cs = b
                    cur += 1
                    if cur > run:
                        run = cur; start = cs
                else:
                    cur = 0
            if run > br:
                br = run; bs = 2 * start + 1
            if run == (1 << k):
                perfect[si] += 1
        best_run[si] = br; best_start[si] = bs
    return best_run, best_start, cover, perfect

def main():
    k = int(sys.argv[1]); n = 7
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    pi = np.array([p[0] for p in pairs], np.int64); pj = np.array([p[1] for p in pairs], np.int64)
    H, S = table(n, pi, pj)
    na = len(pairs)
    subsets = np.array(list(itertools.combinations(range(na), k)), np.int64)
    Xs = list(range(24, 32))
    RX = np.zeros(len(Xs), np.uint64)
    for i, X in enumerate(Xs):
        r = 0
        for v in range(X | 1, 2 * X, 2):
            r |= 1 << ((v - 1) // 2)
        RX[i] = np.uint64(r)
    t0 = time.time()
    br, bs, cover, perfect = scan(H, subsets, k, na, RX)
    j = int(np.argmax(br))
    print(f"order 7, k = {k}: {len(subsets)} arc sets x {1 << (na - k)} bases ({time.time()-t0:.1f}s)")
    print(f"   longest run of consecutive odd values inside one {k}-cube: {int(br.max())} (start {int(bs[j])}, arcs "
          f"{[pairs[a] for a in subsets[j]]}); perfect cubes (2^k consecutive odd values): {int(perfect.sum())}")
    for i, X in enumerate(Xs):
        print(f"   cubes containing all odd numbers of [{X}, {2*X}]: {int(cover[:, i].sum())}")

if __name__ == '__main__':
    main()
