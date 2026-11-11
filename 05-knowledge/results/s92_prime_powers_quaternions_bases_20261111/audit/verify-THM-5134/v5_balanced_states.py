#!/usr/bin/env python3
"""Verifier of THM-5134 (G), balanced states: x_b(h) = y(h), i.e. b psi(h) = 2h.

Own code. q(h) is computed for every h < 2^30 from the closed sum of THM-5128 (B3),
q(h) = sum over the one-bits i of h of ceil(floor(h / 2^(i+1)) / 2), without any table; psi = h - q.
Base 3: b psi = 2h is h = 3 q(h). Bases >= 4: needs psi(h) <= h/2, excluded by 2 psi(h) > h.
"""
import sys

import numpy as np
from numba import njit, prange


@njit(cache=False)
def q_closed(h):
    s = 0
    i = 0
    t = h
    while t:
        if t & 1:
            H = h >> (i + 1)
            s += (H + 1) >> 1
        t >>= 1
        i += 1
    return s


@njit(parallel=True, cache=False)
def scan(lo, hi, out3, minslack):
    nb = (hi - lo) >> 20
    for blk in prange(nb):
        a = lo + (blk << 20)
        cnt = 0
        mn = 1 << 60
        for h in range(a, a + (1 << 20)):
            qq = q_closed(h)
            if h == 3 * qq:
                out3[blk, cnt] = h
                cnt += 1
            if h >= 1:
                sl = 2 * (h - qq) - h
                if sl < mn:
                    mn = sl
        minslack[blk] = mn


def psi_rec(h):
    if h == 0:
        return 0
    k = h >> 1
    return psi_rec(k) + (k if h % 2 == 0 else (k >> 1) + 1)


def main():
    ok = True
    # the closed sum against the recursion
    good = all(q_closed(h) == h - psi_rec(h) for h in list(range(5000)) + [659903829, 13456725, (1 << 30) - 1])
    print("closed sum q(h) = h - psi(h) (recursion) on 5003 values: %s" % good)
    ok &= good
    out3 = np.full((1 << 10, 8), -1, dtype=np.int64)
    mins = np.zeros(1 << 10, dtype=np.int64)
    scan(0, 1 << 30, out3, mins)
    sol = sorted(int(v) for v in out3.ravel() if v >= 0)
    print("balanced states of base 3 below 2^30: %d: %s" % (len(sol), sol))
    for K in (24, 26, 28, 30):
        print("   below 2^%d: %d" % (K, sum(1 for v in sol if v < (1 << K))))
    want = [0, 3, 45, 117, 309, 9045, 62805, 111957, 382293, 2839893, 4412757, 13456725, 659903829]
    ok &= sol == want
    h = 659903829
    print("h = 659,903,829: psi = %d, 1 + 3 psi = %d, 2h + 1 = %d" % (psi_rec(h), 1 + 3 * psi_rec(h), 2 * h + 1))
    ok &= 1 + 3 * psi_rec(h) == 2 * h + 1
    print("min of 2 psi(h) - h over 1 <= h < 2^30: %d (so psi(h) > h/2, and no balanced state h >= 1 in any base "
          "b >= 4 below 2^30)" % int(mins.min()))
    ok &= int(mins.min()) >= 1
    print("12 below 2^24, 12 below 2^28, 13 below 2^30, the thirteenth is 659,903,829: %s" % ok)
    print("ALL CHECKS PASSED" if ok else "FAILURE")
    sys.exit(0 if ok else 1)


main()
