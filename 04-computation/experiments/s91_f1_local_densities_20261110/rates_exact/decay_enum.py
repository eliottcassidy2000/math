#!/usr/bin/env python3
"""S91 rates_exact (B): dyadic-block sums of the limiting F1 weights over all rooted states and over the p-divisible ones.

x(h) = y(h) - 2 q(h), y(h) = 2h + 1, rho = x / y.  Limiting weights of a first state h on a window [X, RX), per unit of
(R - 1) X (window_split.py: W(h, s) -> (B - A) omega):
        sheet 0 (s = 0):   omega_0(h) = (1/x) ln( 2 (1 + rho) / (2 + rho) )
        sheet 1 (s >= 1):  omega_1(h) = (1/x) ln( 1 + rho / 2 )
For every dyadic block n (2^n <= h < 2^(n+1), n < NB) this script accumulates
        S_all[n, sheet]      = sum of omega_sheet(h) over the block,
        S_div[p, n, sheet]   = the same over the h with p | gcd(x(h), y(h)),
        count[p, n]          = number of such h,
for the primes of PRIMES, by direct enumeration of all h < 2^NB (no big table: q(P 2^r + L) = q(P) + q(L) +
P rbar(L) + L_top ceil(P/2), rbar(L) = sum over one-bits i <= r - 2 of L of 2^(r-2-i)).

Checks (exit 1 on failure): the prefix-suffix formula for q against the recursion for all h < 2^20; the block sums
S_all[n] -> (m_sheet / 2) ln 2 with m_0 = 0.3230778345, m_1 = 0.4250278185 (constants.py); h = 2^j is divisible for
every prime factor of 2^(j+1) + 1.

Usage: python3 decay_enum.py [NB] [OUT.npz]      (defaults 26, decay_enum_26.npz; NB = 30 takes a few minutes)
"""
import math
import sys
import time

import numpy as np
from numba import njit, prange, set_num_threads

PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 257]
RLOW = 12
M0 = 0.3230778345
M1 = 0.4250278185


@njit(cache=False)
def q_table(size):
    q = np.zeros(size, np.int64)
    for h in range(1, size):
        half = h >> 1
        q[h] = q[half] + (((half + 1) >> 1) if (h & 1) else 0)
    return q


@njit(cache=False)
def low_tables(r):
    size = 1 << r
    rbar = np.zeros(size, np.int64)
    top = np.zeros(size, np.int64)
    for L in range(size):
        acc = 0
        for i in range(r - 1):
            if (L >> i) & 1:
                acc += 1 << (r - 2 - i)
        rbar[L] = acc
        top[L] = (L >> (r - 1)) & 1
    return rbar, top


@njit(parallel=True, cache=False)
def enumerate_blocks(nb, r, qP, qL, rbar, top, primes, nchunks):
    npr = primes.shape[0]
    nP = 1 << (nb - r)
    nL = 1 << r
    S_all = np.zeros((nchunks, nb, 2))
    S_div = np.zeros((nchunks, npr, nb, 2))
    cnt = np.zeros((nchunks, npr, nb), np.int64)
    step = nP // nchunks + 1
    for c in prange(nchunks):
        for P in range(c * step, min(nP, (c + 1) * step)):
            nbase = -1
            if P > 0:
                nbase = r
                t = P
                while t > 1:
                    t >>= 1
                    nbase += 1
            cP = (P + 1) >> 1
            for L in range(nL):
                h = (P << r) + L
                if h == 0:
                    continue
                if P > 0:
                    n = nbase
                else:
                    n = 0
                    t = L
                    while t > 1:
                        t >>= 1
                        n += 1
                q = qP[P] + qL[L] + P * rbar[L] + top[L] * cP
                y = 2 * h + 1
                x = y - 2 * q
                w0 = math.log((x + y) / (0.5 * x + y)) / x
                w1 = math.log1p(0.5 * x / y) / x
                S_all[c, n, 0] += w0
                S_all[c, n, 1] += w1
                for j in range(npr):
                    p = primes[j]
                    if y % p == 0 and x % p == 0:
                        S_div[c, j, n, 0] += w0
                        S_div[c, j, n, 1] += w1
                        cnt[c, j, n] += 1
    return S_all.sum(axis=0), S_div.sum(axis=0), cnt.sum(axis=0)


def main():
    nb = int(sys.argv[1]) if len(sys.argv) > 1 else 26
    out = sys.argv[2] if len(sys.argv) > 2 else f"decay_enum_{nb}.npz"
    set_num_threads(4)
    t0 = time.time()
    r = RLOW
    qfull = q_table(1 << max(20, nb - r))
    qL = qfull[: 1 << r].copy()
    qP = qfull[: 1 << (nb - r)].copy()
    rbar, top = low_tables(r)
    ok = True
    # check of the prefix-suffix formula on all h < 2^20
    hh = np.arange(1 << 20, dtype=np.int64)
    P, L = hh >> r, hh & ((1 << r) - 1)
    ok &= bool(np.array_equal(qfull[: 1 << 20], qfull[P] + qL[L] + P * rbar[L] + top[L] * ((P + 1) >> 1)))
    print("prefix-suffix formula for q on all h < 2^20:", ok)
    primes = np.array(PRIMES, dtype=np.int64)
    S_all, S_div, cnt = enumerate_blocks(nb, r, qP, qL, rbar, top, primes, 256)
    # h = 0: x = y = 1, rho = 1 (never divisible)
    w00 = math.log(4.0 / 3.0)
    w01 = math.log(1.5)
    np.savez(out, nb=nb, primes=primes, S_all=S_all, S_div=S_div, count=cnt, w_h0=np.array([w00, w01]))
    print(f"enumerated all h < 2^{nb} in {time.time() - t0:.0f} s; written {out}")
    print("block sums of the limiting weights against (m_sheet / 2) ln 2:")
    for n in list(range(0, nb, 4)) + [nb - 1]:
        print(f"   n = {n:2d}: sheet 0 {S_all[n, 0]:.9f} ({M0 / 2 * math.log(2):.9f})   sheet 1 {S_all[n, 1]:.9f} ({M1 / 2 * math.log(2):.9f})")
    ok &= abs(S_all[nb - 1, 0] - M0 / 2 * math.log(2)) < 1e-6 and abs(S_all[nb - 1, 1] - M1 / 2 * math.log(2)) < 1e-6
    print("number of p-divisible states in the last block times p^2 / 2^n:")
    print("   " + "  ".join(f"{p}: {cnt[j, nb - 1] * p * p / 2.0 ** (nb - 1):.4f}" for j, p in enumerate(PRIMES)))
    # h = 2^j is divisible by every prime factor of 2^(j+1) + 1: the first divisible state is at most 2^(mm - 1)
    for j, p in enumerate(PRIMES):
        mm = next((e for e in range(1, 2 * p) if pow(2, e, p) == p - 1), 0)
        if mm and mm - 1 < nb:
            ok &= int(cnt[j, : mm].sum()) >= 1
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
