#!/usr/bin/env python3
"""Independent check (verifier C) of two purely combinatorial steps of THM-5141 (C).

(R) Run-free lemma, on ALL ordered pairs h != h' of n-bit words (top bit 1), n <= NR (not only colliding pairs):
    with m = the longest run of equal bits in h or h', the windows satisfy
        |X|, |Xi| >= 2^(w-1-m) + 1   if w >= m + 2,       |X|, |Xi| >= 1 >= 2^(w-1-m)   if w <= m + 1
    (X = x - x', Xi = rev_w x - rev_w x').  Reports the minimum of |X| - 2^(w-1-m) over the pairs with w >= m + 2
    (the lemma claims >= 1) and an extremal pair.
(K) Carries, abstract system.  The proof of |s - t'| <= rho + 2 uses only:
        c X - P Xi = gamma, |gamma| <= 2^w, X, Xi odd, 2^(w-rho) <= |X|, |Xi| < 2^w,
        P in [2^(t'-1), 2^t'), c in (2^s, 2^(s+1)].
    For all (w, t', s) with w <= WK, t' <= TK, s <= SK the integer solutions of this system are enumerated and the
    extreme values of t' - s - rho and s - t' - rho are reported (the lemma: <= 2 and <= -1), to see whether the
    two bounds are sharp for the hypotheses that the proof uses.
"""
import sys
import time
import numpy as np
from numba import njit

NR = int(sys.argv[1]) if len(sys.argv) > 1 else 12
WK = 7
TK = 9
SK = 6


@njit(cache=False)
def bitlen(a):
    b = 0
    while a > 0:
        a >>= 1
        b += 1
    return b


@njit(cache=False)
def revw(x, w):
    r = 0
    for i in range(w):
        if (x >> i) & 1:
            r |= 1 << (w - 1 - i)
    return r


@njit(cache=False)
def maxrun(h, n):
    best = 1
    cur = 1
    prev = h & 1
    for i in range(1, n):
        b = (h >> i) & 1
        if b == prev:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 1
            prev = b
    return best


@njit(cache=False)
def runlemma(n, res):
    """res[0] pairs, res[1] violations, res[2] min over pairs with w >= m+2 of (min(|X|,|Xi|) - 2^(w-1-m)),
    res[3], res[4] an extremal pair, res[5] pairs with w >= m + 2, res[6] number of pairs attaining the minimum"""
    lo = 1 << (n - 1)
    hi = 1 << n
    mr = np.empty(hi, dtype=np.int64)
    for h in range(lo, hi):
        mr[h] = maxrun(h, n)
    res[2] = 1 << 60
    for h in range(lo, hi):
        for hp in range(lo, hi):
            if h == hp:
                continue
            d = h ^ hp
            top = bitlen(d) - 1
            s = 0
            while ((d >> s) & 1) == 0:
                s += 1
            w = top - s + 1
            x = (h >> s) & ((1 << w) - 1)
            xp = (hp >> s) & ((1 << w) - 1)
            aX = abs(x - xp)
            aXi = abs(revw(x, w) - revw(xp, w))
            m = max(mr[h], mr[hp])
            res[0] += 1
            mn = min(aX, aXi)
            if w >= m + 2:
                res[5] += 1
                slack = mn - (1 << (w - 1 - m))
                if slack < 1:
                    res[1] += 1
                if slack < res[2]:
                    res[2] = slack
                    res[3] = h
                    res[4] = hp
                    res[6] = 1
                elif slack == res[2]:
                    res[6] += 1
            else:
                if mn < 1 or (aX & 1) == 0 or (aXi & 1) == 0:
                    res[1] += 1


@njit(cache=False)
def carries_abstract(w, tp, s, res):
    """res[0] solutions, res[1] max (tp - s - rho), res[2] max (s - tp - rho) (over solutions, with rho the least
    integer with 2^(w-rho) <= |X|, |Xi|), res[3..7] / res[8..12] examples (X, Xi, P, c, rho)"""
    for X in range(1, 1 << w, 2):                 # X > 0 without loss of generality (change the signs of X, Xi, gamma)
        for Xi in range(-(1 << w) + 1, 1 << w, 2):
            aXi = abs(Xi)
            mn = min(X, aXi)
            rho = w - (bitlen(mn) - 1)
            for P in range(1 << (tp - 1), 1 << tp):
                for c in range((1 << s) + 1, (1 << (s + 1)) + 1):
                    gam = c * X - P * Xi
                    if abs(gam) > (1 << w):
                        continue
                    res[0] += 1
                    a = tp - s - rho
                    b = s - tp - rho
                    if a > res[1]:
                        res[1] = a
                        res[3] = X
                        res[4] = Xi
                        res[5] = P
                        res[6] = c
                        res[7] = rho
                    if b > res[2]:
                        res[2] = b
                        res[8] = X
                        res[9] = Xi
                        res[10] = P
                        res[11] = c
                        res[12] = rho


def main():
    t0 = time.time()
    print("(R) run-free lemma on all ordered pairs of n-bit words", flush=True)
    for n in range(2, NR + 1):
        res = np.zeros(8, dtype=np.int64)
        runlemma(n, res)
        ex = ""
        if res[5] > 0:
            ex = "; min of min(|X|,|Xi|) - 2^(w-1-m) over the %d pairs with w >= m + 2: %d (lemma: >= 1), attained by %d pairs, e.g. h = %s, h' = %s" % (
                res[5], res[2], res[6], format(int(res[3]), "0%db" % n), format(int(res[4]), "0%db" % n))
        print("n = %2d: ordered pairs %9d, violations %d%s" % (n, res[0], res[1], ex), flush=True)
    print()
    print("(K) abstract carries system, w <= %d, t' <= %d, s <= %d" % (WK, TK, SK), flush=True)
    gmaxa = -99
    gmaxb = -99
    exa = None
    exb = None
    tot = 0
    for w in range(1, WK + 1):
        for tp in range(1, TK + 1):
            for s in range(0, SK + 1):
                res = np.full(16, -99, dtype=np.int64)
                res[0] = 0
                carries_abstract(w, tp, s, res)
                tot += int(res[0])
                if res[0] > 0:
                    if res[1] > gmaxa:
                        gmaxa = int(res[1])
                        exa = (w, tp, s) + tuple(int(v) for v in res[3:8])
                    if res[2] > gmaxb:
                        gmaxb = int(res[2])
                        exb = (w, tp, s) + tuple(int(v) for v in res[8:13])
    print("solutions examined: %d" % tot)
    print("max (t' - s - rho) = %d (lemma: <= 2); example (w, t', s, X, Xi, P, c, rho) = %s" % (gmaxa, exa))
    if exa:
        w, tp, s, X, Xi, P, c, rho = exa
        print("    check: c X - P Xi = %d, 2^w = %d, 2^(w-rho) = %s" % (c * X - P * Xi, 1 << w, 2.0 ** (w - rho)))
    print("max (s - t' - rho) = %d (lemma: <= -1); example (w, t', s, X, Xi, P, c, rho) = %s" % (gmaxb, exb))
    if exb:
        w, tp, s, X, Xi, P, c, rho = exb
        print("    check: c X - P Xi = %d, 2^w = %d, 2^(w-rho) = %s" % (c * X - P * Xi, 1 << w, 2.0 ** (w - rho)))
    print("time %.1f s" % (time.time() - t0))


main()
