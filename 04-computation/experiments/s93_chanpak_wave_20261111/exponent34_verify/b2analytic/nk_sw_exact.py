#!/usr/bin/env python3
"""
Exact computation of the height histogram of the odd window pairs, of N(K) and of the gcd sum S_w of THM-5141,
and comparison with every inequality of the paragraph "The sum" of the proof of (B2)   (verifier: b2analytic).

Definitions (theorem file): for w-bit words x, x':  X = x - x', Xi = rev_w x - rev_w x'; the ordered pair is odd if
X and Xi are odd; height M = max(|X|, |Xi|)/gcd(X, Xi); N(K) = #{odd pairs : M <= K}; S_w = sum over odd pairs 1/M.

Method: histogram hist[M] over the pairs with top bit (x) = 1, top bit (x') = 0 and different bottom bits (for
w >= 2 these are exactly half of the odd pairs; the other half is (x', x), with the same height), integer gcd.
S_w is then an exact rational (fractions.Fraction).  For w <= WBRUTE the histogram is recomputed by a plain-Python
brute force straight from the definition (all 4^w ordered pairs, oddness tested on X and Xi) and compared.

Inequalities tested (s1 = ceil(w/2), s2 = floor(w/2)):
  (T1) N(K)/K <= C 8^s1 for 1 <= K <= 2^s2            -> report max over ALL integers K in [1, 2^s2] of N(K)/(K 8^s1)
  (T2) N(K) <= 3 * 6^s1 2^s1 (3 + 2 ln K), 1 <= K <= 2^s2 (explicit, from (D2))  -> max ratio must be <= 1
  (T3) N(K) <= C 6^s1 [K (1+ln K)^2 + K (1+ln K) 2^(s1/3) + K^(5/2) 2^(-s2/2)], all K >= 1  -> report max ratio (= needed C)
  (T4) N(K)/K <= C 6^s1 [(1+s1)^2 + (1+s1) 2^(s1/3) + (1+s1)^(3/5) 2^(2 s1/5)], K <= 2^s2 -> report max ratio
  (T5) S_w <= sum_{i>=0} 2^(-i) N(2^(i+1))   (inequality (2))  -> exact rational comparison
  (T6) S_w / (w 2^(3w/2)),  S_w / 6^(w/2)
usage: python nk_sw_exact.py [wmin] [wmax]
"""
import math
import sys
import time
from fractions import Fraction

import numpy as np
from numba import njit, prange

WBRUTE = 8


@njit(cache=False)
def gcd_int(a, b):
    while b:
        a, b = b, a % b
    return a


@njit(parallel=True, cache=False)
def height_hist(w, rev):
    n = 1 << w
    half = n >> 1
    nch = 64
    hist = np.zeros((nch, n), dtype=np.int64)
    for c in prange(nch):
        for x in range(half + c, n, nch):
            rx = rev[x]
            b0 = x & 1
            for xp in range(1 - b0, half, 2):      # bottom bit of x' differs from that of x
                X = x - xp
                Xi = rx - rev[xp]
                aXi = Xi if Xi >= 0 else -Xi
                g = gcd_int(X, aXi)
                M = (X if X >= aXi else aXi) // g
                hist[c, M] += 1
    tot = np.zeros(n, dtype=np.int64)
    for c in range(nch):
        for i in range(n):
            tot[i] += hist[c, i]
    return tot


def rev_table(w):
    n = 1 << w
    rev = np.zeros(n, dtype=np.int64)
    for x in range(n):
        r = 0
        for i in range(w):
            if (x >> i) & 1:
                r |= 1 << (w - 1 - i)
        rev[x] = r
    return rev


def brute_hist(w):
    n = 1 << w
    rev = rev_table(w)
    hist = [0] * n
    for x in range(n):
        for xp in range(n):
            X = x - xp
            Xi = int(rev[x]) - int(rev[xp])
            if X % 2 == 0 or Xi % 2 == 0:
                continue
            g = math.gcd(X, Xi)
            hist[max(abs(X), abs(Xi)) // g] += 1
    return hist


def main():
    wmin = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    wmax = int(sys.argv[2]) if len(sys.argv) > 2 else 14
    print("w  s1 s2 | #odd pairs | S_w (exact rational, float) | S_w/6^(w/2) | S_w/(w 2^(3w/2)) | hist[1]/6^(w/2) |"
          " T1 max N(K)/(K 8^s1) [argmax K] | T2 max [K] | T3 max [K] | T4 max [K] | T5 RHS/S_w")
    summary = []
    for w in range(wmin, wmax + 1):
        t0 = time.time()
        rev = rev_table(w)
        half_hist = height_hist(w, rev)
        hist = [2 * int(v) for v in half_hist]
        n = 1 << w
        npairs = sum(hist)
        assert npairs == 4 ** w // 4, (w, npairs)
        if w <= WBRUTE:
            bh = brute_hist(w)
            assert bh == hist, f"brute-force mismatch at w = {w}"
        S = Fraction(0)
        for M in range(1, n):
            if hist[M]:
                S += Fraction(hist[M], M)
        s1 = (w + 1) // 2
        s2 = w // 2
        # cumulative N(K) for integer K
        Ncum = [0] * n
        acc = 0
        for M in range(1, n):
            acc += hist[M]
            Ncum[M] = acc

        def N(K):
            K = int(K)
            if K >= n:
                return acc
            return Ncum[K] if K >= 1 else 0

        t1 = (0.0, 0)
        t2 = (0.0, 0)
        t4 = (0.0, 0)
        for K in range(1, (1 << s2) + 1):
            v = N(K)
            r = v / (K * 8.0 ** s1)
            if r > t1[0]:
                t1 = (r, K)
            r = v / (3.0 * 6.0 ** s1 * 2.0 ** s1 * (3.0 + 2.0 * math.log(K)))
            if r > t2[0]:
                t2 = (r, K)
            r = (v / K) / (6.0 ** s1 * ((1 + s1) ** 2 + (1 + s1) * 2.0 ** (s1 / 3.0)
                                        + (1 + s1) ** 0.6 * 2.0 ** (0.4 * s1)))
            if r > t4[0]:
                t4 = (r, K)
        t3 = (0.0, 0)
        for K in range(1, n):
            l = 1.0 + math.log(K)
            den = 6.0 ** s1 * (K * l * l + K * l * 2.0 ** (s1 / 3.0) + K ** 2.5 * 2.0 ** (-s2 / 2.0))
            r = N(K) / den
            if r > t3[0]:
                t3 = (r, K)
        # (2): S_w <= sum_i 2^-i N(2^(i+1))
        rhs = Fraction(0)
        i = 0
        while True:
            rhs += Fraction(N(2 ** (i + 1)), 2 ** i)
            if 2 ** (i + 1) >= n:
                # all later terms: N = total
                rhs += Fraction(acc, 2 ** i)      # sum_{i' > i} 2^-i' = 2^-i
                break
            i += 1
        assert S <= rhs
        Sf = float(S)
        print(f"{w:2d} {s1:2d} {s2:2d} | {npairs} | {S.numerator}/{S.denominator if S.denominator < 10**12 else '...'} = {Sf:.6f}"
              f" | {Sf / 6.0 ** (w / 2.0):.4f} | {Sf / (w * 2.0 ** (1.5 * w)):.6f} | {hist[1] / 6.0 ** (w / 2.0):.4f} |"
              f" {t1[0]:.5f} [{t1[1]}] | {t2[0]:.5f} [{t2[1]}] | {t3[0]:.5f} [{t3[1]}] | {t4[0]:.5f} [{t4[1]}] |"
              f" {float(rhs) / Sf:.3f}   ({time.time() - t0:.1f}s)")
        # dyadic profile of N(K)/K
        prof = []
        for i in range(0, w + 1):
            K = 2 ** i
            prof.append(f"{N(K) / K / 6.0 ** s1:.3f}")
        print(f"      N(2^i)/(2^i 6^s1), i = 0..{w}: " + " ".join(prof))
        sys.stdout.flush()
        summary.append((w, Sf, Sf / 6.0 ** (w / 2.0), t1[0], t2[0], t3[0], t4[0]))
    print("summary: max over w of T1 = %.5f, T2 = %.5f, T3 = %.5f, T4 = %.5f; S_w/6^(w/2) range [%.4f, %.4f] for w >= 4"
          % (max(s[3] for s in summary), max(s[4] for s in summary), max(s[5] for s in summary),
             max(s[6] for s in summary),
             min([s[2] for s in summary if s[0] >= 4] or [0]), max([s[2] for s in summary if s[0] >= 4] or [0])))


if __name__ == "__main__":
    main()
