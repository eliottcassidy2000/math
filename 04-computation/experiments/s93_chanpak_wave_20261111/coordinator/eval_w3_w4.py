#!/usr/bin/env python3
"""S93, evaluation of the pre-registered W3 and W4 (PREREGISTRATION.md, pushed in commit 3a675e5b7a before this run).

W4: G(tau) = sum over ordered pairs (T, T') of tau-bit strings with different last bits of
    gcd(alpha, beta) / max(|alpha|, |beta|), alpha = rev(T) - rev(T'), beta = T - T'.
    Registered: G(13)/6^(13/2) and G(14)/6^7 lie in [1.2, 2.2].  (In sample: tau <= 12.)
W3: one-apex words of n bits (top bit 1), x(h) = 2h + 1 - 2 q(h): colliding pairs per word with common prefix of
    length >= t.  Registered for n = 28 and n = 30: within [0.118, 0.185] (t = 4), [0.015, 0.024] (t = 7),
    [0.0037, 0.0058] (t = 9), and the share of distinct (prefix, value) classes at t = ceil(n/4) at least 0.95.
Usage: python3 eval_w3_w4.py w4 | w3 N
"""
import math
import sys

import numpy as np
from numba import njit, prange


@njit(cache=True)
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


@njit(parallel=True, cache=True)
def gsum(tau):
    size = 1 << tau
    rev = np.zeros(size, np.int64)
    for t in range(size):
        r = 0
        for i in range(tau):
            if (t >> i) & 1:
                r |= 1 << (tau - 1 - i)
        rev[t] = r
    acc = np.zeros(size, np.float64)
    for t in prange(size):
        s = 0.0
        for u in range((t & 1) ^ 1, size, 2):       # different last bit
            be = t - u
            al = rev[t] - rev[u]
            if be < 0:
                be = -be
            if al < 0:
                al = -al
            g = gcd(al, be)
            m = al if al > be else be
            s += g / m
        acc[t] = s
    return acc.sum()


@njit(parallel=True, cache=True)
def keys_for(n, out):
    """out[i] = x(h) << 30 | (h - 2^(n-1)) for h = 2^(n-1) + i."""
    lo = 1 << (n - 1)
    for i in prange(out.shape[0]):
        h = lo + i
        q = 0
        for b in range(n):
            if (h >> b) & 1:
                pre = h >> (b + 1)
                q += (pre + 1) >> 1
        x = 2 * h + 1 - 2 * q
        out[i] = (x << 30) | i


@njit(cache=True)
def scan(keys, n, t0, hist):
    """hist[cp] = number of colliding unordered pairs with common prefix length exactly cp; returns the number of
    distinct (prefix of length t0, value) classes."""
    N = keys.shape[0]
    mask = (1 << 30) - 1
    lo = 1 << (n - 1)
    classes = 0
    i = 0
    while i < N:
        x = keys[i] >> 30
        j = i + 1
        while j < N and (keys[j] >> 30) == x:
            j += 1
        # group [i, j)
        for a in range(i, j):
            ha = lo + (keys[a] & mask)
            newclass = True
            for b in range(i, a):
                hb = lo + (keys[b] & mask)
                d = ha ^ hb
                bl = 0
                while d:
                    bl += 1
                    d >>= 1
                cp = n - bl
                hist[cp] += 1
                if cp >= t0:
                    newclass = False
            if newclass:
                classes += 1
        i = j
    return classes


def main():
    mode = sys.argv[1]
    if mode == "w4":
        ok = True
        for tau in (10, 11, 12, 13, 14):
            G = gsum(tau)
            ratio = G / 6 ** (tau / 2)
            tag = ""
            if tau >= 13:
                inside = 1.2 <= ratio <= 2.2
                ok &= inside
                tag = f"  registered [1.2, 2.2]: {'inside' if inside else 'OUTSIDE'}"
            print(f"tau = {tau}: G = {G:.3f} = {G / 2 ** tau:.3f} * 2^tau; G / 6^(tau/2) = {ratio:.4f}{tag}", flush=True)
        print("W4:", "CONFIRMED" if ok else "REFUTED")
    else:
        n = int(sys.argv[2])
        N = 1 << (n - 1)
        keys = np.empty(N, np.int64)
        keys_for(n, keys)
        keys.sort()
        t0 = math.ceil(n / 4)
        hist = np.zeros(n + 1, np.int64)
        classes = scan(keys, n, t0, hist)
        tot_pairs = int(hist.sum())
        print(f"n = {n}: {N} words; colliding pairs per word: {tot_pairs / N:.4f}")
        bounds = {4: (0.118, 0.185), 7: (0.015, 0.024), 9: (0.0037, 0.0058)}
        ok = True
        for t in (1, 3, 4, 5, 6, 7, 8, 9, t0):
            v = int(hist[t:].sum()) / N
            tag = ""
            if t in bounds and n in (28, 30):
                inside = bounds[t][0] <= v <= bounds[t][1]
                ok &= inside
                tag = f"  registered {list(bounds[t])}: {'inside' if inside else 'OUTSIDE'}"
            print(f"  common prefix >= {t}: {v:.5f} pairs per word{tag}")
        share = classes / N
        if n in (28, 30):
            ok &= share >= 0.95
        print(f"  distinct (prefix of length {t0}, value) classes: {classes} of {N} ({share:.4f}; registered at least 0.95 for n = 28, 30)")
        if n in (28, 30):
            print(f"W3 at n = {n}:", "CONFIRMED" if ok else "REFUTED")


if __name__ == "__main__":
    main()
