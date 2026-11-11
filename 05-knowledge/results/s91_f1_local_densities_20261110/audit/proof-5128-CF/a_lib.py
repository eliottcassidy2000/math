#!/usr/bin/env python3
"""Auditor proof-5128-CF: own implementation of the F1 census (no code shared with the session's scripts).

psi(0) = 0, psi(2h) = psi(h) + h, psi(2h+1) = psi(h) + floor(h/2) + 1;  x = 2 psi + 1, y = 2h + 1.
N(h1, h2, s) = x1 x2 + 2^s y1 y2 - [s = 0], ordered triples.
"""
import numpy as np
from numba import njit


@njit(cache=True)
def psi_table(H):
    psi = np.zeros(H, dtype=np.int64)
    for h in range(1, H):
        k = h >> 1
        if h & 1:
            psi[h] = psi[k] + (k >> 1) + 1
        else:
            psi[h] = psi[k] + k
    return psi


def psi_slow(h):
    """Plain recursion (python ints), used to cross-check the table."""
    if h == 0:
        return 0
    k = h >> 1
    if h & 1:
        return psi_slow(k) + (k >> 1) + 1
    return psi_slow(k) + k


@njit(cache=True)
def kbin(N, kmin):
    """Smallest k >= kmin with N <= 2^k."""
    k = kmin
    while (np.int64(1) << k) < N:
        k += 1
    return k


@njit(cache=True)
def count_by_padding(K, x, smax_rec):
    """T[s, k] = number of ORDERED pairs (h1, h2) with N(h1, h2, s) <= 2^k, for s <= smax_rec (larger s pooled in
    row smax_rec + 1), k = 0..K.  Also diag[s, k] (h1 = h2)."""
    X = np.int64(1) << K
    T = np.zeros((smax_rec + 2, K + 1), dtype=np.int64)
    diag = np.zeros((smax_rec + 2, K + 1), dtype=np.int64)
    for s in range(0, K + 1):
        B = X + 1 if s == 0 else X
        row = s if s <= smax_rec else smax_rec + 1
        h1 = 0
        while True:
            y1 = 2 * h1 + 1
            if (y1 * y1) << s >= B:
                break
            x1 = np.int64(x[h1])
            by = y1 << s
            h2max = ((B - 1) // by - 1) // 2
            # diagonal
            M = x1 * x1 + by * y1
            if M <= B:
                N = M - 1 if s == 0 else M
                diag[row, kbin(N, 0)] += 1
                T[row, kbin(N, 0)] += 1
            for h2 in range(h1 + 1, h2max + 1):
                M = x1 * np.int64(x[h2]) + by * (2 * h2 + 1)
                if M <= B:
                    N = M - 1 if s == 0 else M
                    T[row, kbin(N, 0)] += 2
            h1 += 1
    # cumulative in k
    for r in range(smax_rec + 2):
        for k in range(1, K + 1):
            T[r, k] += T[r, k - 1]
            diag[r, k] += diag[r, k - 1]
    return T, diag


@njit(cache=True)
def census_mod(K, kmin, x, m, cls, ncls):
    """hist[sheet, k - kmin, cls_small, cls_large, N mod m] over UNORDERED pairs h1 < h2 (each stands for two ordered
    triples), cumulative in k afterwards is done by the caller; diag[sheet, k - kmin, cls, N mod m] for h1 = h2."""
    X = np.int64(1) << K
    nk = K - kmin + 1
    hist = np.zeros((2, nk, ncls, ncls, m), dtype=np.int64)
    diag = np.zeros((2, nk, ncls, m), dtype=np.int64)
    for s in range(0, K + 1):
        B = X + 1 if s == 0 else X
        sheet = 0 if s == 0 else 1
        h1 = 0
        while True:
            y1 = 2 * h1 + 1
            if (y1 * y1) << s >= B:
                break
            x1 = np.int64(x[h1])
            c1 = cls[h1]
            by = y1 << s
            h2max = ((B - 1) // by - 1) // 2
            M = x1 * x1 + by * y1
            if M <= B:
                N = M - 1 if s == 0 else M
                diag[sheet, kbin(N, kmin) - kmin, c1, N % m] += 1
            for h2 in range(h1 + 1, h2max + 1):
                M = x1 * np.int64(x[h2]) + by * (2 * h2 + 1)
                if M <= B:
                    N = M - 1 if s == 0 else M
                    hist[sheet, kbin(N, kmin) - kmin, c1, cls[h2], N % m] += 1
            h1 += 1
    return hist, diag
