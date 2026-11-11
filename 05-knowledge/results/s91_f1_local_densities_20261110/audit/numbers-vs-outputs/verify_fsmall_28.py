#!/usr/bin/env python3
"""Auditor numbers-vs-outputs: exact #T_sigma(X) and D_sigma(X) (triples whose smaller state is p-divisible) at
X = 2^26 and 2^28 by an own enumeration, giving (F_small - 1/p^2) #T / (2X) to full precision; compared with
first_order_exact_28.out (3 decimals of F_small p^2) and the 4-decimal values quoted in THM-5128 (D4)."""
import sys

import numpy as np
from numba import njit

PR = np.array([3, 5, 7, 11, 13, 17, 31], dtype=np.int64)


@njit(cache=False)
def build(size, primes):
    q = np.zeros(size, np.int32)
    mask = np.zeros(size, np.int32)
    for h in range(size):
        if h > 0:
            half = h >> 1
            q[h] = q[half] + (((half + 1) >> 1) if (h & 1) else 0)
        y = 2 * h + 1
        x = y - 2 * q[h]
        m = 0
        for j in range(primes.shape[0]):
            if y % primes[j] == 0 and x % primes[j] == 0:
                m |= (1 << j)
        mask[h] = m
    return q, mask


@njit(cache=False)
def count(q, mask, X, K, nmask):
    cnt = np.zeros((2, nmask), np.int64)
    H = q.shape[0]
    for s in range(K + 1):
        sheet = 0 if s == 0 else 1
        d0 = 1 if s == 0 else 0
        for h1 in range(H):
            y1 = 2 * h1 + 1
            Y = y1 << s
            if Y > X:
                break
            x1 = y1 - 2 * q[h1]
            y2max = (2 * (X + d0)) // (2 * Y + x1)
            if y2max < 1:
                continue
            h2max = (y2max - 1) // 2 + 1
            if h2max >= H:
                h2max = H - 1
            m1 = mask[h1]
            for h2 in range(0, h2max + 1):
                y2 = 2 * h2 + 1
                x2 = y2 - 2 * q[h2]
                if x1 * x2 + Y * y2 - d0 <= X:
                    if h1 <= h2:
                        cnt[sheet, m1] += 1
                    else:
                        cnt[sheet, mask[h2]] += 1
    return cnt


for K in (26, 28):
    X = 1 << K
    q, mask = build(X // 2 + 2, PR)
    cnt = count(q, mask, X, K, 1 << len(PR))
    tot = cnt.sum(axis=1)
    print(f"X = 2^{K}: triples sheet 0: {int(tot[0])}, sheet 1: {int(tot[1])}")
    for j, p in enumerate(PR.tolist()):
        row = []
        for sheet in (0, 1):
            d = int(sum(cnt[sheet, m] for m in range(1 << len(PR)) if (m >> j) & 1))
            fs = d / tot[sheet]
            row.append(f"sheet {sheet}: F_small p^2 = {fs * p * p:.5f}, (F_small - 1/p^2) #T/(2X) = {(fs - 1.0 / (p * p)) * tot[sheet] / (2.0 * X):+.6f}")
        print(f"   p = {p:2d}: " + ";  ".join(row))
    del q, mask
