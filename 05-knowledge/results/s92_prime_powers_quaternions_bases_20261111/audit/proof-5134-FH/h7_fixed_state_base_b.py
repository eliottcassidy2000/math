#!/usr/bin/env python3
"""Auditor proof-5134-FH.  THM-5134 (A) in base b, tested directly for fixed smaller states (own code).

For fixed (b, h1, s) and an odd modulus m prime to b, d = gcd(X1, y1, m):
   e(N_b) = N_b (s >= 1) or 2 N_b + b (s = 0) is claimed to be uniform on the multiples of d modulo m over
   {h2 > h1 : N_b <= X}, with a power saving:  #{e = c'} = [d | c'] (d/m) G + O((X / (b 2^s y1))^(1 - theta_m)).
Reported:  G,  the number of h2 with d not dividing e (must be 0),  and
   dev = max over the multiples c' of d of | #{e = c'} / ((d/m) G) - 1 |,   against the sampling scale sqrt(m/(d G)).
psi is computed from the closed sum q(h) = sum over one-bits i of ceil(floor(h/2^(i+1))/2), no table.
"""
import numpy as np
from math import gcd
from numba import njit


@njit(cache=True)
def q_sum(h):
    tot = 0
    i = 0
    while (h >> i) != 0:
        if (h >> i) & 1:
            tot += ((h >> (i + 1)) + 1) >> 1
        i += 1
    return tot


@njit(cache=True)
def scan(b, h1, s, X, m, cnt):
    y1 = 2 * h1 + 1
    X1 = 1 + b * (h1 - q_sum(h1))
    G = 0
    h2 = h1 + 1
    while True:
        y2 = 2 * h2 + 1
        if s == 0:
            seam = b * (y1 * y2 - 1) // 2
        else:
            seam = (b << (s - 1)) * y1 * y2
        if X1 * (1 + b * ((h2 + 1) // 2)) + seam > X:
            break
        N = X1 * (1 + b * (h2 - q_sum(h2))) + seam
        if N <= X:
            G += 1
            e = N if s >= 1 else 2 * N + b
            cnt[e % m] += 1
        h2 += 1
    return G


def psi(h):
    return h - q_sum(h)


def main():
    cases = [
        (3, 1, 1, [7, 25, 49, 35]),          # X1 = 4, y1 = 3: type 1
        (3, 3, 1, [7, 49, 35, 245]),         # X1 = y1 = 7: type 7 (balanced state of base 3)
        (3, 3, 0, [7, 49, 35]),
        (3, 45, 1, [91, 637, 49, 169]),      # X1 = y1 = 91 = 7 * 13
        (10, 1, 1, [7, 49, 9, 21]),          # X1 = 11, y1 = 3
        (10, 3, 1, [7, 49, 21, 147]),        # X1 = 21, y1 = 7: type 7, also 3 | X1 only
        (10, 3, 0, [7, 49, 21]),
        (5, 6, 1, [13, 169, 39]),            # X1 = 26, y1 = 13: type 13
        (5, 6, 0, [13, 169, 39]),
        (4, 2, 1, [5, 25, 15]),              # X1 = 9, y1 = 5: type 1 at 5
    ]
    print("THM-5134 (A) in base b for fixed smaller states (auditor proof-5134-FH)")
    for (b, h1, s, mods) in cases:
        X1, y1 = 1 + b * psi(h1), 2 * h1 + 1
        print(f"\nb = {b}, h1 = {h1} (X1 = {X1}, y1 = {y1}), s = {s}:")
        for m in mods:
            assert gcd(m, 2 * b) == 1
            d = gcd(gcd(X1, y1), m)
            row = []
            for K in (20, 24, 28, 31):
                X = 1 << K
                cnt = np.zeros(m, dtype=np.int64)
                G = scan(b, h1, s, X, m, cnt)
                if G == 0:
                    continue
                off = int(sum(int(cnt[c]) for c in range(m) if c % d != 0))
                exp = d * G / m
                dev = max(abs(int(cnt[c]) / exp - 1) for c in range(0, m, d))
                row.append(f"2^{K}: G = {G}, off-type {off}, dev = {dev:.4f} (noise {((m / d) / G) ** 0.5:.4f})")
            print(f"  m = {m:4d}, d = {d:3d}: " + "; ".join(row))


if __name__ == "__main__":
    main()
