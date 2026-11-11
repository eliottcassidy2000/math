#!/usr/bin/env python3
"""Auditor bridges-note, script A8 (independent code): COMPLETE sheet-0 census in base beta,
N_beta = X1 X2 + beta (y1 y2 - 1)/2 <= 2^22 over ordered triples (h1, h2, s = 0), X = 1 + beta psi(h), y = 2h + 1,
against bridges_exact.out S10 (whose state table stops at h < 2^12).

Also reported: the same census restricted to h1, h2 < 4096 (the session's box), to identify its numbers.
"""
import numpy as np
from numba import njit


@njit(cache=True)
def census(psi, beta, X, primes, cnt, bad, hmax):
    H = min(psi.shape[0], hmax)
    tot = 0
    for h1 in range(H):
        y1 = 2 * h1 + 1
        X1 = 1 + beta * psi[h1]
        # N >= X1 * 1 + beta (y1 y2 - 1)/2 with y2 >= y1
        if X1 + beta * (y1 * y1 - 1) // 2 > X:
            break
        for h2 in range(h1, H):
            y2 = 2 * h2 + 1
            if X1 + beta * (y1 * y2 - 1) // 2 > X:
                break
            N = X1 * (1 + beta * psi[h2]) + beta * (y1 * y2 - 1) // 2
            if N > X:
                continue
            w = 1 if h2 == h1 else 2
            tot += w
            for j in range(primes.shape[0]):
                p = primes[j]
                t = 1 if (X1 % p == 0 and y1 % p == 0) else 0
                cnt[j, t, N % p] += w
                if t == 1 and (2 * N + beta) % p != 0:
                    bad[j] += 1
    return tot


def run(psi, beta, X, hmax, label):
    primes = np.array([p for p in (3, 5, 7, 11, 13) if beta % p], np.int64)
    cnt = np.zeros((len(primes), 2, 16), np.int64)
    bad = np.zeros(len(primes), np.int64)
    tot = census(psi, beta, X, primes, cnt, bad, hmax)
    cells = []
    for j, p in enumerate(primes):
        p = int(p)
        c0 = (-beta * pow(2, -1, p)) % p
        D1, Dp = int(cnt[j, 0].sum()), int(cnt[j, 1].sum())
        c = cnt[j, 0] + cnt[j, 1]
        rem = max(abs(int(c[r]) - (D1 / p + (Dp if r == c0 else 0))) for r in range(p)) * p / tot
        fac0 = int(c[c0]) * p / tot
        cells.append("p=%d: class %d factor %.3f, F p^2 = %.2f, rem %.3f, type-p off class: %d" % (p, c0, fac0, Dp / tot * p * p, rem, int(bad[j])))
    print("%s beta = %2d: %d triples; %s" % (label, beta, tot, "; ".join(cells)))


def main():
    X = 1 << 22
    size = 1 << 22
    psi = np.zeros(size, np.int64)
    for h in range(1, size):
        k = h >> 1
        psi[h] = psi[k] + k if h % 2 == 0 else psi[k] + k // 2 + 1
    for beta in (2, 3, 5, 10):
        run(psi, beta, X, 1 << 12, "box h < 2^12 (session S10):")
        run(psi, beta, X, size, "COMPLETE N_beta <= 2^22:    ")


if __name__ == "__main__":
    main()
