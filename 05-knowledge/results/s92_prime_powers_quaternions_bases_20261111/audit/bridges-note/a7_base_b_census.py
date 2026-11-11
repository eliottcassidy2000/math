#!/usr/bin/env python3
"""Auditor bridges-note, script A7 (independent code): spot replication of the base-b census quoted in section 5 of
the note and in THM-5134 (H) (every_base.out (4)): sheet s >= 1, N_b = X1 X2 + b 2^(s-1) y1 y2 <= 2^24 over ordered
triples (h1, h2, s), X = 1 + b psi(h), y = 2h + 1.

Reported: number of triples, N_b = 1 (mod b), class factor p * share at the class 0, range of the other classes,
share F of the triples whose smaller state is p-divisible (p | X and p | y), the exactness of "type p => class 0",
and the remainder of the first-order law  #{N = c} = (D_1 + [c = 0] p D_p) / p.
"""
import sys
import numpy as np
from numba import njit


@njit(cache=True)
def census(psi, b, X, primes, cnt, bad, digit):
    H = psi.shape[0]
    tot = 0
    notone = 0
    s = 1
    while b * (1 << (s - 1)) + 1 <= X:
        w0 = b * (1 << (s - 1))
        for h1 in range(H):
            y1 = 2 * h1 + 1
            X1 = 1 + b * psi[h1]
            if w0 * y1 * y1 + X1 > X:
                # h2 >= h1 so the seam term alone already exceeds X for larger h1 as well
                if w0 * y1 * y1 > X:
                    break
            for h2 in range(h1, H):
                y2 = 2 * h2 + 1
                if w0 * y1 * y2 + X1 > X:
                    break
                N = X1 * (1 + b * psi[h2]) + w0 * y1 * y2
                if N > X:
                    continue
                w = 1 if h2 == h1 else 2
                tot += w
                if N % b != 1:
                    notone += 1
                digit[((N - 1) // b) % b] += w
                for j in range(primes.shape[0]):
                    p = primes[j]
                    t = 1 if (X1 % p == 0 and y1 % p == 0) else 0
                    cnt[j, t, N % p] += w
                    if t == 1 and N % p != 0:
                        bad[j] += 1
        s += 1
    return tot, notone


def main():
    X = 1 << 24
    size = 1 << 23
    psi = np.zeros(size, np.int64)
    for h in range(1, size):
        k = h >> 1
        psi[h] = psi[k] + k if h % 2 == 0 else psi[k] + k // 2 + 1
    for b in (2, 3, 4, 5, 6, 10):
        primes = np.array([p for p in (3, 5, 7, 11, 13) if b % p], np.int64)
        # the session's table stops at h < 2^21 (every_base.py: q32 = q[: (1 << 21)]): reproduce its triple count
        cnt = np.zeros((len(primes), 2, 16), np.int64)
        bad = np.zeros(len(primes), np.int64)
        digit = np.zeros(b, np.int64)
        tot21, _ = census(psi[: (1 << 21)].copy(), b, X, primes, cnt, bad, digit)
        print("b = %2d: triples with h1, h2 < 2^21 (the session's table): %d" % (b, tot21))
        cnt = np.zeros((len(primes), 2, 16), np.int64)
        bad = np.zeros(len(primes), np.int64)
        digit = np.zeros(b, np.int64)
        tot, notone = census(psi, b, X, primes, cnt, bad, digit)
        flat = max(abs(int(digit[d]) * b / tot - 1) for d in range(b))
        line = ["b = %2d: %d triples; N_b != 1 (mod b): %d; second digit flat to %.3f" % (b, tot, notone, flat)]
        for j, p in enumerate(primes):
            p = int(p)
            c = cnt[j, 0] + cnt[j, 1]
            fac = [int(c[r]) * p / tot for r in range(p)]
            D1, Dp = int(cnt[j, 0].sum()), int(cnt[j, 1].sum())
            rem = max(abs(int(c[r]) - (D1 / p + (Dp if r == 0 else 0))) for r in range(p)) * p / tot
            line.append("p=%d: class 0 %.3f, others [%.3f, %.3f], F p^2 = %.2f, type-p triples off class 0: %d, rem %.3f"
                        % (p, fac[0], min(fac[1:]), max(fac[1:]), Dp / tot * p * p, int(bad[j]), rem))
        print("; ".join(line))


if __name__ == "__main__":
    main()
