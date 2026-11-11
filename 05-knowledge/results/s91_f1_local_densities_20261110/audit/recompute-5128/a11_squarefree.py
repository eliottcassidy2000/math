# Auditor recompute-5128, extra: the squarefree form of (C) (stated in THM-5128 (C) without a numerical check):
# for a squarefree odd modulus m the class c has share
#     sum over d | m of F^(d)(X) * prod_{p | d} [c = c_sigma (mod p)] * prod_{p | m/d} (1/p),
# F^(d) = share of triples whose SMALLER state has gcd(g, m) = d   (triples with h1 = h2 are put in d = 1 here).
# Exact census of all triples with N < 2^B; prints m * max_c |share(c) - predicted(c)| per sheet.
import sys
import numpy as np
from numba import njit
from a0_defs import psi_table32


@njit(cache=True)
def census_m(psi, B, m, primes):
    npr = primes.shape[0]
    hsize = psi.shape[0]
    code = np.zeros(hsize, dtype=np.uint8)
    for h in range(hsize):
        x = 2 * np.int64(psi[h]) + 1
        y = 2 * np.int64(h) + 1
        c = 0
        for j in range(npr):
            if y % primes[j] == 0 and x % primes[j] == 0:
                c |= 1 << j
        code[h] = c
    cnt = np.zeros((2, 1 << npr, m), dtype=np.int64)
    for s in range(0, B):
        sheet = 0 if s == 0 else 1
        corr = 1 if s == 0 else 0
        Xs = (np.int64(1) << B) - 1 + corr
        h1 = np.int64(0)
        while True:
            y1 = 2 * h1 + 1
            y2crude = (4 * Xs) // ((1 + (np.int64(1) << (s + 2))) * y1)
            if (y2crude - 1) // 2 < h1:
                break
            x1 = 2 * np.int64(psi[h1]) + 1
            b = y1 << s
            y2max = (2 * Xs) // (x1 + 2 * b)
            h2max = (y2max - 1) // 2
            c1 = code[h1]
            for h2 in range(h1, h2max + 1):
                M = x1 * (2 * np.int64(psi[h2]) + 1) + b * (2 * h2 + 1)
                if M <= Xs:
                    N = M - corr
                    if h2 > h1:
                        cnt[sheet, c1, N % m] += 2
                    else:
                        cnt[sheet, 0, N % m] += 1
            h1 += 1
    return cnt


def main():
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    psi = psi_table32(((1 << B) + 1) // 3 + 4)
    for primes in ([3, 5], [5, 11], [5, 17], [3, 5, 11]):
        m = int(np.prod(primes))
        cnt = census_m(psi, B, m, np.array(primes, dtype=np.int64))
        print(f"== m = {m} = {' * '.join(map(str, primes))}, N < 2^{B} ==")
        for sheet in (0, 1):
            T = cnt[sheet].sum()
            shares = cnt[sheet].sum(axis=0) / T
            F = cnt[sheet].sum(axis=1) / T  # indexed by the bitmask of primes dividing g(smaller state)
            pred = np.zeros(m)
            for c in range(m):
                tot = 0.0
                for mask in range(1 << len(primes)):
                    term = F[mask]
                    for j, p in enumerate(primes):
                        cs = (p - 1) if sheet == 0 else 0
                        if (mask >> j) & 1:
                            term *= 1.0 if c % p == cs else 0.0
                        else:
                            term *= 1.0 / p
                    tot += term
                pred[c] = tot
            # the fully distinguished class
            cdist = [c for c in range(m) if all(c % p == ((p - 1) if sheet == 0 else 0) for p in primes)][0]
            flat = np.full(m, 1.0 / m)
            print(f"   sheet {sheet}: T = {int(T)}; m*share at the fully distinguished class {cdist}: measured {m * shares[cdist]:.5f}, product law {m * pred[cdist]:.5f} (flat = 1); "
                  f"m*max_c|share - product law| = {m * np.abs(shares - pred).max():.5f}; m*max_c|share - 1/m| = {m * np.abs(shares - flat).max():.5f}; sum of predictions = {pred.sum():.6f}")


if __name__ == "__main__":
    main()
