#!/usr/bin/env python3
"""clock_duality_oct15's odd-source fold R = I + B (marked_contacts_root_current_20261017e.md, R1) on the FULL
ROOT current, every octave (S40).

With chi(m) = 0, 1, -1 for m = 0, 1, 2 (mod 3), w(m) = 3[m = 2 (mod 3)] - 1 and f(m) = (-1)^d(m) (every m in range
reaches ROOT 1, so the full profile is known), R1 says

    a_N = H(4N) - H(N) = I(N) + B(N),   I(N) = -3 sum_(m <= 2N, m odd) chi(m) f(m),
                                        B(N) = sum_(2N < m <= 4N, m odd) w(m) f(m).

This script checks the identity at every row of every octave 2^KMIN..2^KMAX and reports the normalized energies
||I||^2/X^2, ||B||^2/X^2, <I,B>/X^2, ||R||^2/X^2 and the correlation <I,B>/(||I|| ||B||), extending the four staged
rows of the source note (k = 10, 14, 18) to the final state at every k. All energies are exact integer numerators over
X^2 (int64 accumulation, asserted below 2^62); the printed decimals are their roundings.
Usage: python collatz_root_current_fold_full_20261017.py KMIN KMAX
Session opus-2026-10-17-S40.
"""
import sys, os, json
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collatz_root_current_exact_20261017 import parity_clock


@njit(cache=True)
def odd_prefixes(d):
    """Pchi[x] = sum_(m <= x, m odd) chi(m) f(m); Pw[x] = sum_(m <= x, m odd) w(m) f(m); H[x] = sum_(n <= x) w(n) f(n)."""
    L = d.shape[0]
    Pchi = np.zeros(L, dtype=np.int32)
    Pw = np.zeros(L, dtype=np.int32)
    H = np.zeros(L, dtype=np.int32)
    sc = 0; sw = 0; sh = 0
    for n in range(1, L):
        f = 1 - 2 * d[n]
        r = n % 3
        w = 2 if r == 2 else -1
        sh += w * f
        if n % 2 == 1:
            sw += w * f
            if r == 1:
                sc += f
            elif r == 2:
                sc -= f
        Pchi[n] = sc; Pw[n] = sw; H[n] = sh
    return Pchi, Pw, H


@njit(cache=True)
def octave(Pchi, Pw, H, X):
    I = np.empty(X, dtype=np.int64); B = np.empty(X, dtype=np.int64); bad = 0
    for i in range(X):
        N = X + i
        I[i] = -3 * np.int64(Pchi[2 * N])
        B[i] = np.int64(Pw[4 * N]) - np.int64(Pw[2 * N])
        if I[i] + B[i] != np.int64(H[4 * N]) - np.int64(H[N]):
            bad += 1
    return I, B, bad


def main():
    KMIN, KMAX = int(sys.argv[1]), int(sys.argv[2])
    L = 4 * 2 ** (KMAX + 1) + 1
    d = parity_clock(L)
    Pchi, Pw, H = odd_prefixes(d)
    del d
    rows = []
    print(' k   ||I||^2/X^2  ||B||^2/X^2  <I,B>/X^2  ||R||^2/X^2   corr(I,B)  R1 failures')
    for k in range(KMIN, KMAX + 1):
        X = 2 ** k
        I, B, bad = octave(Pchi, Pw, H, X)
        mx = max(int(np.max(np.abs(I))), int(np.max(np.abs(B))))
        assert X * mx * mx < 2 ** 62                       # no int64 wraparound in the sums below
        ii_n = int(np.sum(I * I)); bb_n = int(np.sum(B * B)); ib_n = int(np.sum(I * B))     # exact int64 numerators
        rr_n = ii_n + bb_n + 2 * ib_n
        ii, bb, ib, rr = ii_n / X ** 2, bb_n / X ** 2, ib_n / X ** 2, rr_n / X ** 2
        corr = ib / np.sqrt(ii * bb) if ii > 0 and bb > 0 else float('nan')
        rows.append(dict(k=k, II=ii, BB=bb, IB=ib, RR=rr, corr=corr, failures=int(bad), II_num=ii_n, BB_num=bb_n,
                         IB_num=ib_n, RR_num=rr_n, denominator=X ** 2))
        print('%2d  %11.6f  %11.6f  %10.6f  %11.6f  %9.4f  %d' % (k, ii, bb, ib, rr, corr, bad), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(rows, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'collatz_root_current_fold_full_20261017.json'), 'w'),
              indent=1)


if __name__ == '__main__':
    main()
