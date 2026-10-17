#!/usr/bin/env python3
"""The exact ROOT current a_N = H(4N) - H(N) of clock_duality_oct15's tail-debt programme, at full resolution (S40).

Definitions (escape_renewal_idoneal_20261017d.md, golden_tail_idoneal_20261017c.md): T the shortcut map
(n/2, (3n+1)/2), d(n) the number of T-steps from n to 1 (ROOT), w(n) = 3[n = 2 mod 3] - 1, H(M) = sum_(n <= M) w(n)
(-1)^d(n), and on the octave X <= N < 2X the current a_N = H(4N) - H(N). Every n in range reaches 1, so the full
current is known: U = 0 and a_N = R(N) in the notation of D1 there.

Outputs, for every octave X = 2^k, k = KMIN..KMAX (exact, all rows): the energy sum a_N^2 / X^2, the tail constant
K_X = max_t t^2 #{N : |a_N| >= t} / X^2 (target (T) of graded_tail_carrier_20261017b.md), the fraction of positive rows,
the number of sign changes, min |a_N| and max |a_N|; and a comparison with the census H at the grid points.
Usage: python collatz_root_current_exact_20261017.py KMIN KMAX
Session opus-2026-10-17-S40.
"""
import sys, json
import numpy as np
from numba import njit


@njit(cache=True)
def parity_clock(L):
    """d(n) mod 2 for 1 <= n < L (d(1) = 0), shortcut map; only the parity is stored."""
    d = np.zeros(L, dtype=np.int8)
    for n in range(2, L):
        m = n
        steps = 0
        while m >= n:
            if m % 2 == 0:
                m //= 2
            else:
                m = (3 * m + 1) // 2
            steps += 1
        d[n] = (steps + d[m]) & 1
    return d


@njit(cache=True)
def prefix_H(d):
    L = d.shape[0]
    H = np.zeros(L, dtype=np.int32)
    s = 0
    big = 0
    for n in range(1, L):
        w = 2 if n % 3 == 2 else -1
        s += w * (1 - 2 * d[n])
        if abs(s) > big:
            big = abs(s)
        H[n] = s
    if big >= 2 ** 31 - 1:
        raise ValueError('H overflow')
    return H


@njit(cache=True)
def octave_stats(H, X):
    a = np.empty(X, dtype=np.int64)
    for i in range(X):
        N = X + i
        a[i] = np.int64(H[4 * N]) - np.int64(H[N])
    return a


def main():
    KMIN, KMAX = int(sys.argv[1]), int(sys.argv[2])
    L = 4 * 2 ** (KMAX + 1) + 1
    print('table size', L, flush=True)
    d = parity_clock(L)
    H = prefix_H(d)
    rows = []
    print(' k        X     energy/X^2  K_X(exact)  frac>0  signchg   min|a|   max|a|   max|a|/X^(2/3)')
    for k in range(KMIN, KMAX + 1):
        X = 2 ** k
        a = octave_stats(H, X).astype(np.float64)
        e = float(np.sum(a * a) / X ** 2)
        s = np.sort(np.abs(a))[::-1]
        K = float(np.max(s ** 2 * np.arange(1, X + 1)) / X ** 2)
        pos = float(np.mean(a > 0)); chg = int(np.sum(np.diff(np.sign(a)) != 0))
        rows.append(dict(k=k, energy=e, K=K, frac_pos=pos, sign_changes=chg, min_abs=float(np.min(np.abs(a))),
                         max_abs=float(np.max(np.abs(a)))))
        print('%2d %9d   %9.4f   %9.4f   %6.3f  %6d  %8d  %8d   %7.3f' % (k, X, e, K, pos, chg, np.min(np.abs(a)), np.max(np.abs(a)),
                                                                     np.max(np.abs(a)) / X ** (2 / 3)), flush=True)
    # census cross-check at grid points
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    Z = np.load(os.path.join(here, '..', '..', '05-knowledge', 'results', 'collatz_census_dense_1024_2p40_20261016.npz'))
    g, Hc = Z['grid'].astype(np.int64), Z['H'].astype(np.int64)
    mask = g < L
    agree = bool(np.all(Hc[mask] == H[g[mask]]))
    print('census H at %d grid points below %d equals this H: %s' % (int(mask.sum()), L, agree))
    json.dump(dict(rows=rows, census_agree=agree, max_H=int(np.max(np.abs(H)))), open(os.path.join(here, '..', '..', '05-knowledge', 'results',
              'collatz_root_current_exact_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
