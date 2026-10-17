#!/usr/bin/env python3
"""Rotation (Beatty) orders of the Paley tournament: rank x in Z/p by the fractional part {alpha x}; alpha = a/p gives
the affine orders (excess +-p h(-p), THM-4666). Which rational alpha = k/(m p) with small m give the largest excess,
and do those excesses follow class numbers of other discriminants? (S40, after the independent audit found rotation
orders beating the natural order at p = 167, 263.)

For each prime p = 3 (mod 4) in the list, scan alpha = k/(m p) for 1 <= m <= MMAX, 1 <= k < m p, gcd(k, m) = 1 (the
order only depends on {alpha x} for integer x, so k is taken mod m p), and report the best excess per m, the best overall,
and the class numbers h(-p), h(-3p), h(-4p) (= h(-p) or 3h(-p)), h(-8p) for comparison.
Usage: python paley_rotation_orders_20261017.py MMAX p1 p2 ...
Session opus-2026-10-17-S40.
"""
import sys, os, json, math
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paley_class_number_ranking_20261017 import kronecker_table, class_number


@njit(cache=True)
def excess_of(order, chi, p):
    s = 0
    n = order.shape[0]
    for i in range(n):
        u = order[i]
        for j in range(i + 1, n):
            s += chi[(order[j] - u) % p]
    return s


@njit(cache=True)
def scan(chi, p, m):
    """Best excess over k = 1..m p - 1 for alpha = k/(m p); ties in {alpha x} broken by x (stable)."""
    best = -10 ** 9; bestk = -1
    q = m * p
    keys = np.empty(p, dtype=np.int64)
    for k in range(1, q):
        for x in range(p):
            keys[x] = (k * x) % q              # {k x / (m p)} * m p
        order = np.argsort(keys, kind='mergesort')
        e = excess_of(order, chi, p)
        if e > best:
            best = e; bestk = k
    return best, bestk


def main():
    MMAX = int(sys.argv[1])
    primes = [int(a) for a in sys.argv[2:]]
    out = []
    for p in primes:
        chi = kronecker_table(-p)
        h = class_number(-p)
        hs = {d: class_number(-d * p) if (-d * p) % 4 in (0, 1) else None for d in (3, 4, 8, 11)}
        rows = []
        for m in range(1, MMAX + 1):
            b, k = scan(chi, p, m)
            rows.append((m, int(b), int(k)))
        bm = max(rows, key=lambda r: r[1])
        print('p = %d: natural p h = %d (h = %d); class numbers h(-3p) = %s, h(-4p) = %s, h(-8p) = %s, h(-11p) = %s'
              % (p, p * h, h, hs[3], hs[4], hs[8], hs[11]))
        print('   best excess by m: %s' % [(m, b) for m, b, k in rows])
        print('   best overall: m = %d, k = %d, excess %d (= %.3f p h, %.4f p^1.5)' % (bm[0], bm[2], bm[1], bm[1] / (p * h), bm[1] / p ** 1.5),
              flush=True)
        out.append(dict(p=p, h=h, hs={str(k): v for k, v in hs.items()}, rows=rows))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'paley_rotation_orders_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
