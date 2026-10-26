#!/usr/bin/env python3
"""Drift of class-wise untouchable shares toward the Pollack-Pomerance values (opus-2026-10-26-S57).
For even k in windows (2^(e-1), 2^e], e = 14..30, and classes a = gcd(k, 2520): empirical zero share vs PP's exp(-a/s(a)).
Usage: python pp_class_drift_20261026.py bigdir
"""
import sys, os, json, math
import numpy as np
import numba as nb
HERE = os.path.dirname(os.path.abspath(__file__))
A = 2520
CLASSES = [2, 4, 6, 8, 10, 12, 14, 24, 30, 40, 120, 840]

def sigma(n):
    return sum(d for d in range(1, n + 1) if n % d == 0)

@nb.njit(cache=False)
def window(mult, j0, j1, A, idx_of, cnt, zr):
    for j in range(j0, j1):
        k = 2 * j
        a = k
        b = A
        while b:
            a, b = b, a % b
        c = idx_of[a]
        if c >= 0:
            cnt[c] += 1
            if mult[j] == 0:
                zr[c] += 1

if __name__ == '__main__':
    mult = np.load(os.path.join(sys.argv[1], 'mult_even_1000000000.npy'))
    idx_of = np.full(A + 1, -1, np.int64)
    for i, d in enumerate(CLASSES):
        idx_of[d] = i
    pp = [math.exp(-d / (sigma(d) - d)) for d in CLASSES]
    print('class a      : ' + ' '.join(f'{d:>7}' for d in CLASSES))
    print('PP exp(-a/s) : ' + ' '.join(f'{v:7.4f}' for v in pp))
    rows = {}
    for e in range(14, 31):
        hi = min(2 ** e, 10 ** 9)
        lo = 2 ** (e - 1)
        if lo >= 10 ** 9:
            break
        cnt = np.zeros(len(CLASSES), np.int64)
        zr = np.zeros(len(CLASSES), np.int64)
        window(mult, lo // 2 + 1, hi // 2 + 1, A, idx_of, cnt, zr)
        sh = [zr[i] / cnt[i] for i in range(len(CLASSES))]
        rows[e] = sh
        print(f'window 2^{e:<2}  : ' + ' '.join(f'{v:7.4f}' for v in sh), flush=True)
    json.dump({'classes': CLASSES, 'pp': pp, 'windows': {str(k): v for k, v in rows.items()}},
              open(os.path.join(HERE, 'pp_class_drift.json'), 'w'), indent=1)
