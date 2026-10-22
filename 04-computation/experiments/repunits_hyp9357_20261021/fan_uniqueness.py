"""Least-complexity realizations of the Mersenne numbers (session opus-2026-10-21-S51; THM-4748 (6)).
For 5 <= k <= 11, the strong classes with H = 2^k - 1 and feedback arc number <= 3 (the minimum, by THM-4748 (3)),
from the complete strong inventory below 2100, compared with the Fibonacci fan F_k and its converse.
Prints ALL CHECKS PASSED."""
import os, sys, time
from multiprocessing import Pool
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'h_spectrum_multiplicities_20261021', 'audit_S50'))
import alib as L                                       # noqa: E402
from b3_ear_census_625 import grow                     # noqa: E402

FAIL = []
T0 = time.time()


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


@nb.njit(cache=True)
def fa_one(n, adj):
    N = 1 << n
    f = np.full(N, 10 ** 9, dtype=np.int64)
    f[0] = 0
    for S in range(1, N):
        best = 10 ** 9
        x = S
        while x:
            v = 0
            while not (x >> v) & 1:
                v += 1
            x &= x - 1
            T = S & ~(1 << v)
            m = adj[v] & T
            cnt = 0
            while m:
                m &= m - 1
                cnt += 1
            if f[T] + cnt < best:
                best = f[T] + cnt
        f[S] = best
    return f[N - 1]


def fan_adj(k):
    a = k // 2; b = k + 1 - a; m = k + 2
    out = [0] * m
    for i in range(1, m):
        for j in range(i + 1, m):
            out[i] |= 1 << j
    for j in range(1, m):
        if j in (a, b, a + b):
            out[j] |= 1
        else:
            out[0] |= 1 << j
    return out


if __name__ == '__main__':
    LIM = 2100
    inv = {3: {L.canon(L.chain((2,))): 3}}
    with Pool(12) as pool:
        n = 4
        while True:
            par = [c for c, h in inv[n - 1].items() if h <= LIM]
            got, _ = grow(pool, n, par, cap=LIM) if par else ({}, None)
            got = {c: h for c, h in got.items() if h <= LIM}
            if not got:
                break
            inv[n] = got
            n += 1
    print(f'     strong inventory below {LIM}: {sum(len(v) for v in inv.values())} classes [{time.time() - T0:.0f}s]')
    res = {}
    for k in range(5, 12):
        h = 2 ** k - 1
        F = fan_adj(k)
        check(L.H(F) == h, f'k = {k}: H(F_k) = {h}')
        cf, cc = L.canon(F), L.canon(L.converse(F))
        low = []
        for m in sorted(inv):
            for code, hh in inv[m].items():
                if hh == h:
                    f = int(fa_one(m, np.array(L.decode(m, code), dtype=np.int64)))
                    if f <= 3:
                        low.append((m, f, m == k + 2 and code in (cf, cc)))
        res[k] = low
        print(f'     k = {k}, H = {h}: {len(low)} strong classes with fa <= 3 (orders {sorted(set(x[0] for x in low))}); '
              f'the fan or its converse among them: {sum(1 for x in low if x[2])}; fan self-converse: {cf == cc}')
    check(all(min(f for _, f, _ in res[k]) == 3 for k in res), 'the least feedback arc number is 3 for every 5 <= k <= 11')
    only = [k for k in res if all(x[2] for x in res[k])]
    check(only == [6, 9, 11], f'the fan and its converse are the only least-complexity realizations exactly for k in {only}')
    print()
    print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
