#!/usr/bin/env python3
"""Frontier energy E_Y(N)/N for every N < 2^(KMAX+1), with int64 overflow guards, for KMAX up to 29 (S42).

Same algorithm and output format as root_transport_frontier_energy_big_20261017.py (which it reuses for the main loop),
but the clock parities and the frontier map are computed by guarded copies: before every 3x + 1 the value x is compared
with (2^63 - 2)/3 and a flag is raised if it is exceeded, so a silent int64 wrap-around is impossible. (For starting
values below 2^32 the largest Collatz excursion is about 7.7e18, under 2^63, so no flag is expected at KMAX = 29.)
Usage: python root_transport_frontier_energy_huge_20261018.py KMAX Y1 Y2 ...
Session opus-2026-10-17-S42.
"""
import sys, os, json, time
import numpy as np
from numba import njit, prange
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from root_transport_frontier_energy_big_20261017 import run

LIM = (2 ** 63 - 2) // 3


@njit(cache=True)
def clock_parity8_guarded(M):
    par = np.zeros(M + 1, dtype=np.int8)
    bad = 0
    for n in range(2, M + 1):
        x = n
        st = 0
        over = False
        while x >= n:
            if x % 2 == 0:
                x //= 2
            else:
                if x > LIM:
                    over = True
                    break
                x = (3 * x + 1) // 2
            st += 1
        if over:
            bad += 1
            par[n] = 0
        else:
            par[n] = (st + par[x]) & 1
    return par, bad


@njit(parallel=True, cache=True)
def frontier_map_guarded(M, Y):
    n = (M + 1) // 2
    fu = np.zeros(n, dtype=np.int32)
    fs = np.zeros(n, dtype=np.int8)
    badv = np.zeros(n, dtype=np.int8)
    for i in prange(n):
        s = 2 * i + 1
        if s % 8 == 5:
            continue
        sign = 1
        x = s
        while x > Y:
            if x > LIM:
                badv[i] = 1
                break
            u = 3 * x + 1
            a = 0
            while u % 2 == 0:
                u //= 2
                a += 1
            while u % 8 == 5:
                u = (u - 1) // 4
            x = u
            if a != 2:
                sign = -sign
        fu[i] = x if x <= Y else 0
        fs[i] = sign
    return fu, fs, badv.sum()


def main():
    KMAX = int(sys.argv[1]); Ys = [int(a) for a in sys.argv[2:]]
    Nmax = 2 ** (KMAX + 1); M = 4 * Nmax
    t0 = time.time()
    par, bad = clock_parity8_guarded(M)
    print('clock parities to %d [%.0f s]; overflow flags: %d' % (M, time.time() - t0, bad), flush=True)
    res = {}
    for Y in Ys:
        fu, fs, badf = frontier_map_guarded(M, Y)
        ok, emax, eargmax, esum, emin, rmax, rargmax, cnt = run(Nmax, Y, fu, fs, par)
        FY = sum(1 for s in range(1, Y + 1, 2) if s % 8 != 5)
        rows = []
        for k in range(8, KMAX + 1):
            rows.append(dict(k=k, EY_max=float(emax[k]), EY_argmax=int(eargmax[k]), EY_mean=float(esum[k] / cnt[k]),
                             EY_min=float(emin[k]), R2N_max=float(rmax[k]), R2N_argmax=int(rargmax[k])))
        res[Y] = dict(F_Y=FY, invariant_ok=bool(ok), overflow_flags=int(bad + badf), octaves=rows)
        print('Y = %d (F_Y = %d; R(N) invariant at every N: %s; overflow flags %d) [%.0f s]'
              % (Y, FY, ok, bad + badf, time.time() - t0), flush=True)
        for r in rows:
            print('   octave 2^%2d: E_Y/N mean %.4f  max %.4f (N = %d)  min %.4f | max R^2/N %.3f (N = %d)'
                  % (r['k'], r['EY_mean'], r['EY_max'], r['EY_argmax'], r['EY_min'], r['R2N_max'], r['R2N_argmax']),
                  flush=True)
        del fu, fs
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump({str(k): v for k, v in res.items()}, open(os.path.join(here, '..', '..', '05-knowledge', 'results',
              'root_transport_frontier_energy_huge_%d_20261018.json' % KMAX), 'w'), indent=1)


if __name__ == '__main__':
    main()
