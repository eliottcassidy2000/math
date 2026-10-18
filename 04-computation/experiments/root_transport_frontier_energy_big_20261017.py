#!/usr/bin/env python3
"""Frontier energy E_Y(N)/N for every N < 2^(KMAX+1): memory-lean version for the HYP-9290 pre-registered scan (S42).

Same quantities as root_transport_frontier_energy_20261017.py (stop every source at its first seed <= Y; c_Y(u) frontier
coefficients; E_Y(N) = sum_u c_Y(u)^2; R(N) = sum_u f(u) c_Y(u) checked against H(4N) - H(N) at every N), but:
clock parities are stored as int8, the frontier map as int32 and computed in parallel, H(4N) - H(N) is accumulated
incrementally, and per-octave statistics (max / mean / min of E_Y/N with arg-max, max of R^2/N with arg-max) are
accumulated inside the loop. Memory is about 4 bytes per integer up to 4 * 2^(KMAX+1).
Usage: python root_transport_frontier_energy_big_20261017.py KMAX Y1 Y2 ...
Session opus-2026-10-17-S42.
"""
import sys, os, json, time
import numpy as np
from numba import njit, prange


@njit(cache=True)
def clock_parity8(M):
    par = np.zeros(M + 1, dtype=np.int8)
    for n in range(2, M + 1):
        x = n
        st = 0
        while x >= n:
            if x % 2 == 0:
                x //= 2
            else:
                x = (3 * x + 1) // 2
            st += 1
        par[n] = (st + par[x]) & 1
    return par


@njit(parallel=True, cache=True)
def frontier_map(M, Y):
    n = (M + 1) // 2
    fu = np.zeros(n, dtype=np.int32)
    fs = np.zeros(n, dtype=np.int8)
    for i in prange(n):
        s = 2 * i + 1
        if s % 8 == 5:
            continue
        sign = 1
        x = s
        while x > Y:
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
        fu[i] = x
        fs[i] = sign
    return fu, fs


@njit(cache=True)
def seed(u):
    while u % 8 == 5:
        u = (u - 1) // 4
    return u


@njit(cache=True)
def wfun(m):
    return 2 if m % 3 == 2 else -1


@njit(cache=True)
def thr(m):
    r = m % 3
    return 1 if r == 0 else (-2 if r == 1 else 1)


@njit(cache=True)
def run(Nmax, Y, fu, fs, par):
    c = np.zeros(Y + 1, dtype=np.int64)
    E = 0
    R = 0
    H4 = 0
    H1 = 0
    ok = True
    K = 64
    emax = np.zeros(K); eargmax = np.zeros(K, dtype=np.int64); esum = np.zeros(K); emin = np.full(K, 1e300)
    rmax = np.zeros(K); rargmax = np.zeros(K, dtype=np.int64); cnt = np.zeros(K, dtype=np.int64)
    for N in range(1, Nmax):
        m1 = 4 * N - 3
        m2 = 4 * N - 1
        m3 = 2 * N - 1
        for t in range(3):
            if t == 0:
                m = m1; amt = wfun(m1)
            elif t == 1:
                m = m2; amt = wfun(m2)
            else:
                m = m3; amt = thr(m3)
            s = seed(m)
            i = (s - 1) // 2
            u = fu[i]
            delta = fs[i] * amt
            old = c[u]
            E += 2 * old * delta + delta * delta
            c[u] = old + delta
            R += (1 - 2 * par[u]) * delta
        for k in range(4 * N - 3, 4 * N + 1):
            H4 += wfun(k) * (1 - 2 * par[k])
        H1 += wfun(N) * (1 - 2 * par[N])
        if R != H4 - H1:
            ok = False
        k = 0
        x = N
        while x > 1:
            x >>= 1
            k += 1
        e = E / N
        r2 = R * R / N
        if e > emax[k]:
            emax[k] = e; eargmax[k] = N
        if e < emin[k]:
            emin[k] = e
        esum[k] += e
        cnt[k] += 1
        if r2 > rmax[k]:
            rmax[k] = r2; rargmax[k] = N
    return ok, emax, eargmax, esum, emin, rmax, rargmax, cnt


def main():
    KMAX = int(sys.argv[1]); Ys = [int(a) for a in sys.argv[2:]]
    Nmax = 2 ** (KMAX + 1); M = 4 * Nmax
    t0 = time.time()
    par = clock_parity8(M)
    print('clock parities to %d [%.0f s]' % (M, time.time() - t0), flush=True)
    res = {}
    for Y in Ys:
        fu, fs = frontier_map(M, Y)
        ok, emax, eargmax, esum, emin, rmax, rargmax, cnt = run(Nmax, Y, fu, fs, par)
        FY = sum(1 for s in range(1, Y + 1, 2) if s % 8 != 5)
        rows = []
        for k in range(8, KMAX + 1):
            rows.append(dict(k=k, EY_max=float(emax[k]), EY_argmax=int(eargmax[k]), EY_mean=float(esum[k] / cnt[k]),
                             EY_min=float(emin[k]), R2N_max=float(rmax[k]), R2N_argmax=int(rargmax[k])))
        res[Y] = dict(F_Y=FY, invariant_ok=bool(ok), octaves=rows)
        print('Y = %d (F_Y = %d; R(N) invariant at every N: %s) [%.0f s]' % (Y, FY, ok, time.time() - t0), flush=True)
        for r in rows:
            print('   octave 2^%2d: E_Y/N mean %.4f  max %.4f (N = %d)  min %.4f | max R^2/N %.3f (N = %d)'
                  % (r['k'], r['EY_mean'], r['EY_max'], r['EY_argmax'], r['EY_min'], r['R2N_max'], r['R2N_argmax']),
                  flush=True)
        del fu, fs
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump({str(k): v for k, v in res.items()}, open(os.path.join(here, '..', '..', '05-knowledge', 'results',
              'root_transport_frontier_energy_big_%d_20261017.json' % KMAX), 'w'), indent=1)


if __name__ == '__main__':
    main()
