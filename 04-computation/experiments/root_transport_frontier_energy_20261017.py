#!/usr/bin/env python3
"""Penrose-type frontier for the ROOT transport: stop every source at its first seed <= Y (S42).

Idea (pattern-matched from the conformal-flow proof of the Riemannian Penrose inequality, which replaces a singular
horizon by an outer frontier and finishes with an enclosing comparison): instead of transporting every source to ROOT,
stop each source x at the first seed u <= Y on its odd-transport path (sources <= Y stop at once). The invariant
sum_u c(u) f(u) = R(N) survives stopping, so with frontier coefficients c_Y(u) = sum_{x stops at u} sigma(x -> u) c_0(x),

    R(N) = sum_{u <= Y} f(u) c_Y(u),      R(N)^2 <= F_Y * E_Y(N),      E_Y(N) = sum_u c_Y(u)^2,

by Cauchy-Schwarz, where F_Y is the number of seeds <= Y. The coherent merging near ROOT (the termination regime of
THM-4667 / HYP-9288) is replaced by freezing; moving/frozen arrivals still enter E_Y.
Scope repair: first-hit maps below are computed only for the finite sources whose paths complete;
the mathematical frontier identity uses f=0 on nonhitting ROOT-complement sources, not universal entry. The script computes E_Y(N)/N for every
N < 2^(KMAX+1) incrementally (three source updates per N), for several fixed Y, and checks R(N) against H(4N) - H(N).
Usage: python root_transport_frontier_energy_20261017.py KMAX Y1 Y2 ...
Session opus-2026-10-17-S42.
"""
import sys, os, json, time
import numpy as np
from numba import njit


@njit(cache=True)
def seed(u):
    while u % 8 == 5:
        u = (u - 1) // 4
    return u


@njit(cache=True)
def step(s):
    """Odd transport of a seed s > 1: returns (P(s), sign)."""
    u = 3 * s + 1
    a = 0
    while u % 2 == 0:
        u //= 2
        a += 1
    return seed(u), (1 if a == 2 else -1)


@njit(cache=True)
def frontier_map(M, Y):
    """For every odd seed s <= M (index (s-1)//2), the first seed <= Y on its path and the accumulated sign."""
    n = (M + 1) // 2
    fu = np.zeros(n, dtype=np.int64)
    fs = np.zeros(n, dtype=np.int8)
    for i in range(n):
        s = 2 * i + 1
        if s % 8 == 5:
            continue
        sign = 1
        x = s
        while x > Y:
            x, sg = step(x)
            sign *= sg
        fu[i] = x
        fs[i] = sign
    return fu, fs


@njit(cache=True)
def clock_parity(M):
    clk = np.zeros(M + 1, dtype=np.int64)
    for n in range(2, M + 1):
        x = n; st = 0
        while x >= n:
            x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
            st += 1
        clk[n] = st + clk[x]
    return clk


@njit(cache=True)
def wfun(m):
    return 2 if m % 3 == 2 else -1


@njit(cache=True)
def thr(m):
    r = m % 3
    return 1 if r == 0 else (-2 if r == 1 else 1)


@njit(cache=True)
def run(Nmax, Y, fu, fs, fY):
    """fY[u] = f(u) for u <= Y. Returns E_Y(N)/N and R_Y(N) (should equal R(N)) for all N."""
    c = np.zeros(Y + 1, dtype=np.int64)
    E = 0
    out = np.zeros(Nmax, dtype=np.float64)
    Rv = np.zeros(Nmax, dtype=np.int64)
    R = 0
    for N in range(1, Nmax):
        m1 = 4 * N - 3; m2 = 4 * N - 1; m3 = 2 * N - 1
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
            R += fY[u] * delta
        out[N] = E / N
        Rv[N] = R
    return out, Rv


def main():
    KMAX = int(sys.argv[1]); Ys = [int(a) for a in sys.argv[2:]]
    Nmax = 2 ** (KMAX + 1); M = 4 * Nmax
    t0 = time.time()
    clk = clock_parity(M)
    f = 1 - 2 * (clk & 1)
    n = np.arange(M + 1)
    wv = np.where(n % 3 == 2, 2, -1); wv[0] = 0
    H = np.cumsum(wv * f)
    Ns = np.arange(1, Nmax)
    Rdir = H[4 * Ns] - H[Ns]
    print('clock parities [%.0f s]' % (time.time() - t0), flush=True)
    res = {}
    for Y in Ys:
        fu, fs = frontier_map(M, Y)
        fY = np.zeros(Y + 1, dtype=np.int64)
        fY[1:Y + 1] = f[1:Y + 1]
        EY, RY = run(Nmax, Y, fu, fs, fY)
        ok = np.array_equal(RY[1:], Rdir)
        FY = sum(1 for s in range(1, Y + 1, 2) if s % 8 != 5)
        rows = []
        for k in range(8, KMAX + 1):
            sl = slice(2 ** k, 2 ** (k + 1))
            e = EY[sl]
            r2 = Rdir[2 ** k - 1:2 ** (k + 1) - 1].astype(np.float64) ** 2 / Ns[2 ** k - 1:2 ** (k + 1) - 1]
            ratio = r2 / np.maximum(FY * e, 1e-300)
            rows.append(dict(k=k, EY_mean=float(e.mean()), EY_max=float(e.max()), EY_min=float(e.min()),
                             R2N_max=float(r2.max()), CS_ratio_max=float(ratio.max())))
        res[Y] = dict(F_Y=FY, invariant_ok=bool(ok), octaves=rows)
        print('Y = %d (F_Y = %d seeds; invariant R(N) check: %s) [%.0f s]' % (Y, FY, ok, time.time() - t0), flush=True)
        for r in rows:
            print('   octave 2^%2d: E_Y/N mean %.4f  max %.4f  min %.4f | max R^2/N %.3f | max R^2/(F_Y E_Y) %.4f'
                  % (r['k'], r['EY_mean'], r['EY_max'], r['EY_min'], r['R2N_max'], r['CS_ratio_max']), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump({str(k): v for k, v in res.items()}, open(os.path.join(here, '..', '..', '05-knowledge', 'results',
                                                                        'root_transport_frontier_energy_%d_20261017.json' % KMAX), 'w'), indent=1)


if __name__ == '__main__':
    main()
