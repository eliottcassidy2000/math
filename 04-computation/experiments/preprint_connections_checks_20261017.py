#!/usr/bin/env python3
"""Checks for preprint_connections_solfeggio_digits_20261017.md (S42).

1. L2: for odd n > 1, (-1)^(sigma_T(n)) equals the product of the Legendre symbols (u/3) over the odd elements u != n of
   the T-orbit of n (T = shortcut map), for every odd n < 10^6 (exact).
2. Fibonacci digital roots mod 9 (Pisano period 24), complementary halves, and positions of 3, 6, 9.
3. L1 (Lipschitz amplification): per octave X = 2^k, compares max |R(N)| with 4 D(X)^(1/3), D(X) = sum R(N)^2, and
   records the observed max |R(N+1) - R(N)| (proved bound 8).
4. Character structure of the frontier coefficients: at given N and Y, correlation of c_Y(u) with f(u), with (u/3),
   and the energy share of c_Y by residue class of u mod 9 (non-units, Solfeggio fibres, golden classes) against the
   share of frontier seeds in each class.
Usage: python preprint_connections_checks_20261017.py KMAX "N1,N2,..." "Y1,Y2,..."
Session opus-2026-10-17-S42.
"""
import sys, os, json
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from root_transport_frontier_energy_big_20261017 import clock_parity8, frontier_map, seed, wfun, thr


@njit(cache=True)
def legendre_product_check(M):
    bad = 0
    for n in range(3, M, 2):
        # sigma_T(n) and product of (u/3) over odd orbit elements after n
        x = n
        st = 0
        prod = 1
        while x != 1:
            if x % 2 == 0:
                x //= 2
            else:
                x = (3 * x + 1) // 2
            st += 1
            if x % 2 == 1:
                r = x % 3
                if r == 2:
                    prod = -prod
                elif r == 0:
                    prod = 0
        if (1 - 2 * (st & 1)) != prod:
            bad += 1
    return bad


@njit(cache=True)
def octave_stats(Nmax, par):
    K = 64
    D = np.zeros(K)
    mx = np.zeros(K)
    H4 = 0
    H1 = 0
    Rprev = 0
    lip = 0
    for N in range(1, Nmax):
        for k in range(4 * N - 3, 4 * N + 1):
            H4 += wfun(k) * (1 - 2 * par[k])
        H1 += wfun(N) * (1 - 2 * par[N])
        R = H4 - H1
        if N > 1 and abs(R - Rprev) > lip:
            lip = abs(R - Rprev)
        Rprev = R
        k = 0
        x = N
        while x > 1:
            x >>= 1
            k += 1
        D[k] += R * R
        if abs(R) > mx[k]:
            mx[k] = abs(R)
    return D, mx, lip


@njit(cache=True)
def frontier_vector(N, Y, fu, fs):
    c = np.zeros(Y + 1, dtype=np.int64)
    for Np in range(1, N + 1):
        for t in range(3):
            if t == 0:
                m = 4 * Np - 3; amt = wfun(m)
            elif t == 1:
                m = 4 * Np - 1; amt = wfun(m)
            else:
                m = 2 * Np - 1; amt = thr(m)
            s = seed(m)
            i = (s - 1) // 2
            c[fu[i]] += fs[i] * amt
    return c


def main():
    KMAX = int(sys.argv[1]); Ns = [int(a) for a in sys.argv[2].split(',')]; Ys = [int(a) for a in sys.argv[3].split(',')]
    out = {}
    bad = legendre_product_check(10 ** 6)
    print('L2: odd n < 10^6 violating (-1)^sigma_T(n) = prod (u/3): %d' % bad, flush=True)
    out['L2_violations_below_1e6'] = int(bad)
    F = [0, 1]
    for i in range(60):
        F.append(F[-1] + F[-2])
    dr = [((F[n] - 1) % 9) + 1 if F[n] else 9 for n in range(1, 25)]
    print('Fibonacci digital roots F_1..F_24:', dr)
    comp = all((F[n] + F[n + 12]) % 9 == 0 for n in range(1, 40))
    pos369 = [n for n in range(1, 25) if F[n] % 3 == 0]
    print('F_n + F_(n+12) = 0 mod 9 for n < 40: %s; positions of 3,6,9 in the cycle: %s' % (comp, pos369), flush=True)
    out['fib_digital_roots'] = dr; out['fib_complement'] = comp; out['fib_pos369'] = pos369
    Nmax = 2 ** (KMAX + 1); M = 4 * max(Nmax, max(Ns) + 1)
    par = clock_parity8(M)
    D, mx, lip = octave_stats(Nmax, par)
    print('L1: observed max |R(N+1) - R(N)| for N < 2^%d: %d (proved bound 8)' % (KMAX + 1, lip))
    rows = []
    for k in range(6, KMAX + 1):
        X = 2 ** k
        b = 4 * D[k] ** (1 / 3)
        rows.append(dict(k=k, maxR=float(mx[k]), bound=float(b), DX_over_X2=float(D[k] / X ** 2), ratio=float(b / mx[k])))
        print('   octave 2^%2d: max|R| = %8d   4 D(X)^(1/3) = %10.1f   D(X)/X^2 = %.3f   bound/max = %.1f   X^(1/6) = %.1f'
              % (k, mx[k], b, D[k] / X ** 2, b / mx[k], X ** (1 / 6)), flush=True)
    out['L1'] = dict(lip=int(lip), octaves=rows)
    f = (1 - 2 * par.astype(np.int64))
    res = []
    for Y in Ys:
        fu, fs = frontier_map(M, Y)
        us = np.array([u for u in range(1, Y + 1) if u % 2 == 1 and u % 8 != 5])
        chi = np.array([0 if u % 3 == 0 else (1 if u % 3 == 1 else -1) for u in us])
        cls = {'nonunit {0,3,6}': [u for u in us if u % 3 == 0],
               'fibre {1,4,7}': [u for u in us if u % 3 == 1], 'fibre {2,5,8}': [u for u in us if u % 3 == 2],
               'golden +-1': [u for u in us if u % 9 in (1, 8)], 'golden +-2': [u for u in us if u % 9 in (2, 7)],
               'golden +-4': [u for u in us if u % 9 in (4, 5)]}
        for N in Ns:
            c = frontier_vector(N, Y, fu, fs)
            cu = c[us].astype(np.float64)
            E = float(cu @ cu)
            R = int((c[us] * f[us]).sum())
            kk = np.arange(N + 1, 4 * N + 1)
            Rdir = int((np.where(kk % 3 == 2, 2, -1) * f[kk]).sum())
            FY = len(us)
            rho_f = R / np.sqrt(FY * E)
            rho_chi = float((cu * chi).sum()) / np.sqrt(FY * E)
            shares = {}
            for name, lst in cls.items():
                lst = np.array(lst)
                shares[name] = (float((c[lst].astype(np.float64) ** 2).sum()) / E, len(lst) / FY)
            res.append(dict(N=N, Y=Y, E_over_N=E / N, R=R, R_direct=Rdir, rho_f=float(rho_f), rho_chi=rho_chi, shares=shares))
            print('N = %d Y = %d: E_Y/N = %.4f  R = %d (direct %d)  corr with f %.4f  corr with (u/3) %.4f' %
                  (N, Y, E / N, R, Rdir, rho_f, rho_chi))
            print('      energy share / seed share: ' + '  '.join('%s %.3f/%.3f' % (k, v[0], v[1]) for k, v in shares.items()),
                  flush=True)
    out['frontier_characters'] = res
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'preprint_connections_checks_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
