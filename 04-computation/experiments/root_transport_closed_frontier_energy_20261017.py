#!/usr/bin/env python3
"""Frontier energy on clock_duality_oct15's certified forward-closed frontiers V (R6 of root_sparse_boundary_module_20261017i) (S42).

R6 (clock_duality): if V is finite, contains 1, P(V) is contained in V and a rank strictly descends to ROOT, then stopping on
first entry into V gives the same outside vector as the unstopped transport at every depth (no return kernel), and
R(N)^2 <= |V| E_V(N) with E_V(N) = sum_(u in V) c_V(u)^2. Their certified sets for heights Y = 100 and 1000 have 74 and
627 seeds (largest 2051 and 63665). This script computes E_V(N)/N for every N < 2^(KMAX+1), with the invariant
R(N) = sum_(u in V) f(u) c_V(u) checked against H(4N) - H(N) at every N, and per-octave statistics, for the closed sets
in the receipt 04-computation/experiments/energy_extraction_20261017i/frontier_renewal.json.
Usage: python root_transport_closed_frontier_energy_20261017.py KMAX
Session opus-2026-10-17-S42.
"""
import sys, os, json, time
import numpy as np
from numba import njit, prange
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from root_transport_frontier_energy_big_20261017 import clock_parity8, seed, wfun, thr


@njit(parallel=True, cache=True)
def closed_frontier_map(M, inV, maxV):
    n = (M + 1) // 2
    fu = np.zeros(n, dtype=np.int32)
    fs = np.zeros(n, dtype=np.int8)
    for i in prange(n):
        s = 2 * i + 1
        if s % 8 == 5:
            continue
        sign = 1
        x = s
        while not (x <= maxV and inV[x]):
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
def run(Nmax, maxV, fu, fs, par):
    c = np.zeros(maxV + 1, dtype=np.int64)
    E = 0
    R = 0
    H4 = 0
    H1 = 0
    ok = True
    K = 64
    emax = np.zeros(K); eargmax = np.zeros(K, dtype=np.int64); esum = np.zeros(K); emin = np.full(K, 1e300)
    cnt = np.zeros(K, dtype=np.int64)
    for N in range(1, Nmax):
        for t in range(3):
            if t == 0:
                m = 4 * N - 3; amt = wfun(m)
            elif t == 1:
                m = 4 * N - 1; amt = wfun(m)
            else:
                m = 2 * N - 1; amt = thr(m)
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
        if e > emax[k]:
            emax[k] = e; eargmax[k] = N
        if e < emin[k]:
            emin[k] = e
        esum[k] += e
        cnt[k] += 1
    return ok, emax, eargmax, esum, emin, cnt


def main():
    KMAX = int(sys.argv[1])
    here = os.path.dirname(os.path.abspath(__file__))
    rec = json.load(open(os.path.join(here, 'energy_extraction_20261017i', 'frontier_renewal.json')))
    Nmax = 2 ** (KMAX + 1); M = 4 * Nmax
    t0 = time.time()
    par = clock_parity8(M)
    print('clock parities to %d [%.0f s]' % (M, time.time() - t0), flush=True)
    res = {}
    for cf in rec['closed_frontiers']:
        if cf['Y'] < 100:
            continue
        V = [v[0] for v in cf['vertices']]
        maxV = max(V)
        inV = np.zeros(maxV + 1, dtype=np.bool_)
        inV[V] = True
        fu, fs = closed_frontier_map(M, inV, maxV)
        ok, emax, eargmax, esum, emin, cnt = run(Nmax, maxV, fu, fs, par)
        rows = []
        for k in range(8, KMAX + 1):
            rows.append(dict(k=k, EV_max=float(emax[k]), EV_argmax=int(eargmax[k]), EV_mean=float(esum[k] / cnt[k]),
                             EV_min=float(emin[k])))
        res[cf['Y']] = dict(size=len(V), max_seed=maxV, invariant_ok=bool(ok), octaves=rows)
        print('closed frontier from Y = %d: |V| = %d, max seed %d (R(N) invariant at every N: %s) [%.0f s]'
              % (cf['Y'], len(V), maxV, ok, time.time() - t0), flush=True)
        for r in rows:
            print('   octave 2^%2d: E_V/N mean %.4f  max %.4f (N = %d)  min %.4f'
                  % (r['k'], r['EV_mean'], r['EV_max'], r['EV_argmax'], r['EV_min']), flush=True)
    json.dump({str(k): v for k, v in res.items()}, open(os.path.join(here, '..', '..', '05-knowledge', 'results',
              'root_transport_closed_frontier_energy_%d_20261017.json' % KMAX), 'w'), indent=1)


if __name__ == '__main__':
    main()
