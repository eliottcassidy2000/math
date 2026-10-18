#!/usr/bin/env python3
"""Per-basin square-root test: max over each octave of R_u(N)^2/N for every frontier seed u (S42).

By the basin form (root_transport_basin_currents_check_20261018.py), c_Y(u) = f(u) R_u(N), where R_u is the ROOT current
restricted to the first-entry basin of u. HYP-9290 asks sum_u R_u(N)^2 <= C N. This script asks the finer question
whether every single basin current is O(sqrt(N)): for each u <= Y and each octave 2^k <= N < 2^(k+1), it records
max R_u(N)^2/N. (A maximum of c[u]^2/N over N is attained at an N where c[u] changes, so it suffices to evaluate at the
three updated entries per N.) It prints, per octave, the largest per-basin maximum and the seed attaining it, and the
number of seeds whose octave maximum exceeds 1.
Usage: python root_transport_basin_maxima_20261018.py KMAX Y
Session opus-2026-10-17-S42.
"""
import sys, os, json, time
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from root_transport_frontier_energy_big_20261017 import clock_parity8, frontier_map, seed, wfun, thr


@njit(cache=True)
def run(Nmax, Y, fu, fs):
    c = np.zeros(Y + 1, dtype=np.int64)
    K = 64
    mx = np.zeros((K, Y + 1))
    for N in range(1, Nmax):
        k = 0
        x = N
        while x > 1:
            x >>= 1
            k += 1
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
            c[u] += fs[i] * amt
            v = c[u] * c[u] / N
            if v > mx[k, u]:
                mx[k, u] = v
        if N > 1 and (N & (N - 1)) == 0:
            # first N of a new octave: record every coordinate once (covers seeds not updated in the octave)
            for u in range(Y + 1):
                v = c[u] * c[u] / N
                if v > mx[k, u]:
                    mx[k, u] = v
    return mx


def main():
    KMAX = int(sys.argv[1]); Y = int(sys.argv[2])
    Nmax = 2 ** (KMAX + 1); M = 4 * Nmax
    fu, fs = frontier_map(M, Y)
    mx = run(Nmax, Y, fu, fs)
    seeds = [u for u in range(1, Y + 1) if u % 2 == 1 and u % 8 != 5]
    rows = []
    for k in range(10, KMAX + 1):
        vals = [(mx[k, u], u) for u in seeds]
        top = max(vals)
        rows.append(dict(k=k, max_basin=float(top[0]), argmax_seed=int(top[1]), n_above_1=int(sum(1 for v, u in vals if v > 1)),
                         top5=[[int(u), float(v)] for v, u in sorted(vals, reverse=True)[:5]]))
        print('octave 2^%2d: largest per-basin max R_u^2/N = %.3f (seed %d); seeds with max > 1: %d; top: %s'
              % (k, top[0], top[1], rows[-1]['n_above_1'], ', '.join('%d:%.2f' % (u, v) for v, u in sorted(vals, reverse=True)[:5])),
              flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(rows, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'root_transport_basin_maxima_%d_%d_20261018.json' % (KMAX, Y)), 'w'), indent=1)


if __name__ == '__main__':
    main()
