#!/usr/bin/env python3
"""Incremental energy profiles of the exact ROOT transport for EVERY N below 2^(KMAX+1) (S41).

The initial source vector changes by three entries when N -> N+1 (odd m = 4N+1, 4N+3 enter the w-term and m = 2N+1
the threshold term of clock_duality_oct15's fold R1), and the odd transport maps a point mass at a seed to a signed
point mass. So the transported vectors c_d(N), d = 0, 1, ..., can be maintained for all depths at once: each update
walks one path of the odd-transport tree (nodes = (depth, seed)), and the energies E_d(N) = ||c_d(N)||^2 are updated
in O(1) per node. Sources that reach ROOT stay there; ROOT is kept as a separate column per depth.

Outputs, for every N: sup_d E_d(N)/N, the depth where it is attained, the final R(N)^2/N (checked against the direct
current H(4N) - H(N) from the clock parities), and the climb factor max_d E_d / min_(d' <= d) E_d'. Octave traces
T_d(X) = sum_(X <= N < 2X) E_d(N) for every octave X = 2^k, k <= KMAX, are accumulated exactly.
Usage: python root_transport_incremental_trace_20261017.py KMAX
Session opus-2026-10-17-S41.
"""
import sys, os, json, time
import numpy as np
from numba import njit


def seed_np(u):
    u = u.copy()
    mask = u % 8 == 5
    while mask.any():
        u[mask] = (u[mask] - 1) // 4
        mask = u % 8 == 5
    return u


def build_tree(M):
    """Levels of the odd-transport tree for all seeds s in [3, M]; returns succ (global id or -1 = ROOT), sgn, sizes."""
    lvl = np.arange(3, M + 1, 2, dtype=np.int64)
    lvl = lvl[lvl % 8 != 5]
    level0 = lvl
    succ_parts, sgn_parts, sizes = [], [], []
    offset = 0
    while len(lvl):
        u = 3 * lvl + 1
        a = np.where(u % 4 == 0, 2, 1)
        u = seed_np(u >> a)
        sg = np.where(a == 2, 1, -1).astype(np.int8)
        nxt = np.unique(u[u != 1])
        loc = np.searchsorted(nxt, u)
        succ = np.where(u == 1, -1, offset + len(lvl) + loc).astype(np.int64)
        succ_parts.append(succ); sgn_parts.append(sg); sizes.append(len(lvl))
        offset += len(lvl)
        lvl = nxt
    return level0, np.concatenate(succ_parts), np.concatenate(sgn_parts), sizes


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
def apply(s, amt, level0, succ, sgn, V, Vroot, E, DMAX):
    d = 0; sign = 1
    if s != 1:
        node = np.searchsorted(level0, s)
        while node >= 0:
            delta = sign * amt
            old = V[node]
            E[d] += 2 * old * delta + delta * delta
            V[node] = old + delta
            sign *= sgn[node]
            node = succ[node]
            d += 1
    delta = sign * amt
    for dd in range(d, DMAX):
        old = Vroot[dd]
        E[dd] += 2 * old * delta + delta * delta
        Vroot[dd] = old + delta


@njit(cache=True)
def run(Nmax, KMAX, level0, succ, sgn, DMAX):
    V = np.zeros(len(succ), dtype=np.int64)
    Vroot = np.zeros(DMAX, dtype=np.int64)
    E = np.zeros(DMAX, dtype=np.int64)
    T = np.zeros((KMAX + 1, DMAX), dtype=np.float64)
    supE = np.zeros(Nmax, dtype=np.float64)
    dsup = np.zeros(Nmax, dtype=np.int32)
    climb = np.zeros(Nmax, dtype=np.float64)
    E6 = np.zeros(Nmax, dtype=np.float64)
    Rf = np.zeros(Nmax, dtype=np.int64)
    k = 0
    for N in range(1, Nmax):
        m1 = 4 * N - 3; m2 = 4 * N - 1; m3 = 2 * N - 1
        for m, amt in ((m1, wfun(m1)), (m2, wfun(m2)), (m3, thr(m3))):
            s = m
            while s % 8 == 5:
                s = (s - 1) // 4
            apply(s, amt, level0, succ, sgn, V, Vroot, E, DMAX)
        while (1 << (k + 1)) <= N:
            k += 1
        if k <= KMAX:
            for d in range(DMAX):
                T[k, d] += E[d]
        best = E[0]; bd = 0; mn = E[0]; cl = 1.0
        for d in range(DMAX):
            e = E[d]
            if e > best:
                best = e; bd = d
            if e < mn:
                mn = e
            if mn > 0:
                r = e / mn
                if r > cl:
                    cl = r
        supE[N] = best / N; dsup[N] = bd; climb[N] = cl
        E6[N] = E[6] / N if DMAX > 6 else 0.0
        Rf[N] = Vroot[DMAX - 1]
    return T, supE, dsup, climb, E6, Rf


def main():
    KMAX = int(sys.argv[1])
    Nmax = 2 ** (KMAX + 1)
    M = 4 * Nmax
    t0 = time.time()
    level0, succ, sgn, sizes = build_tree(M)
    DMAX = len(sizes) + 1
    print('tree: %d sources, %d nodes, %d levels (DMAX %d) [%.0f s]' % (len(level0), len(succ), len(sizes), DMAX,
                                                                         time.time() - t0), flush=True)
    clk = clock_parity(M)
    f = 1 - 2 * (clk & 1)
    n = np.arange(M + 1)
    wv = np.where(n % 3 == 2, 2, -1); wv[0] = 0
    H = np.cumsum(wv * f)
    print('clock parities done [%.0f s]' % (time.time() - t0), flush=True)
    T, supE, dsup, climb, E6, Rf = run(Nmax, KMAX, level0, succ, sgn, DMAX)
    print('transport run done [%.0f s]' % (time.time() - t0), flush=True)
    Ns = np.arange(1, Nmax)
    Rdir = H[4 * Ns] - H[Ns]
    ok = np.array_equal(Rdir, Rf[1:])
    print('invariant check R(N) = H(4N) - H(N) for all N < 2^%d: %s' % (KMAX + 1, ok), flush=True)
    rows = []
    for kk in range(1, KMAX + 1):
        X = 2 ** kk
        prof = T[kk] / X ** 2
        sl = slice(X, 2 * X)
        DX = float(np.sum(Rdir[X - 1:2 * X - 1].astype(np.float64) ** 2))
        imin = int(np.argmin(prof))
        row = dict(k=kk, T0=prof[0], T6=prof[6] if len(prof) > 6 else None, Tmax=float(prof.max()),
                   d_Tmax=int(prof.argmax()), Tmin=float(prof.min()), d_Tmin=imin, Tfinal=float(prof[-1]),
                   DX=DX / X ** 2, sup_pointwise=float(supE[sl].max()), N_at_sup=int(X + np.argmax(supE[sl])),
                   frac_sup_at_0=float(np.mean(dsup[sl] == 0)), climb_max=float(climb[sl].max()),
                   climb_mean=float(climb[sl].mean()), E6_mean=float(E6[sl].mean()),
                   R2N_max=float((Rdir[X - 1:2 * X - 1].astype(np.float64) ** 2 / Ns[X - 1:2 * X - 1]).max()),
                   profile=prof[:min(len(prof), 260)].tolist())
        rows.append(row)
        print('X = 2^%2d: T_0 %.4f  T_6 %.4f  min %.4f (d %3d)  max %.4f (d %3d)  final %.4f = D_X %.4f | pointwise sup_d E_d/N max '
              '%.3f (N = %d), at d = 0: %.3f; climb max %.2f mean %.3f; E_6/N mean %.4f; max R^2/N %.3f'
              % (kk, row['T0'], row['T6'] or 0, row['Tmin'], imin, row['Tmax'], row['d_Tmax'], row['Tfinal'], row['DX'],
                 row['sup_pointwise'], row['N_at_sup'], row['frac_sup_at_0'], row['climb_max'], row['climb_mean'],
                 row['E6_mean'], row['R2N_max']), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    res = os.path.join(here, '..', '..', '05-knowledge', 'results')
    json.dump(dict(KMAX=KMAX, ok=bool(ok), DMAX=DMAX, octaves=rows), open(os.path.join(
        res, 'root_transport_incremental_trace_%d_20261017.json' % KMAX), 'w'), indent=1)
    np.savez_compressed(os.path.join(res, 'root_transport_incremental_trace_%d_20261017.npz' % KMAX),
                        supE=supE.astype(np.float32), dsup=dsup.astype(np.int16), climb=climb.astype(np.float32),
                        R=Rf.astype(np.int32))


if __name__ == '__main__':
    main()
