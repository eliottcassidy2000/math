#!/usr/bin/env python3
"""Cross-term decomposition of the frontier energy E_Y(N) of the stopped ROOT transport (S42).

The stopped process is run synchronously in depth: at each step every active seed s > Y moves to P(s) with its sign;
seeds <= Y are frozen. Contributions landing on one seed create cross terms of three kinds:

    MM_hi : two moving contributions meeting at a seed > Y        (the merges of the ideal regime, above the frontier)
    MM_lo : two moving contributions arriving at the same frozen seed <= Y in the same step
    MF    : a moving contribution arriving at a frozen seed that already carries current (arrival against the frontier)

and E_Y(N) = E_0(N) + 2 (MM_hi + MM_lo + MF) exactly (checked). The invariant sum_u f(u) c_Y(u) = R(N) is checked too.
This separates the part of the frontier energy that is a statement about the ideal regime (MM_hi) from the arrival
terms that only exist because the process is stopped.
Usage: python root_transport_frontier_decomposition_20261017.py "N1,N2,..." "Y1,Y2,..."
Session opus-2026-10-17-S42.
"""
import sys, os, json
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from root_transport_energy_profile_20261017 import initial


@njit(cache=True)
def move(S):
    n = S.shape[0]
    T = np.empty(n, dtype=np.int64)
    G = np.empty(n, dtype=np.int64)
    for i in range(n):
        u = 3 * S[i] + 1
        a = 0
        while u % 2 == 0:
            u //= 2
            a += 1
        while u % 8 == 5:
            u = (u - 1) // 4
        T[i] = u
        G[i] = 1 if a == 2 else -1
    return T, G


def parity(n):
    c = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        c += 1
    return 1 - 2 * (c & 1)


def run(N, Y):
    S, C = initial(N)
    E0 = int((C * C).sum())
    fz = np.zeros(Y + 1, dtype=np.int64)
    low = S <= Y
    np.add.at(fz, S[low], C[low])
    S, C = S[~low], C[~low]
    MMhi = MMlo = MF = 0
    rows = []
    d = 0
    while len(S):
        d += 1
        T, G = move(S)
        A = G * C
        hi = T > Y
        Th, Ah = T[hi], A[hi]
        v, inv = np.unique(Th, return_inverse=True)
        sums = np.zeros(len(v), dtype=np.int64); np.add.at(sums, inv, Ah)
        sq = np.zeros(len(v), dtype=np.int64); np.add.at(sq, inv, Ah * Ah)
        mmhi = int((sums * sums - sq).sum()) // 2
        Tl, Al = T[~hi], A[~hi]
        sl = np.zeros(Y + 1, dtype=np.int64); np.add.at(sl, Tl, Al)
        ql = np.zeros(Y + 1, dtype=np.int64); np.add.at(ql, Tl, Al * Al)
        mmlo = int((sl * sl - ql).sum()) // 2
        mf = int((fz * sl).sum())
        fz += sl
        MMhi += mmhi; MMlo += mmlo; MF += mf
        keep = sums != 0
        S, C = v[keep], sums[keep]
        rows.append((d, mmhi, mmlo, mf, int(len(S)), int(np.abs(Al).sum())))
    EY = int((fz * fz).sum())
    assert EY == E0 + 2 * (MMhi + MMlo + MF)
    R = sum(int(fz[u]) * parity(u) for u in range(1, Y + 1) if fz[u] != 0)
    return dict(N=N, Y=Y, E0=E0 / N, EY=EY / N, MMhi=2 * MMhi / N, MMlo=2 * MMlo / N, MF=2 * MF / N,
                depth=d, R2N=R * R / N, rows=rows)


def main():
    Ns = [int(a) for a in sys.argv[1].split(',')]
    Ys = [int(a) for a in sys.argv[2].split(',')]
    out = []
    for N in Ns:
        for Y in Ys:
            r = run(N, Y)
            out.append({k: v for k, v in r.items() if k != 'rows'})
            # depth profile of the three kinds, in thirds of the stopped depth
            rows = r['rows']; D = len(rows)
            parts = []
            for a, b in ((0, D // 3), (D // 3, 2 * D // 3), (2 * D // 3, D)):
                parts.append(tuple(2 * sum(x[i] for x in rows[a:b]) / N for i in (1, 2, 3)))
            print('N = %8d  Y = %6d: E_0/N = %.3f -> E_Y/N = %.3f = E_0 + MM_hi %+.3f + MM_lo %+.3f + MF %+.3f'
                  '  (depth %d; R^2/N = %.3f)' % (N, Y, r['E0'], r['EY'], r['MMhi'], r['MMlo'], r['MF'], D, r['R2N']))
            print('      by depth thirds (MM_hi, MM_lo, MF): ' + '  '.join('(%+.3f, %+.3f, %+.3f)' % p for p in parts),
                  flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results',
                                     'root_transport_frontier_decomposition_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
