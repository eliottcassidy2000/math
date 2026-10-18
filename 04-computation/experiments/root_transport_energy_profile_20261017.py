#!/usr/bin/env python3
"""Energy profile of the exact ROOT transport over many N (S41, vectorised).

For each N: initial coefficients on primitive odd seeds (as in root_transport_energy_ledger_20261017.py), then repeated
odd transports s -> seed((3s+1)/2^a) with sign (-1)^a until every source sits at ROOT. Records, per N, the energy
profile E_d(N)/N, its maximum over d, the depth at which it is attained, the l1 profile, and the final R(N)^2/N. The
uniform-control target is sup_d E_d(N) <= C N, which would give |R(N)| <= sqrt(C N) for every N.
Usage: python root_transport_energy_profile_20261017.py KMIN KMAX [SAMPLES_PER_OCTAVE]
Session opus-2026-10-17-S41.
"""
import sys, os, json, math
import numpy as np


def initial(N):
    m = np.arange(1, 4 * N + 1, 2, dtype=np.int64)
    r = m % 3
    w = np.where(r == 2, 2, -1).astype(np.int64)
    chi = np.where(r == 0, 0, np.where(r == 1, 1, -1)).astype(np.int64)
    k = w + np.where(m < 2 * N, -3 * chi - w, 0)
    s = m.copy()
    mask = s % 8 == 5
    while mask.any():
        s[mask] = (s[mask] - 1) // 4
        mask = s % 8 == 5
    u, inv = np.unique(s, return_inverse=True)
    c = np.zeros(len(u), dtype=np.int64)
    np.add.at(c, inv, k)
    keep = c != 0
    return u[keep], c[keep]


def step(S, C):
    root = S == 1
    u = 3 * S + 1
    a = np.where(u % 4 == 0, 2, 1)
    u = np.where(root, 1, u >> a)
    sign = np.where(root, 1, np.where(a == 2, 1, -1))
    mask = u % 8 == 5
    while mask.any():
        u[mask] = (u[mask] - 1) // 4
        mask = u % 8 == 5
    v, inv = np.unique(u, return_inverse=True)
    c = np.zeros(len(v), dtype=np.int64)
    np.add.at(c, inv, sign * C)
    keep = c != 0
    return v[keep], c[keep]


def profile(N):
    S, C = initial(N)
    E = []; L1 = []
    while True:
        Cf = C.astype(np.float64)
        E.append(float(Cf @ Cf) / N); L1.append(float(np.abs(Cf).sum()) / N)
        if len(S) == 0 or (len(S) == 1 and S[0] == 1):
            break
        S, C = step(S, C)
    R = int(C[0]) if len(S) else 0
    return R, E, L1


def main():
    KMIN, KMAX = int(sys.argv[1]), int(sys.argv[2])
    SAMP = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    rng = np.random.default_rng(20261017)
    rows = []
    for k in range(KMIN, KMAX + 1):
        X = 2 ** k
        Ns = [X] + [int(x) for x in rng.integers(X, 2 * X, SAMP - 1)]
        for N in Ns:
            R, E, L1 = profile(N)
            dmax = int(np.argmax(E))
            rows.append(dict(N=N, R=R, R2_over_N=R * R / N, depth=len(E) - 1, Emax=max(E), d_at_max=dmax, E=E[:12],
                             E_tail_max=max(E[6:]) if len(E) > 6 else None))
            print('N = %9d (2^%d octave): R = %7d  R^2/N = %6.3f  depth %3d  max_d E_d/N = %6.3f at d = %3d  (E_0..E_6: %s)'
                  % (N, k, R, R * R / N, len(E) - 1, max(E), dmax, ' '.join('%.2f' % e for e in E[:7])), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(rows, open(os.path.join(here, '..', '..', '05-knowledge', 'results',
                                      'root_transport_energy_profile_%d_%d_20261017.json' % (KMIN, KMAX)), 'w'), indent=1)


if __name__ == '__main__':
    main()
