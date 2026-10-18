#!/usr/bin/env python3
"""S41 audit A12: own direct transport at the committed octave maximisers beyond the range of A3/A11
(N = 6904411 for the octave 2^22, sup_d E_d/N = 6.2698 claimed; N = 2109637 for 2^21, 4.8572 claimed), recording
sup_d E_d(N)/N, its depth, R(N)^2/N and the depth to ROOT. R(N) is also computed directly as H(4N) - H(N) from own
clock parities (numba).
Usage: python a12_direct_large_N.py [N ...]
"""
import os, sys, time, json
sys.dont_write_bytecode = True
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(os.environ.get('TEMP', '/tmp'), 's41_audit_numba_cache'))
import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))


def seed_vec(u):
    u = u.copy()
    m = (u % 8) == 5
    while m.any():
        u[m] = (u[m] - 1) // 4
        m = (u % 8) == 5
    return u


def aggregate(keys, vals):
    order = np.argsort(keys, kind='stable')
    k = keys[order]
    v = vals[order]
    starts = np.flatnonzero(np.r_[True, k[1:] != k[:-1]])
    s = np.add.reduceat(v, starts)
    keep = s != 0
    return k[starts][keep], s[keep]


@njit(cache=True)
def current(N):
    """R(N) = sum_(N < n <= 4N) w(n) (-1)^clock(n), clock by direct iteration (first hit of 1)."""
    tot = 0
    for n in range(N + 1, 4 * N + 1):
        x = n
        c = 0
        while x != 1:
            if x % 2 == 0:
                x //= 2
            else:
                x = (3 * x + 1) // 2
            c += 1
        w = 2 if n % 3 == 2 else -1
        tot += w if c % 2 == 0 else -w
    return tot


def main():
    Ns = [int(a) for a in sys.argv[1:]] or [2109637, 6904411]
    out = {}
    for N in Ns:
        t0 = time.time()
        m = np.arange(1, 4 * N + 1, 2, dtype=np.int64)
        r = m % 3
        k = np.where(m < 2 * N, -3 * np.where(r == 0, 0, np.where(r == 1, 1, -1)), np.where(r == 2, 2, -1)).astype(np.int64)
        S, C = aggregate(seed_vec(m), k)
        del m, r, k
        E = []
        while True:
            E.append(int(np.dot(C, C)))
            if len(S) == 1 and S[0] == 1:
                break
            y = 3 * S + 1
            a = np.where((y & (-y)) == 4, 2, 1)
            S, C = aggregate(seed_vec(y >> a), np.where(a == 2, C, -C))
        R = int(C[0])
        Rd = int(current(N))
        sup = max(E) / N
        out[N] = dict(sup=sup, d_sup=int(np.argmax(E)), R=R, R_direct=Rd, R2N=R * R / N, depth=len(E) - 1, E0=E[0] / N)
        print('N = %d: sup_d E_d/N = %.6f at d = %d; E_0/N = %.4f; R = %d (direct H(4N)-H(N): %d); R^2/N = %.4f; depth %d [%.0f s]'
              % (N, sup, out[N]['d_sup'], E[0] / N, R, Rd, R * R / N, len(E) - 1, time.time() - t0), flush=True)
    json.dump(out, open(os.path.join(HERE, 'a12_direct_large_N.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
