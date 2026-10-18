#!/usr/bin/env python3
"""S42 audit a2: independent re-implementation of the frontier-energy octave scan (no session code imported).

For every N < 2^(KMAX+1) and each Y: frontier coefficients c_Y updated by the three source entries of N
(odd m = 4N-3, 4N-1 with weight w(m), and m = 2N-1 with weight -3 chi(m) - w(m)), each pushed to its first seed <= Y
along the odd transport (own frontier map). Octave statistics of E_Y(N)/N (max, argmax, mean, min) and R(N)^2/N
(max, argmax); R(N) = sum_u f(u) c_Y(u) checked against H(4N) - H(N) at every N (own clock-parity table).
The three-entry update is re-derived in REPORT.md (K_N(m) = w(m)[m <= 4N] + (-3 chi(m) - w(m))[m < 2N]); a1 checked it
against the brute-force definition.
Usage: python a2_independent_incremental.py KMAX Y1 Y2 ...
"""
import os, sys, json, time
import numpy as np
from numba import njit, prange

HERE = os.path.dirname(os.path.abspath(__file__))


@njit(cache=False)
def parity_table(M):
    # par[n] = clock(n) mod 2, clock = first-hit shortcut steps to 1; par[1] = 0.
    par = np.zeros(M + 1, dtype=np.int8)
    for n in range(2, M + 1):
        if n % 2 == 0:
            par[n] = par[n // 2] ^ 1
        else:
            x = n
            st = 0
            while x >= n:
                if x % 2 == 0:
                    x = x // 2
                else:
                    x = (3 * x + 1) // 2
                st += 1
            par[n] = (st + par[x]) & 1
    return par


@njit(parallel=True, cache=False)
def fmap(M, Y):
    n_odd = (M + 1) // 2
    dest = np.zeros(n_odd, dtype=np.int32)
    sgn = np.zeros(n_odd, dtype=np.int8)
    for i in prange(n_odd):
        x = 2 * i + 1
        while x % 8 == 5:          # every odd m is mapped through its seed (seed preserves f)
            x = (x - 1) // 4
        s = 1
        while x > Y:
            v = 3 * x + 1
            a = 0
            while v % 2 == 0:
                v = v // 2
                a += 1
            if a & 1:
                s = -s
            while v % 8 == 5:
                v = (v - 1) // 4
            x = v
        dest[i] = x
        sgn[i] = s
    return dest, sgn


@njit(cache=False)
def scan(Nmax, Y, dest, sgn, par):
    c = np.zeros(Y + 1, dtype=np.int64)
    E = np.int64(0)
    R = np.int64(0)
    H4 = np.int64(0)
    H1 = np.int64(0)
    bad = 0
    K = 40
    emax = np.zeros(K); earg = np.zeros(K, dtype=np.int64); esum = np.zeros(K); emin = np.full(K, 1e300)
    rmax = np.zeros(K); rarg = np.zeros(K, dtype=np.int64); cnt = np.zeros(K, dtype=np.int64)
    for N in range(1, Nmax):
        for t in range(3):
            if t == 0:
                m = 4 * N - 3
            elif t == 1:
                m = 4 * N - 1
            else:
                m = 2 * N - 1
            r = m % 3
            w = 2 if r == 2 else -1
            if t == 2:
                chi = 0 if r == 0 else (1 if r == 1 else -1)
                amt = -3 * chi - w
            else:
                amt = w
            i = (m - 1) // 2
            u = dest[i]
            dlt = sgn[i] * amt
            E += dlt * (2 * c[u] + dlt)
            c[u] += dlt
            R += (1 - 2 * par[u]) * dlt
        for kk in range(4 * N - 3, 4 * N + 1):
            H4 += (2 if kk % 3 == 2 else -1) * (1 - 2 * par[kk])
        H1 += (2 if N % 3 == 2 else -1) * (1 - 2 * par[N])
        if R != H4 - H1:
            bad += 1
        k = 63 - 0
        # floor(log2 N)
        k = 0
        x = N
        while x > 1:
            x >>= 1
            k += 1
        e = E / N
        if e > emax[k]:
            emax[k] = e
            earg[k] = N
        if e < emin[k]:
            emin[k] = e
        esum[k] += e
        cnt[k] += 1
        r2 = (R * R) / N
        if r2 > rmax[k]:
            rmax[k] = r2
            rarg[k] = N
    return bad, emax, earg, esum, emin, rmax, rarg, cnt


def main():
    KMAX = int(sys.argv[1]); Ys = [int(a) for a in sys.argv[2:]]
    Nmax = 2 ** (KMAX + 1); M = 4 * Nmax
    t0 = time.time()
    par = parity_table(M)
    print('parity table to %d [%.0f s]' % (M, time.time() - t0), flush=True)
    out = {}
    for Y in Ys:
        dest, sgn = fmap(M, Y)
        bad, emax, earg, esum, emin, rmax, rarg, cnt = scan(Nmax, Y, dest, sgn, par)
        rows = []
        for k in range(8, KMAX + 1):
            rows.append(dict(k=k, max=float(emax[k]), argmax=int(earg[k]), mean=float(esum[k] / cnt[k]), min=float(emin[k]),
                             R2N_max=float(rmax[k]), R2N_argmax=int(rarg[k])))
        out[Y] = dict(invariant_failures=int(bad), octaves=rows)
        print('Y = %d: invariant failures %d [%.0f s]' % (Y, bad, time.time() - t0), flush=True)
        for r in rows:
            print('   2^%2d: max %.6f (N=%d) mean %.6f min %.6f | R^2/N max %.4f (N=%d)'
                  % (r['k'], r['max'], r['argmax'], r['mean'], r['min'], r['R2N_max'], r['R2N_argmax']), flush=True)
        del dest, sgn
    json.dump({str(k): v for k, v in out.items()}, open(os.path.join(HERE, 'a2_independent_incremental_%d.json' % KMAX), 'w'),
              indent=1)


if __name__ == '__main__':
    main()
