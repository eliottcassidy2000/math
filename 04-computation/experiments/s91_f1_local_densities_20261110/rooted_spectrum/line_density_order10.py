#!/usr/bin/env python3
"""S91 (HYP-9661 (A), pre-registered): how many one-vertex sources lie on each line D + 2U = N.

For odd N let r_m(N) be the number of odd D <= N/3 such that (D, (N - D)/2) is a pair (H(B), U_B(y) 2^j) for a rooted
tournament (B, y) of order <= m and some j >= 0 (j sinks appended to B double the slope and keep the offset).  Then
N = H(B') + 2 U_{B'}(y) is the Hamiltonian path count of the two-step tournament of THM-4865 / the near-join of
HYP-9453.  A line has about N/6 lattice points (D odd, D <= N/3, since U >= D).

Every rooted tournament of order n + 1 is (base of order n, cut of the new vertex) with the new vertex as root; the
bases are the class representatives reps_le8.npz (orders 1..8) and the list codes9.npy (labelled tournaments of
order 9 containing every isomorphism class; written by census.py, not kept in the repository).

Usage: python3 line_density_order10.py [MAXBASE]     (MAXBASE = 8: rooted order <= 9; 9: rooted order <= 10)
Output: the odd N < 2^15 with r(N) = 0; minimum, median and mean of r on dyadic blocks; the test r(N) >= N/64 on
[129, 2^14); the share of lattice points that are sources by decile of 3D/N.  Pre-registered predictions:
PREREGISTRATION.md, section "Line density at rooted order 10".
"""
import os
import sys
import time

import numpy as np
from numba import njit, prange

HERE = os.path.dirname(os.path.abspath(__file__))
NMAX = 1 << 15


@njit(parallel=True, cache=True)
def scan(n, reps, M, nmax):
    full = (1 << n) - 1
    for r in prange(reps.shape[0]):
        f = np.zeros((1 << n, n), np.int64)
        g = np.zeros((1 << n, n), np.int64)
        hS = np.zeros(1 << n, np.int64)
        out = np.zeros(n, np.int64)
        Q = np.zeros((n, n), np.int64)
        ym = np.zeros(n, np.int64)
        code = reps[r]
        for j in range(1, n):
            base = j * (j - 1) // 2
            for i in range(j):
                if (code >> (base + i)) & 1:
                    out[i] |= (1 << j)
                else:
                    out[j] |= (1 << i)
        hS[0] = 1
        for S in range(1, full + 1):
            tot = 0
            for v in range(n):
                if ((S >> v) & 1) == 0:
                    continue
                R = S ^ (1 << v)
                if R == 0:
                    f[S, v] = 1
                    g[S, v] = 1
                    tot += 1
                    continue
                a = 0
                b = 0
                for w in range(n):
                    if (R >> w) & 1:
                        if (out[w] >> v) & 1:
                            a += f[R, w]
                        else:
                            b += g[R, w]
                f[S, v] = a
                g[S, v] = b
                tot += a
            hS[S] = tot
        for S in range(1, full + 1):
            C = full ^ S
            for a in range(n):
                if (S >> a) & 1:
                    ym[a] += g[S, a] * hS[C]
                    if f[S, a]:
                        for b in range(n):
                            if (C >> b) & 1:
                                Q[a, b] += f[S, a] * g[C, b]
        HB = hS[full]
        for cut in range(1 << n):
            H = 0
            U = HB
            for b in range(n):
                if (cut >> b) & 1:
                    H += g[full, b]
                    U += ym[b]
                else:
                    H += f[full, b]
                    for c in range(n):
                        if (cut >> c) & 1:
                            H += Q[b, c]
            N = H + 2 * U
            while N < nmax:
                M[N >> 1, H >> 1] = 1
                U *= 2
                N = H + 2 * U


def main():
    maxbase = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    t0 = time.time()
    M = np.zeros((NMAX >> 1, (NMAX // 3 >> 1) + 2), dtype=np.uint8)
    # rooted order 1: (D, U) = (1, 1)
    U = 1
    while 1 + 2 * U < NMAX:
        M[(1 + 2 * U) >> 1, 0] = 1
        U *= 2
    z = np.load(os.path.join(HERE, "reps_le8.npz"))
    for n in range(1, min(maxbase, 8) + 1):
        scan(n, z["n%d" % n].astype(np.int64), M, NMAX)
    if maxbase >= 9:
        codes9 = np.load(os.path.join(HERE, "codes9.npy")).astype(np.int64)
        scan(9, codes9, M, NMAX)
        print(f"bases of order 9: {len(codes9)} labelled tournaments (every class)  ({time.time() - t0:.0f} s)")
    m = maxbase + 1
    r = M.sum(axis=1).astype(np.int64)          # r[N >> 1] for odd N
    odd = np.arange(1, NMAX, 2)
    rr = r[odd >> 1]
    miss = odd[rr == 0]
    print(f"rooted order <= {m}: odd N < {NMAX} with r(N) = 0: {miss.tolist()}")
    ok = miss.tolist() == [1, 7, 11, 13, 21, 47, 67]
    print("block: minimum of r (at N), median, mean, mean share of the N/6 lattice points")
    for k in range(6, 15):
        sel = (odd >= (1 << k)) & (odd < (1 << (k + 1)))
        c = rr[sel]
        nn = odd[sel]
        print(f"  [2^{k}, 2^{k + 1}): min {c.min():5d} (N = {nn[c.argmin()]}), median {int(np.median(c)):5d}, "
              f"mean {c.mean():8.1f}, share {np.mean(c / (nn / 6.0)):.3f}")
    sel = (odd >= 129) & (odd < (1 << 14))
    viol = odd[sel][rr[sel] * 64 < odd[sel]]
    print(f"r(N) >= N/64 for every odd N in [129, 2^14): {'yes' if len(viol) == 0 else 'NO'}; violations {len(viol)}"
          + (f", first {viol[:8].tolist()}" if len(viol) else ""))
    sel = (odd >= (1 << 12)) & (odd < (1 << 13))
    print(f"mean of r on [2^12, 2^13): {rr[sel].mean():.1f}")
    for lo, hi in ((1 << 10, 1 << 11), (1 << 12, 1 << 13)):
        num = np.zeros((2, 10))
        den = np.zeros((2, 10))
        for N in range(lo + 1, hi, 2):
            D = np.arange(1, N // 3 + 1, 2)
            Uv = (N - D) // 2
            k = np.minimum(9, (30 * D) // N)
            par = Uv & 1
            hit = M[N >> 1, D >> 1]
            np.add.at(den, (par, k), 1)
            np.add.at(num, (par, k), hit)
        print(f"N in [{lo}, {hi}): share of lattice points (D, U) that are sources, by decile of 3D/N")
        for par in (1, 0):
            print(f"   U {'odd ' if par else 'even'}: " + " ".join(f"{num[par, k] / max(den[par, k], 1):.2f}" for k in range(10))
                  + f"   overall {num[par].sum() / den[par].sum():.3f}")
    print(f"({time.time() - t0:.0f} s)")
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
