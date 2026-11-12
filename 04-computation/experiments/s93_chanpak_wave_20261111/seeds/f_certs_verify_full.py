#!/usr/bin/env python3
"""S93 lane seeds, FINISHER: independent recount of the completed certificate files certs_plus_full.txt and
certs_minus_full.txt (written by f_cert_coverage.py write, after the interrupted re-run of b_census.py was completed
with the sub-parts b_census_part_P_32.npz, P = 4..7 mod 8).

Each line: c v order colcode Omask M  (a pair (B, O), offset +1 or -1, slope M = c 2^v).
B: vertices 0..n-1, column-major code (pairs (i, j), i < j, ordered by j then i; bit = 1 iff i -> j).
Closed tournament T_t (re-derived here from the definition, no code shared with the census or with seedlib):
vertices of B, then a transitive tail x_1 -> ... -> x_t, then the closer w;  every vertex of B beats every x_k;
every x_k beats w;  w beats the vertices of O and loses to the other vertices of B.
Check: the number of Hamiltonian paths of T_t (subset recursion by end vertex) equals sign + 2^t M for t = 1 and t = 2.
"""
import sys

import numpy as np
from numba import njit, prange


@njit(cache=True)
def count_closed(n, code, O, t, dp):
    N = n + t + 1
    out = np.zeros(N, np.int64)
    b = 0
    for j in range(1, n):
        for i in range(j):
            if (code >> b) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            b += 1
    w = n + t
    for x in range(n):
        for k in range(t):
            out[x] |= 1 << (n + k)
        if (O >> x) & 1:
            out[w] |= 1 << x
        else:
            out[x] |= 1 << w
    for k in range(t):
        for l in range(k + 1, t):
            out[n + k] |= 1 << (n + l)
        out[n + k] |= 1 << w
    full = (1 << N) - 1
    for S in range(full + 1):
        for v in range(N):
            dp[S, v] = 0
    for v in range(N):
        dp[1 << v, v] = 1
    for S in range(1, full + 1):
        for v in range(N):
            c = dp[S, v]
            if c:
                rest = out[v] & ~S & full
                for u in range(N):
                    if (rest >> u) & 1:
                        dp[S | (1 << u), u] += c
    tot = 0
    for v in range(N):
        tot += dp[full, v]
    return tot


@njit(parallel=True, cache=True)
def verify(rows, sign, bad):
    nthr = 16
    m = rows.shape[0]
    for ch in prange(nthr):
        dp = np.zeros((1 << 13, 13), np.int64)
        for idx in range(ch, m, nthr):
            c, v, n, code, O, M = rows[idx, 0], rows[idx, 1], rows[idx, 2], rows[idx, 3], rows[idx, 4], rows[idx, 5]
            ok = (M == (c << v)) and (c & 1) == 1 and O > 0 and O < (1 << n)
            for t in (1, 2):
                if count_closed(n, code, O, t, dp) != sign + (1 << t) * M:
                    ok = False
            if not ok:
                bad[ch] += 1


def main():
    for sign, fn, gap in ((1, "certs_plus_full.txt", 53559), (-1, "certs_minus_full.txt", 59083)):
        rows = np.array([[int(x) for x in line.split()] for line in open(fn) if not line.startswith("#")], dtype=np.int64)
        assert rows[:, 2].max() <= 10
        bad = np.zeros(16, np.int64)
        verify(rows, sign, bad)
        have = set(rows[:, 0].tolist())
        first = min(c for c in range(1, 65537, 2) if c not in have)
        below = rows[rows[:, 0] < gap]
        print(f"sign {sign:+d}: {rows.shape[0]} certificates in {fn}; closed tournaments recounted: {2 * rows.shape[0]} (t = 1, 2); certificates failing: {int(bad.sum())}")
        print(f"   first odd c without a certificate: {first} (census gap: {gap}); odd c < {gap} certified: {below.shape[0]} of {gap // 2}; largest v among them: {int(below[:, 1].max())}; "
              f"orders used: {sorted(set(rows[:, 2].tolist()))}; certificates with v = 0 and c <= 32767: {int(((rows[:, 1] == 0) & (rows[:, 0] <= 32767)).sum())}")
        if bad.sum():
            sys.exit(1)


if __name__ == "__main__":
    main()
