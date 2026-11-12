#!/usr/bin/env python3
"""Verifier library (independent of the lane's code): H tables, strong test by reachability, f/g/Q statistics."""
import numpy as np
from numba import njit, prange, set_num_threads
set_num_threads(4)


def pairs_of(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def pair_arrays(n):
    p = pairs_of(n)
    return np.array([a for a, b in p], np.int64), np.array([b for a, b in p], np.int64)


@njit(cache=True)
def adj_from_mask(n, m, pi, pj, adj, radj):
    for v in range(n):
        adj[v] = 0
        radj[v] = 0
    na = n * (n - 1) // 2
    for k in range(na):
        if (m >> k) & 1:
            a = pi[k]; b = pj[k]
        else:
            a = pj[k]; b = pi[k]
        adj[a] |= 1 << b
        radj[b] |= 1 << a


@njit(cache=True)
def hp(n, adj, dp):
    """number of Hamiltonian paths; dp[S, v] = paths covering S ending at v (pull form, different from the lane's push)."""
    full = (1 << n) - 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S >> v) & 1:
                dp[S, v] = 0
                continue
            R = S ^ (1 << v)
            if R == 0:
                dp[S, v] = 1
                continue
            t = 0
            for w in range(n):
                if (R >> w) & 1 and (adj[w] >> v) & 1:
                    t += dp[R, w]
            dp[S, v] = t
    tot = 0
    for v in range(n):
        tot += dp[full, v]
    return tot


@njit(cache=True)
def is_strong(n, adj, radj):
    full = (1 << n) - 1
    reach = 1
    ch = True
    while ch:
        ch = False
        for v in range(n):
            if (reach >> v) & 1:
                nr = reach | adj[v]
                if nr != reach:
                    reach = nr; ch = True
    if reach != full:
        return 0
    back = 1
    ch = True
    while ch:
        ch = False
        for v in range(n):
            if (back >> v) & 1:
                nr = back | radj[v]
                if nr != back:
                    back = nr; ch = True
    if back != full:
        return 0
    return 1


@njit(parallel=True, cache=True)
def build_table(n, pi, pj):
    na = n * (n - 1) // 2
    N = 1 << na
    H = np.zeros(N, np.int32)
    S = np.zeros(N, np.uint8)
    for m in prange(N):
        adj = np.zeros(n, np.int64); radj = np.zeros(n, np.int64)
        dp = np.zeros((1 << n, n), np.int64)
        adj_from_mask(n, m, pi, pj, adj, radj)
        H[m] = hp(n, adj, dp)
        S[m] = is_strong(n, adj, radj)
    return H, S


@njit(cache=True)
def fgQ(n, adj, f, g, Q):
    """f[S,x] = paths covering S ending at x; g[S,y] = paths covering S starting at y;
    Q[x,y] = ordered pairs (P1 ending at x, P2 starting at y), both nonempty, disjoint, covering V."""
    full = (1 << n) - 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S >> v) & 1:
                f[S, v] = 0; g[S, v] = 0
                continue
            R = S ^ (1 << v)
            if R == 0:
                f[S, v] = 1; g[S, v] = 1
                continue
            t = 0; t2 = 0
            for w in range(n):
                if (R >> w) & 1:
                    if (adj[w] >> v) & 1:
                        t += f[R, w]
                    else:
                        t2 += g[R, w]      # v -> w
            f[S, v] = t; g[S, v] = t2
    for x in range(n):
        for y in range(n):
            Q[x, y] = 0
    for S in range(1, full):
        R = full ^ S
        for x in range(n):
            if (S >> x) & 1:
                fx = f[S, x]
                if fx:
                    for y in range(n):
                        if (R >> y) & 1:
                            Q[x, y] += fx * g[R, y]
    h = 0
    for x in range(n):
        h += f[full, x]
    return h


def brute_H(n, adjm):
    """pure python brute force over permutations; adjm[i][j] = 1 iff i -> j."""
    import itertools
    c = 0
    for p in itertools.permutations(range(n)):
        ok = True
        for i in range(n - 1):
            if not adjm[p[i]][p[i + 1]]:
                ok = False; break
        if ok:
            c += 1
    return c


def adjm_from_mask(n, m):
    A = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(pairs_of(n)):
        if (m >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    return A
