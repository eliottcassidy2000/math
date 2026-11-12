#!/usr/bin/env python3
"""S93 lane "seeds": common numba routines (not canon; every statement is re-derived in PROOFS.txt).

Conventions.  A tournament on n <= 16 vertices is an int64 array out[v] = bitmask of the out-neighbours of v.
Notation of THM-5130 / THM-5126 for a pair (B, O):
  fe[S, v] = number of Hamiltonian paths of B[S] ending at v;   gs[S, v] = ... starting at v;   hS[S] = H(B[S]), hS[0] = 1
  Q[e, f]  = number of ordered pairs (P1, P2) of disjoint NONEMPTY paths covering B, P1 ending at e, P2 starting at f
  U[o]     = number of ordered pairs (P1, P2) covering B, P2 starting at o, P1 possibly empty
  W(O) = sum over ordered pairs e != f of O of Q[e, f];  M(O) = sum over o in O of U[o];  r(O) = H(B) - W(O)
Closing law (KNOWN): H(B => TT_t => v, v -> O) = r(O) + 2^t M(O).
"""
import numpy as np
from numba import njit


@njit(cache=True)
def dp_tables(n, out, fe, gs, hS):
    full = (1 << n) - 1
    for S in range(full + 1):
        hS[S] = 0
        for v in range(n):
            fe[S, v] = 0
            gs[S, v] = 0
    hS[0] = 1
    for v in range(n):
        fe[1 << v, v] = 1
        gs[1 << v, v] = 1
    for S in range(1, full + 1):
        if S & (S - 1):
            for v in range(n):
                if (S >> v) & 1:
                    R = S ^ (1 << v)
                    a = 0
                    b = 0
                    ov = out[v]
                    for u in range(n):
                        if (R >> u) & 1:
                            if (ov >> u) & 1:
                                b += gs[R, u]          # v -> u : v first, then a path of R starting at u
                            else:
                                a += fe[R, u]          # u -> v : a path of R ending at u, then v
                    fe[S, v] = a
                    gs[S, v] = b
        t = 0
        for v in range(n):
            t += fe[S, v]
        hS[S] = t


@njit(cache=True)
def pair_mats(n, fe, gs, hS, Q, U):
    full = (1 << n) - 1
    for e in range(n):
        U[e] = 0
        for f in range(n):
            Q[e, f] = 0
    for S in range(0, full):
        C = full ^ S
        h = hS[S]
        for f in range(n):
            if (C >> f) & 1:
                g = gs[C, f]
                if g:
                    U[f] += h * g
                    if S:
                        for e in range(n):
                            if (S >> e) & 1:
                                Q[e, f] += fe[S, e] * g


@njit(cache=True)
def subset_sums(n, Q, U, Wsub, Msub):
    full = (1 << n) - 1
    Wsub[0] = 0
    Msub[0] = 0
    for O in range(1, full + 1):
        low = O & (-O)
        l = 0
        while (low >> l) != 1:
            l += 1
        R = O ^ low
        w = Wsub[R]
        for e in range(n):
            if (R >> e) & 1:
                w += Q[l, e] + Q[e, l]
        Wsub[O] = w
        Msub[O] = Msub[R] + U[l]


@njit(cache=True)
def hpaths(n, out):
    """H of a tournament on n <= 24 vertices by the subset recursion (ending vertex only)."""
    full = (1 << n) - 1
    dp = np.zeros((full + 1, n), np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for S in range(1, full + 1):
        for v in range(n):
            c = dp[S, v]
            if c:
                rest = out[v] & ~S & full
                while rest:
                    low = rest & (-rest)
                    u = 0
                    while (low >> u) != 1:
                        u += 1
                    dp[S | low, u] += c
                    rest ^= low
    t = 0
    for v in range(n):
        t += dp[full, v]
    return t


@njit(cache=True)
def closed_out(n, out, O, t, res):
    """res = out-masks of the closing B => TT_t => v, v -> O, B - O -> v; order n + t + 1; closer = last vertex."""
    N = n + t + 1
    tailmask = 0
    for j in range(t):
        tailmask |= 1 << (n + j)
    vbit = 1 << (n + t)
    for x in range(n):
        m = out[x] | tailmask
        if not ((O >> x) & 1):
            m |= vbit
        res[x] = m
    for j in range(t):
        m = vbit
        for k in range(j + 1, t):
            m |= 1 << (n + k)
        res[n + j] = m
    res[n + t] = O
    return N


def out_from_code(n, code):
    """upper-triangle code, pairs (i, j), i < j in lexicographic order, bit b = 1 iff i -> j (format of kps certs)."""
    out = np.zeros(n, np.int64)
    b = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (code >> b) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            b += 1
    return out


@njit(cache=True)
def out_from_colcode(n, code, out):
    """column-major code: pairs (i, j), i < j, ordered by j then i; bit = 1 iff i -> j."""
    for v in range(n):
        out[v] = 0
    b = 0
    for j in range(1, n):
        for i in range(j):
            if (code >> b) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            b += 1
