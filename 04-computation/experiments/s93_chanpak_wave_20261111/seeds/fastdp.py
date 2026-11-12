#!/usr/bin/env python3
"""S93 lane "seeds": the subset recursions of b_census.py with a larger bit table (orders up to 18).
dp_range: fe, gs, hS;  wu_tables: W_ef (symmetric) and U;  subset_sums_w: W(O), M(O) for every O."""
import numpy as np
from numba import njit

CTZ = np.zeros((1 << 18) + 1, np.int64)
for _i in range(19):
    CTZ[1 << _i] = _i


@njit(cache=True)
def popcount(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@njit(cache=True)
def decode(n, code, out):
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


@njit(cache=True)
def dp_range(n, out, fe, gs, hS, lo):
    """tables for the subsets S in [lo, 2^n) (lo = 1: all; lo = 2^(n-1): those containing the vertex n - 1)."""
    full = (1 << n) - 1
    if lo <= 1:
        hS[0] = 1
        for v in range(n):
            fe[0, v] = 0
            gs[0, v] = 0
        lo = 1
    for S in range(lo, full + 1):
        for v in range(n):
            fe[S, v] = 0
            gs[S, v] = 0
        if (S & (S - 1)) == 0:
            v = CTZ[S]
            fe[S, v] = 1
            gs[S, v] = 1
            hS[S] = 1
            continue
        t = 0
        rem = S
        while rem:
            lb = rem & (-rem)
            rem ^= lb
            v = CTZ[lb]
            R = S ^ lb
            ov = out[v]
            a = 0
            x = R & ~ov
            while x:
                l2 = x & (-x)
                x ^= l2
                a += fe[R, CTZ[l2]]
            b = 0
            x = R & ov
            while x:
                l2 = x & (-x)
                x ^= l2
                b += gs[R, CTZ[l2]]
            fe[S, v] = a
            gs[S, v] = b
            t += a
        hS[S] = t


@njit(cache=True)
def wu_tables(n, fe, gs, hS, Wm, U):
    """Wm[e, f] = W_ef = Q(e, f) + Q(f, e) (symmetric), U[o]."""
    full = (1 << n) - 1
    top = 1 << (n - 1)
    for e in range(n):
        U[e] = 0
        for f in range(n):
            Wm[e, f] = 0
    for S in range(0, full):
        C = full ^ S
        h = hS[S]
        x = C
        while x:
            lb = x & (-x)
            x ^= lb
            f = CTZ[lb]
            g = gs[C, f]
            U[f] += h * g
        if S & top:                    # each unordered partition {S, C} once (S contains the top vertex, C nonempty)
            x = C
            while x:
                lb = x & (-x)
                x ^= lb
                f = CTZ[lb]
                g = gs[C, f]
                a = fe[C, f]
                y = S
                while y:
                    l2 = y & (-y)
                    y ^= l2
                    e = CTZ[l2]
                    w = fe[S, e] * g + gs[S, e] * a
                    Wm[e, f] += w
                    Wm[f, e] += w


@njit(cache=True)
def subset_sums_w(n, Wm, U, Wsub, Msub):
    full = (1 << n) - 1
    Wsub[0] = 0
    Msub[0] = 0
    for O in range(1, full + 1):
        low = O & (-O)
        l = CTZ[low]
        R = O ^ low
        w = Wsub[R]
        x = R
        while x:
            l2 = x & (-x)
            x ^= l2
            w += Wm[l, CTZ[l2]]
        Wsub[O] = w
        Msub[O] = Msub[R] + U[l]


