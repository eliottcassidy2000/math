#!/usr/bin/env python3
"""S91 (THM-5129 (E)): the head correction is an alternating sum over qualifying paths.

For a rooted tournament (B, y) let U_even, U_odd be the sums of H(B - V(R)) over the directed paths R that start at
y and have an even, odd number of vertices (R = (y) is odd), so U = U_even + U_odd is the slope U_B(y).  With
E = number of Hamiltonian paths of B ending at y:

    E = U_odd - U_even,        H(B + z) = D + 2 U_even = D - E + U        (z -> y, every other vertex -> z),
    H(B => TT_t => z) = D + 2^t U    (t >= 1).

So the bare apex closes only the qualifying paths with an even number of vertices (odd cycles z, R), and each sink
of a tail removes the parity condition.  (The first identity is Session261's marked odd-cycle formula at a
one-vertex out-neighbourhood combined with the pendant law of THM-4865 / THM-5126 (A).)

Checks (FINITE-EXACT, every root of every isomorphism class of order <= 8, 58,653 rooted tournaments):
  * E = U_odd - U_even;
  * H(B + z) computed by its own subset DP equals D + 2 U_even, for every class of order <= 6 and every root;
  * the four attained numbers that are not two-step values, 11, 13, 47, 67, are head values D - E + U: the least
    order of a base is printed with (D, E, U).
Exit status 1 on failure.
"""
import os
import sys

import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))


@njit(cache=True)
def tables(n, code, f, g, hS, out):
    full = (1 << n) - 1
    for i in range(n):
        out[i] = 0
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
            f[S, v] = 0
            g[S, v] = 0
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


@njit(cache=True)
def popcount(x):
    c = 0
    while x:
        c += 1
        x &= x - 1
    return c


@njit(cache=True)
def ham(n, out):
    """Hamiltonian paths of the tournament with out-masks out[0..n-1] (own DP, by last vertex)."""
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for S in range(1, full + 1):
        for v in range(n):
            c = dp[S, v]
            if c:
                m = out[v] & ~S & full
                while m:
                    w = 0
                    while ((m >> w) & 1) == 0:
                        w += 1
                    m &= m - 1
                    dp[S | (1 << w), w] += c
    tot = 0
    for v in range(n):
        tot += dp[full, v]
    return tot


@njit(cache=True)
def run(n, reps, rows, direct):
    """rows[k] = (D, E, U_even, U_odd, direct H(B + z) or -1) for every class and root."""
    full = (1 << n) - 1
    f = np.zeros((1 << n, n), np.int64)
    g = np.zeros((1 << n, n), np.int64)
    hS = np.zeros(1 << n, np.int64)
    out = np.zeros(n, np.int64)
    o2 = np.zeros(n + 1, np.int64)
    k = 0
    for r in range(reps.shape[0]):
        tables(n, reps[r], f, g, hS, out)
        for y in range(n):
            ue = 0
            uo = 0
            for S in range(1, full + 1):
                if (S >> y) & 1:
                    term = g[S, y] * hS[full ^ S]
                    if popcount(S) & 1:
                        uo += term
                    else:
                        ue += term
            rows[k, 0] = hS[full]
            rows[k, 1] = f[full, y]
            rows[k, 2] = ue
            rows[k, 3] = uo
            rows[k, 4] = -1
            if direct:
                for v in range(n):
                    o2[v] = out[v]
                    if v != y:
                        o2[v] |= (1 << n)
                o2[n] = 1 << y
                rows[k, 4] = ham(n + 1, o2)
            k += 1
    return k


def main():
    z = np.load(os.path.join(HERE, "reps_le8.npz"))
    ok = True
    total = 0
    heads = {}
    for n in range(1, 9):
        reps = z["n%d" % n].astype(np.int64)
        rows = np.zeros((len(reps) * n, 5), np.int64)
        k = run(n, reps, rows, n <= 6)
        rows = rows[:k]
        total += k
        D, E, ue, uo, hd = (rows[:, i] for i in range(5))
        a = bool(np.all(E == uo - ue))
        b = bool(np.all((hd < 0) | (hd == D + 2 * ue)))
        ok &= a and b
        print(f"order {n}: {k:6d} rooted tournaments; E = U_odd - U_even: {a}; "
              f"H(B + z) = D + 2 U_even (direct count{'' if n <= 6 else ' not run'}): {b}")
        head = D - E + ue + uo
        for N in (11, 13, 47, 67):
            if N not in heads:
                idx = np.nonzero(head == N)[0]
                if len(idx):
                    i = idx[0]
                    heads[N] = (n, int(D[i]), int(E[i]), int(ue[i] + uo[i]))
    print(f"total rooted tournaments checked: {total}")
    for N in (11, 13, 47, 67):
        if N in heads:
            n, d, e, u = heads[N]
            print(f"  {N} = D - E + U = {d} - {e} + {u}: head value of a base of order {n} (tournament of order {n + 1})")
        else:
            print(f"  {N}: no head representation with a base of order <= 8")
            ok = False
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
