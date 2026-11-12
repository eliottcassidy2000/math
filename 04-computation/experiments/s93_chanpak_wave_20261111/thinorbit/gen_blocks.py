#!/usr/bin/env python3
"""Lane thinorbit (S93).  Independent recomputation of the base-2 port-chain block matrices of THM-5070 / HYP-9640.

Definitions (THM-5070).  Ported tournament (T, u), state s(T, u) = (h, E, O):
    h = H(T) (number of Hamiltonian paths),
    E = sum over directed paths P of T starting at u with |V(P)| even of H(T - V(P)),
    O = the same with |V(P)| odd (the path (u) included);  H(empty) = 1.
Gluing a marked block (B, w, u2): T' = T => B (every vertex of T beats every vertex of B) except w -> u; new port u2.
Claim: s(T', u2) = N(B, w, u2) s(T, u) with N a nonnegative integer 3 x 3 matrix depending only on the block.

Method here (independent of chain_blocks.py and of OCF): Hamiltonian paths are counted directly by subset DP.
For every marked block the matrix is SOLVED from three base ported tournaments with linearly independent states
(det = -1, so the solution is integral iff the law holds) and then CHECKED on further bases.

Usage: gen_blocks.py KMAX OUT          (writes OUT and OUT.reps)
"""
import itertools
import sys
from fractions import Fraction

import numpy as np
from numba import njit


@njit(cache=True)
def state_of(n, adj, port):
    """(h, E, O) of the ported tournament (adj, port) on n vertices; adj[i, j] = 1 iff i -> j."""
    full = (1 << n) - 1
    # dp[mask, v]: Hamiltonian paths of the subtournament on mask ending at v
    dp = np.zeros((1 << n, n), np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    Hs = np.zeros(1 << n, np.int64)
    Hs[0] = 1
    for mask in range(1, 1 << n):
        tot = 0
        for v in range(n):
            c = dp[mask, v]
            if c == 0:
                continue
            tot += c
            for t in range(n):
                if adj[v, t] == 1 and ((mask >> t) & 1) == 0:
                    dp[mask | (1 << t), t] += c
        Hs[mask] = tot
    # f[mask, v]: paths starting at port, vertex set mask, ending at v
    f = np.zeros((1 << n, n), np.int64)
    f[1 << port, port] = 1
    E = 0
    O = 0
    for mask in range(1, 1 << n):
        if ((mask >> port) & 1) == 0:
            continue
        tot = 0
        pc = 0
        for v in range(n):
            if (mask >> v) & 1:
                pc += 1
            c = f[mask, v]
            if c == 0:
                continue
            tot += c
            for t in range(n):
                if adj[v, t] == 1 and ((mask >> t) & 1) == 0:
                    f[mask | (1 << t), t] += c
        if tot:
            if pc % 2 == 0:
                E += tot * Hs[full & ~mask]
            else:
                O += tot * Hs[full & ~mask]
    return Hs[full], E, O


def glue(nT, adjT, u, k, adjB, w):
    n = nT + k
    adj = np.zeros((n, n), np.int64)
    adj[:nT, :nT] = adjT
    adj[nT:, nT:] = adjB
    adj[:nT, nT:] = 1
    adj[u, nT + w] = 0
    adj[nT + w, u] = 1
    return n, adj


def tournaments(k):
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    for bits in range(1 << len(pairs)):
        adj = np.zeros((k, k), np.int64)
        for t, (i, j) in enumerate(pairs):
            if bits >> t & 1:
                adj[i, j] = 1
            else:
                adj[j, i] = 1
        yield adj


def solve3(S, R):
    """rows r with r . S[i] = R[i] for the three base states S (3 x 3 integer, det +-1)."""
    S = [[Fraction(int(x)) for x in row] for row in S]
    # solve S^T-system by Cramer: N S_i = R_i  <=>  N = Rmat * Smat^{-1}, Smat columns = states
    import sympy
    Sm = sympy.Matrix(3, 3, lambda i, j: int(S[j][i]))
    Rm = sympy.Matrix(3, 3, lambda i, j: int(R[j][i]))
    N = Rm * Sm.inv()
    return N


def main():
    kmax = int(sys.argv[1])
    out = sys.argv[2]
    # base ported tournaments: singleton; arc with port at the source; 3-cycle; plus checks
    bases = []
    a1 = np.zeros((1, 1), np.int64)
    bases.append((1, a1, 0))
    a2 = np.array([[0, 1], [0, 0]], np.int64)
    bases.append((2, a2, 0))
    a3 = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], np.int64)
    bases.append((3, a3, 0))
    # checks: every ported tournament of order 3 and a spread of order 4, 5
    rng = np.random.default_rng(20261111)
    checks = []
    for adj in tournaments(3):
        for u in range(3):
            checks.append((3, adj.copy(), u))
    all4 = list(tournaments(4))
    for idx in rng.choice(len(all4), 10, replace=False):
        checks.append((4, all4[idx].copy(), int(rng.integers(0, 4))))
    for _ in range(4):
        adj = np.zeros((5, 5), np.int64)
        for i in range(5):
            for j in range(i + 1, 5):
                if rng.integers(0, 2):
                    adj[i, j] = 1
                else:
                    adj[j, i] = 1
        checks.append((5, adj, int(rng.integers(0, 5))))
    base_states = [state_of(n, a, u) for (n, a, u) in bases]
    check_states = [state_of(n, a, u) for (n, a, u) in checks]
    print("base states:", base_states, file=sys.stderr)
    import sympy
    Sm = sympy.Matrix(3, 3, lambda i, j: int(base_states[j][i]))
    print("det of base states:", Sm.det(), file=sys.stderr)
    Sinv = Sm.inv()
    mats = {}
    nblocks = 0
    nchecks = 0
    for k in range(1, kmax + 1):
        for adjB in tournaments(k):
            for w in range(k):
                for up in range(k):
                    nblocks += 1
                    R = []
                    for (n, a, u) in bases:
                        nn, adj = glue(n, a, u, k, adjB, w)
                        R.append(state_of(nn, adj, n + up))
                    Rm = sympy.Matrix(3, 3, lambda i, j: int(R[j][i]))
                    N = Rm * Sinv
                    M = tuple(int(N[i, j]) for i in range(3) for j in range(3))
                    assert all(N[i, j] == int(N[i, j]) for i in range(3) for j in range(3))
                    # checks
                    for (n, a, u), st in zip(checks, check_states):
                        nn, adj = glue(n, a, u, k, adjB, w)
                        got = state_of(nn, adj, n + up)
                        pred = tuple(M[3 * i] * st[0] + M[3 * i + 1] * st[1] + M[3 * i + 2] * st[2] for i in range(3))
                        nchecks += 1
                        if tuple(int(x) for x in got) != pred:
                            print("TRANSFER FAILS", k, adjB.tolist(), w, up, (n, a.tolist(), u), got, pred)
                            sys.exit(1)
                    if M not in mats:
                        mats[M] = (k, "".join(str(int(adjB[i, j])) for i in range(k) for j in range(k)), w, up)
        print(f"order <= {k}: {nblocks} marked blocks, {len(mats)} distinct matrices, {nchecks} transfer checks passed",
              file=sys.stderr)
    with open(out, "w") as f:
        for M in sorted(mats):
            f.write(" ".join(map(str, M)) + "\n")
    with open(out + ".reps", "w") as f:
        for M in sorted(mats):
            k, rows, w, up = mats[M]
            f.write(" ".join(map(str, M)) + f" | {k} {rows} {w} {up}\n")
    neg = sum(1 for M in mats if min(M) < 0)
    print(f"{nblocks} marked blocks of order <= {kmax}; {len(mats)} distinct matrices; {neg} with a negative entry; "
          f"{nchecks} transfer checks on {len(checks)} further bases, 0 failures")


if __name__ == "__main__":
    main()
