#!/usr/bin/env python3
"""Lane thinorbit (S93).  Fast (numba, parallel) version of gen_blocks.py for one block order K.

Same definition and same method: for every marked block (B, w, u2) of order K the matrix N is solved from the three
base ported tournaments (singleton; arc with port at the source; 3-cycle), states computed by direct Hamiltonian
path DP on the glued tournament.  No further bases are checked here (gen_blocks.py checks the law on 38 further
bases for K <= 5; the law itself is THM-5070 (1), proved).  A spot check on NCHK random marked blocks against a
fourth base (the transitive triple with port at its middle vertex) is included.
Usage: gen_blocks_fast.py K OUT      (writes the distinct matrices of order exactly K, with representatives)
"""
import sys

import numpy as np
from numba import njit, prange


@njit(cache=True)
def state_of(n, adj, port, dp, f, Hs):
    full = (1 << n) - 1
    for mask in range(1 << n):
        Hs[mask] = 0
        for v in range(n):
            dp[mask, v] = 0
            f[mask, v] = 0
    for v in range(n):
        dp[1 << v, v] = 1
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


@njit(cache=True, parallel=True)
def all_blocks(K, base_n, base_adj, base_port, Sinv, out):
    npairs = K * (K - 1) // 2
    ntour = 1 << npairs
    nb = base_n.shape[0]
    for code in prange(ntour):
        adjB = np.zeros((K, K), np.int64)
        t = 0
        for i in range(K):
            for j in range(i + 1, K):
                if (code >> t) & 1:
                    adjB[i, j] = 1
                else:
                    adjB[j, i] = 1
                t += 1
        nmax = K + 3
        dp = np.zeros((1 << nmax, nmax), np.int64)
        f = np.zeros((1 << nmax, nmax), np.int64)
        Hs = np.zeros(1 << nmax, np.int64)
        adj = np.zeros((nmax, nmax), np.int64)
        R = np.zeros((3, nb), np.int64)
        for w in range(K):
            for up in range(K):
                for b in range(nb):
                    nT = base_n[b]
                    n = nT + K
                    for i in range(n):
                        for j in range(n):
                            adj[i, j] = 0
                    for i in range(nT):
                        for j in range(nT):
                            adj[i, j] = base_adj[b, i, j]
                    for i in range(K):
                        for j in range(K):
                            adj[nT + i, nT + j] = adjB[i, j]
                    for i in range(nT):
                        for j in range(K):
                            adj[i, nT + j] = 1
                    u = base_port[b]
                    adj[u, nT + w] = 0
                    adj[nT + w, u] = 1
                    h, E, O = state_of(n, adj, nT + up, dp, f, Hs)
                    R[0, b] = h
                    R[1, b] = E
                    R[2, b] = O
                idx = (code * K + w) * K + up
                for i in range(3):
                    for j in range(3):
                        s = 0
                        for b in range(3):
                            s += R[i, b] * Sinv[b, j]
                        out[idx, i * 3 + j] = s
                # fourth base as a check: predicted state must equal the computed one
                if nb > 3:
                    for i in range(3):
                        s = 0
                        for j in range(3):
                            s += out[idx, i * 3 + j] * base_state3[j]
                        if s != R[i, 3]:
                            out[idx, 0] = -1


base_state3 = np.zeros(3, np.int64)


def main():
    K = int(sys.argv[1])
    outp = sys.argv[2]
    base_n = np.array([1, 2, 3, 3], np.int64)
    base_adj = np.zeros((4, 3, 3), np.int64)
    base_adj[1, 0, 1] = 1
    base_adj[2, 0, 1] = 1
    base_adj[2, 1, 2] = 1
    base_adj[2, 2, 0] = 1
    base_adj[3, 0, 1] = 1
    base_adj[3, 0, 2] = 1
    base_adj[3, 1, 2] = 1
    base_port = np.array([0, 0, 0, 1], np.int64)
    dp = np.zeros((1 << 3, 3), np.int64)
    f = np.zeros((1 << 3, 3), np.int64)
    Hs = np.zeros(1 << 3, np.int64)
    sts = [state_of(int(base_n[b]), base_adj[b][:base_n[b], :base_n[b]].copy(), int(base_port[b]), dp, f, Hs) for b in range(4)]
    print("base states:", sts, file=sys.stderr)
    import sympy as sp
    S = sp.Matrix(3, 3, lambda i, j: int(sts[j][i]))
    assert S.det() in (1, -1)
    Sinv = np.array(S.inv().tolist(), dtype=np.int64)
    global base_state3
    base_state3[:] = [int(x) for x in sts[3]]
    npairs = K * (K - 1) // 2
    out = np.zeros(((1 << npairs) * K * K, 9), np.int64)
    all_blocks(K, base_n, base_adj, base_port, Sinv, out)
    bad = int((out[:, 0] == -1).sum())
    uniq, first = np.unique(out, axis=0, return_index=True)
    neg = int((uniq < 0).any(axis=1).sum())
    with open(outp, "w") as fh:
        for row in uniq:
            fh.write(" ".join(str(int(x)) for x in row) + "\n")
    with open(outp + ".reps", "w") as fh:
        for row, idx in zip(uniq, first):
            code = idx // (K * K)
            w = (idx // K) % K
            up = idx % K
            adjB = [[0] * K for _ in range(K)]
            t = 0
            for i in range(K):
                for j in range(i + 1, K):
                    if (code >> t) & 1:
                        adjB[i][j] = 1
                    else:
                        adjB[j][i] = 1
                    t += 1
            rows = "".join(str(adjB[i][j]) for i in range(K) for j in range(K))
            fh.write(" ".join(str(int(x)) for x in row) + f" | {K} {rows} {w} {up}\n")
    print(f"order {K}: {out.shape[0]} marked blocks, {len(uniq)} distinct matrices, {neg} with a negative entry, "
          f"fourth-base check failures: {bad}")


if __name__ == "__main__":
    main()
