#!/usr/bin/env python3
"""Verifier (lane digits), check V1c: the law s(T', m) = N(X, W, m) s(T, z) by direct count for EVERY labelled block
with |X| = k <= 5 and EVERY labelled marked tournament (T, z) of order <= min(5, 8 - k)  (the coverage the lane's
docstring states; its script tested only the first block of each distinct matrix).  N = formula matrix of V1."""
import time
import numpy as np
from numba import njit, prange
from v1_onemark_matrices import path_table, block_matrix, direct_state


@njit(cache=True)
def all_bases(nmax):
    """labelled marked tournaments of order <= nmax: rows (n, z, adj[0..4]) and states."""
    cnt = 0
    for n in range(1, nmax + 1):
        cnt += (1 << (n * (n - 1) // 2)) * n
    B = np.zeros((cnt, 7), np.int64)
    S = np.zeros((cnt, 4), np.int64)
    idx = 0
    s = np.zeros(4, np.int64)
    for n in range(1, nmax + 1):
        adj = np.zeros(n, np.int64)
        for bits in range(1 << (n * (n - 1) // 2)):
            for i in range(n):
                adj[i] = 0
            t = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if (bits >> t) & 1:
                        adj[i] |= 1 << j
                    else:
                        adj[j] |= 1 << i
                    t += 1
            for z in range(n):
                direct_state(n, adj, z, s)
                B[idx, 0] = n
                B[idx, 1] = z
                for i in range(n):
                    B[idx, 2 + i] = adj[i]
                for q in range(4):
                    S[idx, q] = s[q]
                idx += 1
    return B, S


@njit(parallel=True, cache=True)
def check(k, B, S):
    npairs = k * (k - 1) // 2
    ntour = 1 << npairs
    fails = np.zeros(ntour, np.int64)
    tests = np.zeros(ntour, np.int64)
    for bits in prange(ntour):
        adjX = np.zeros(k, np.int64)
        t = 0
        for i in range(k):
            for j in range(i + 1, k):
                if (bits >> t) & 1:
                    adjX[i] |= 1 << j
                else:
                    adjX[j] |= 1 << i
                t += 1
        P = path_table(k, adjX)
        N = np.zeros((4, 4), np.int64)
        s2 = np.zeros(4, np.int64)
        a = np.zeros(10, np.int64)
        for W in range(1 << k):
            for m in range(-1, k):
                block_matrix(k, P, W, m, N)
                for b in range(B.shape[0]):
                    n = B[b, 0]
                    z = B[b, 1]
                    NN = n + k
                    xmask = ((1 << k) - 1) << n
                    for i in range(n):
                        a[i] = B[b, 2 + i]
                        if i != z:
                            a[i] |= xmask
                        else:
                            a[i] |= (((1 << k) - 1) & ~W) << n
                    for i in range(k):
                        a[n + i] = adjX[i] << n
                        if (W >> i) & 1:
                            a[n + i] |= 1 << z
                    direct_state(NN, a[:NN], z if m < 0 else n + m, s2)
                    tests[bits] += 1
                    for r in range(4):
                        v = 0
                        for q in range(4):
                            v += N[r, q] * S[b, q]
                        if v != s2[r]:
                            fails[bits] += 1
                            break
    return tests.sum(), fails.sum()


t0 = time.time()
tot = 0
totf = 0
for k in range(1, 6):
    nb = min(5, 8 - k)
    B, S = all_bases(nb)
    tests, fails = check(k, B, S)
    tot += tests
    totf += fails
    print(f"|X| = {k}: every labelled block x all {B.shape[0]} labelled marked tournaments of order <= {nb}: {tests} tests, failures {fails}   [{time.time() - t0:.1f} s]", flush=True)
print(f"total {tot} tests, failures {totf}")
