#!/usr/bin/env python3
"""Independent check of A_10(q) (single run of 10 steps, 11 distinct points) by parallel DFS:
A_10 = q (q-1) * sum over injective (v_0 = 0, v_1 = 1, v_2, ..., v_10) of prod chi(v_(i+1) - v_i)   (L even: scaling).
Referee's symbolic result: A_10 = P(q) + 34 q (q-1) (a_q + q), P = 42q^6 - 1320q^5 + 15334q^4 - 78286q^3 + 161554q^2 - 97324q.
Usage: python audit_A10_dfs_parallel.py q
"""
import sys, time
import numpy as np
import numba as nb


@nb.njit(parallel=True, cache=True)
def a10_parallel(q, chi):
    L = 10
    res = np.zeros(q, dtype=np.int64)
    for v2 in nb.prange(2, q):
        used = np.zeros(q, dtype=np.int64)
        seq = np.zeros(L + 1, dtype=np.int64)
        nxt = np.zeros(L + 2, dtype=np.int64)
        sgn = np.zeros(L + 1, dtype=np.int64)
        seq[0] = 0
        seq[1] = 1
        seq[2] = v2
        used[0] = 1
        used[1] = 1
        used[v2] = 1
        sgn[1] = chi[1]
        sgn[2] = sgn[1] * chi[(v2 - 1) % q]
        depth = 3
        nxt[3] = 0
        tot = 0
        while depth >= 3:
            if depth == L + 1:
                tot += sgn[L]
                depth -= 1
                used[seq[depth]] = 0
                nxt[depth] = seq[depth] + 1
                continue
            x = nxt[depth]
            while x < q and used[x] == 1:
                x += 1
            if x >= q:
                depth -= 1
                if depth >= 3:
                    used[seq[depth]] = 0
                    nxt[depth] = seq[depth] + 1
                continue
            seq[depth] = x
            used[x] = 1
            sgn[depth] = sgn[depth - 1] * chi[(x - seq[depth - 1]) % q]
            depth += 1
            nxt[depth] = 0
        res[v2] = tot
    return res.sum() * q * (q - 1)


def chi_table(q):
    chi = np.full(q, -1, dtype=np.int64)
    chi[0] = 0
    for x in range(1, q):
        chi[(x * x) % q] = 1
    return chi


def eta_ap(N):
    P = [0] * (N + 1)
    P[0] = 1
    for n in range(1, N + 1):
        for step in (2 * n, 4 * n):
            if step > N:
                continue
            for _ in range(4):
                for i in range(N, step - 1, -1):
                    P[i] -= P[i - step]
    return [0] + P[:N]


def main():
    a = eta_ap(100)
    for q in [int(x) for x in sys.argv[1].split(',')]:
        t0 = time.time()
        v = int(a10_parallel(q, chi_table(q)))
        P = 42 * q ** 6 - 1320 * q ** 5 + 15334 * q ** 4 - 78286 * q ** 3 + 161554 * q ** 2 - 97324 * q
        pred = P + 34 * q * (q - 1) * (a[q] + q)
        print(f'q = {q}: DFS A_10 = {v}; symbolic P + 34q(q-1)(a_q+q) = {pred} (a_q = {a[q]}); agree {v == pred}; '
              f'pure polynomial P + 34q^2(q-1) would give {P + 34 * q * q * (q - 1)}; {time.time() - t0:.1f} s', flush=True)


if __name__ == '__main__':
    main()
