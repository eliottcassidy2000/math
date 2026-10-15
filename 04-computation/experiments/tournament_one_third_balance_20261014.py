#!/usr/bin/env python3
"""A tournament analogue of the 1/3-2/3 conjecture (Kislitsyn-Fredman-Linial), tested exhaustively.

For a tournament T, take a uniformly random Hamiltonian path and let p(x,y) = Pr[x comes before y].  Define the balance
delta(T) = max_{x != y} min(p(x,y), 1 - p(x,y)).  A transitive tournament (the analogue of a chain) has a unique
Hamiltonian path, so delta = 0.

HYPOTHESIS (tournament 1/3-2/3): every non-transitive tournament has delta(T) >= 1/3.  Equality holds for the cyclic
triangle C_3: its three Hamiltonian paths give p = 2/3 for every arc.

Counting: with F[S][v] = # Hamiltonian paths of T[S] ending at v and B[S][v] = # starting at v,
#{paths with x before y} = sum_{S: x in S, y notin S} sum_{u in S, u -> y} F[S][u] B[V \\ S][y].
Exhaustive over all 2^C(n,2) labelled tournaments for n <= 7 (numba, parallel).  Also the Paley arc balance
rho(q) = Pr[0 before 1] (1 a nonzero square) for q = 3, 7, 11, 19.
Session opus-2026-10-14-S35.
"""
import sys, time
import numpy as np
import numba as nb
from fractions import Fraction

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7


@nb.njit(cache=True)
def fb_tables(n, adj):
    """adj[u] = bitmask of v with u -> v.  F[S][v], B[S][v] as int64"""
    N = 1 << n
    F = np.zeros((N, n), dtype=np.int64)
    B = np.zeros((N, n), dtype=np.int64)
    for v in range(n):
        F[1 << v, v] = 1
        B[1 << v, v] = 1
    for S in range(1, N):
        for v in range(n):
            if not (S >> v) & 1:
                continue
            c = F[S, v]
            if c:
                # extend forward: v -> w
                m = adj[v] & ~S
                while m:
                    w = 0
                    low = m & (-m)
                    while (low >> w) != 1:
                        w += 1
                    F[S | low, w] += c
                    m &= m - 1
            c = B[S, v]
            if c:
                # extend backward: w -> v (w precedes v)
                for w in range(n):
                    if (S >> w) & 1:
                        continue
                    if (adj[w] >> v) & 1:
                        B[S | (1 << w), w] += c
    return F, B


@nb.njit(cache=True)
def balance(n, adj):
    F, B = fb_tables(n, adj)
    full = (1 << n) - 1
    H = 0
    for v in range(n):
        H += F[full, v]
    best_num = 0  # maximise min(c, H - c) over pairs
    for x in range(n):
        for y in range(n):
            if x == y:
                continue
            cnt = 0
            for S in range(1, full + 1):
                if not (S >> x) & 1 or (S >> y) & 1:
                    continue
                comp = full ^ S
                b = B[comp, y]
                if b == 0:
                    continue
                for u in range(n):
                    if (S >> u) & 1 and (adj[u] >> y) & 1:
                        cnt += F[S, u] * b
            m = cnt if cnt < H - cnt else H - cnt
            if m > best_num:
                best_num = m
    return H, best_num


@nb.njit(parallel=True, cache=True)
def sweep(n, pairs_u, pairs_v, nchunk):
    E = pairs_u.shape[0]
    total = 1 << E
    # per chunk: minimum ratio best_num/H over non-transitive tournaments, as (num, den, code)
    out_num = np.zeros(nchunk, dtype=np.int64)
    out_den = np.ones(nchunk, dtype=np.int64)
    out_code = np.zeros(nchunk, dtype=np.int64)
    cnt_viol = np.zeros(nchunk, dtype=np.int64)
    chunk = (total + nchunk - 1) // nchunk
    for c in nb.prange(nchunk):
        bn, bd, bc = 1, 1, -1
        adj = np.zeros(n, dtype=np.int64)
        for code in range(c * chunk, min(total, (c + 1) * chunk)):
            for i in range(n):
                adj[i] = 0
            for e in range(E):
                u = pairs_u[e]; v = pairs_v[e]
                if (code >> e) & 1:
                    adj[u] |= 1 << v
                else:
                    adj[v] |= 1 << u
            H, m = balance(n, adj)
            if H == 1:
                continue  # transitive (unique Hamiltonian path)
            # compare m/H < bn/bd
            if m * bd < bn * H:
                bn, bd, bc = m, H, code
            if 3 * m < H:
                cnt_viol[c] += 1
        out_num[c] = bn; out_den[c] = bd; out_code[c] = bc
    return out_num, out_den, out_code, cnt_viol


def paley_rho(q):
    sq = {(x * x) % q for x in range(1, q)}
    adj = np.zeros(q, dtype=np.int64)
    for u in range(q):
        for v in range(q):
            if u != v and (v - u) % q in sq:
                adj[u] |= 1 << v
    F, B = fb_tables(q, adj)
    full = (1 << q) - 1
    H = int(F[full].sum())
    x, y = 0, 1
    cnt = 0
    for S in range(1, full + 1):
        if not (S >> x) & 1 or (S >> y) & 1:
            continue
        b = int(B[full ^ S, y])
        if b == 0:
            continue
        for u in range(q):
            if (S >> u) & 1 and (adj[u] >> y) & 1:
                cnt += int(F[S, u]) * b
    return H, cnt


def main():
    for n in range(3, NMAX + 1):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        pu = np.array([p[0] for p in pairs], dtype=np.int64)
        pv = np.array([p[1] for p in pairs], dtype=np.int64)
        t0 = time.time()
        num, den, code, viol = sweep(n, pu, pv, 256)
        j = min(range(len(num)), key=lambda k: Fraction(int(num[k]), int(den[k])))
        r = Fraction(int(num[j]), int(den[j]))
        print(f'n = {n}: min over non-transitive tournaments of delta(T) = {r} = {float(r):.4f}  '
              f'(attained at labelled code {int(code[j])});  tournaments with delta < 1/3: {int(viol.sum())};  {time.time() - t0:.1f} s', flush=True)
    for q in (3, 7, 11, 19):
        H, c = paley_rho(q)
        print(f'Paley P_{q}: H = {H}, #paths with 0 before 1 = {c}, rho = {Fraction(c, H)} = {c / H:.6f}')


if __name__ == '__main__':
    main()
