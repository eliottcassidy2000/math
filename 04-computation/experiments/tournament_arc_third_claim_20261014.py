#!/usr/bin/env python3
"""Claim A (a sufficient condition for the tournament 1/3-2/3 theorem, HYP-9277): for every arc x -> y of every
tournament, the tail precedes the head in at least one third of the Hamiltonian paths: 3 N(x before y) >= H.

PROVED companion (Claim B): in a non-transitive tournament some arc x -> y has p(x before y) <= 2/3 (otherwise the
relation p > 2/3 contains every arc, is transitive by p(x,z) >= p(x,y) + p(y,z) - 1, and T is transitive).
So Claim A for all tournaments implies part 1 of HYP-9277 (a balanced arc), and Claim A with strict inequality at
the B-arc would give part 2.
Exhaustive over labelled tournaments, n <= NMAX (numba, parallel): reports min over tournaments and arcs of
N(x before y)/H, and the number of violations.
Session opus-2026-10-14-S35.
"""
import sys, time
import numpy as np
import numba as nb
from fractions import Fraction

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7


@nb.njit(cache=True)
def fb_tables(n, adj):
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
                for w in range(n):
                    if not (S >> w) & 1 and (adj[v] >> w) & 1:
                        F[S | (1 << w), w] += c
            c = B[S, v]
            if c:
                for w in range(n):
                    if not (S >> w) & 1 and (adj[w] >> v) & 1:
                        B[S | (1 << w), w] += c
    return F, B


@nb.njit(cache=True)
def arc_min(n, adj):
    F, B = fb_tables(n, adj)
    full = (1 << n) - 1
    H = 0
    for v in range(n):
        H += F[full, v]
    best_num = H; best_den = H
    for x in range(n):
        for y in range(n):
            if x == y or not (adj[x] >> y) & 1:
                continue
            cnt = 0
            for S in range(1, full + 1):
                if not (S >> x) & 1 or (S >> y) & 1:
                    continue
                b = B[full ^ S, y]
                if b == 0:
                    continue
                for u in range(n):
                    if (S >> u) & 1 and (adj[u] >> y) & 1:
                        cnt += F[S, u] * b
            if cnt * best_den < best_num * H:
                best_num = cnt; best_den = H
    return best_num, best_den


@nb.njit(parallel=True, cache=True)
def sweep(n, pu, pv, nchunk):
    E = pu.shape[0]
    total = 1 << E
    on = np.ones(nchunk, dtype=np.int64); od = np.ones(nchunk, dtype=np.int64)
    oc = np.full(nchunk, -1, dtype=np.int64); viol = np.zeros(nchunk, dtype=np.int64)
    chunk = (total + nchunk - 1) // nchunk
    for c in nb.prange(nchunk):
        adj = np.zeros(n, dtype=np.int64)
        for code in range(c * chunk, min(total, (c + 1) * chunk)):
            for i in range(n):
                adj[i] = 0
            for e in range(E):
                if (code >> e) & 1:
                    adj[pu[e]] |= 1 << pv[e]
                else:
                    adj[pv[e]] |= 1 << pu[e]
            a, b = arc_min(n, adj)
            if 3 * a < b:
                viol[c] += 1
            if a * od[c] < on[c] * b:
                on[c] = a; od[c] = b; oc[c] = code
    return on, od, oc, viol


if __name__ == '__main__':
    for n in range(3, NMAX + 1):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        pu = np.array([p[0] for p in pairs], dtype=np.int64); pv = np.array([p[1] for p in pairs], dtype=np.int64)
        t0 = time.time()
        on, od, oc, viol = sweep(n, pu, pv, 512)
        k = min(range(len(on)), key=lambda i: Fraction(int(on[i]), int(od[i])))
        print(f'n = {n}: min over tournaments and arcs of p(tail before head) = {Fraction(int(on[k]), int(od[k]))} '
              f'= {on[k] / od[k]:.4f} (code {int(oc[k])});  tournaments with an arc below 1/3: {int(viol.sum())}  [{time.time() - t0:.1f} s]', flush=True)
