#!/usr/bin/env python3
"""Adversarial search (simulated annealing over arc flips) for tournaments with small balance, n = 8..11.

Objectives:
  (A) min over arcs x -> y of p(x before y)       -- Claim A (sufficient for HYP-9277 part 1) predicts >= 1/3;
  (S) delta(T) over STRONG tournaments            -- HYP-9277 part 3 predicts >= 4/11; part 1 predicts >= 1/3.
p is over uniformly random Hamiltonian paths (forward/backward subset DP, numba).
Session opus-2026-10-14-S35.
"""
import sys, time, math, random
import numpy as np
import numba as nb
from fractions import Fraction


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
def pair_counts(n, adj):
    F, B = fb_tables(n, adj)
    full = (1 << n) - 1
    H = 0
    for v in range(n):
        H += F[full, v]
    C = np.zeros((n, n), dtype=np.int64)  # C[x, y] = # paths with x before y
    for x in range(n):
        for y in range(n):
            if x == y:
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
            C[x, y] = cnt
    return H, C


@nb.njit(cache=True)
def is_strong(n, adj):
    full = (1 << n) - 1
    seen = 1; fr = 1
    while fr:
        nx = 0
        for v in range(n):
            if (fr >> v) & 1:
                nx |= adj[v]
        nx &= ~seen; seen |= nx; fr = nx
    if seen != full:
        return False
    seen = 1; fr = 1
    while fr:
        nx = 0
        for v in range(n):
            if (fr >> v) & 1:
                for u in range(n):
                    if (adj[u] >> v) & 1:
                        nx |= 1 << u
        nx &= ~seen; seen |= nx; fr = nx
    return seen == full


def objectives(n, adj):
    H, C = pair_counts(n, np.array(adj, dtype=np.int64))
    arc = min(Fraction(int(C[x, y]), int(H)) for x in range(n) for y in range(n) if x != y and adj[x] >> y & 1)
    delta = max(Fraction(min(int(C[x, y]), int(H - C[x, y])), int(H)) for x in range(n) for y in range(x + 1, n))
    return H, arc, delta


def anneal(n, which, iters, seed):
    rng = random.Random(seed)
    # start near-transitive: transitive order plus a few random backward arcs
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            adj[i] |= 1 << j
    for _ in range(rng.randint(1, n)):
        i, j = sorted(rng.sample(range(n), 2))
        adj[i] &= ~(1 << j); adj[j] |= 1 << i
    def score(a):
        if which == 'S' and not is_strong(n, np.array(a, dtype=np.int64)):
            return None
        H, arc, delta = objectives(n, a)
        return float(arc if which == 'A' else delta), arc, delta, H
    cur = score(adj)
    while cur is None:
        i, j = rng.sample(range(n), 2)
        if adj[i] >> j & 1:
            adj[i] &= ~(1 << j); adj[j] |= 1 << i
        cur = score(adj)
    best = (cur, list(adj))
    T0 = 0.05
    for it in range(iters):
        i, j = rng.sample(range(n), 2)
        new = list(adj)
        if new[i] >> j & 1:
            new[i] &= ~(1 << j); new[j] |= 1 << i
        else:
            new[j] &= ~(1 << i); new[i] |= 1 << j
        sc = score(new)
        if sc is None:
            continue
        T = T0 * (1 - it / iters) + 1e-4
        if sc[0] <= cur[0] or rng.random() < math.exp(-(sc[0] - cur[0]) / T):
            adj, cur = new, sc
            if cur[0] < best[0][0]:
                best = (cur, list(adj))
    return best


if __name__ == '__main__':
    NS = [int(a) for a in sys.argv[1:]] or [8, 9, 10]
    for n in NS:
        t0 = time.time()
        for which in ('A', 'S'):
            res = [anneal(n, which, 1500 if n <= 9 else 600, s) for s in range(6 if n <= 9 else 3)]
            (val, arc, delta, H), adj = min(res, key=lambda r: r[0][0])
            label = 'min arc ratio (Claim A)' if which == 'A' else 'min delta over strong tournaments'
            print(f'n = {n}: {label}: arc = {arc} ({float(arc):.4f}), delta = {delta} ({float(delta):.4f}), H = {H}, adj = {adj}  [{time.time() - t0:.0f} s]', flush=True)
