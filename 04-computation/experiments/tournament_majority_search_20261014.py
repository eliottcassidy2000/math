#!/usr/bin/env python3
"""Adversarial search against Conjecture M / F (strong tournaments, n = 8..13).

Objectives (simulated annealing over arc flips, restricted to strong tournaments):
  (I)  maximise the number of inverted arcs (arcs x -> y of T with p(x before y) < 1/2);
  (C)  minimise the number of directed 3-cycles of the majority tournament M(T)  (0 would refute Conjecture M);
  (F)  report whether the majority-forward arc digraph F(T) is strongly connected at the optimum.
Session opus-2026-10-14-S36.
"""
import sys, time, math, random
import numpy as np
import numba as nb
import importlib.util, os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('mc', os.path.join(HERE, 'tournament_majority_cycle_20261014.py'))
mc = importlib.util.module_from_spec(spec)
_argv = sys.argv; sys.argv = ['x']; spec.loader.exec_module(mc); sys.argv = _argv
pair_counts, is_strong = mc.pair_counts, mc.is_strong


@nb.njit(cache=True)
def stats(n, adj):
    H, C = pair_counts(n, adj)
    madj = np.zeros(n, dtype=np.int64)
    ninv = 0
    for x in range(n):
        for y in range(n):
            if x != y and 2 * C[x, y] > H:
                madj[x] |= 1 << y
                if not (adj[x] >> y) & 1:
                    ninv += 1
    tri = 0
    for x in range(n):
        for y in range(n):
            if (madj[x] >> y) & 1:
                for z in range(n):
                    if (madj[y] >> z) & 1 and (madj[z] >> x) & 1:
                        tri += 1
    tri //= 3
    fadj = np.zeros(n, dtype=np.int64)
    for x in range(n):
        fadj[x] = adj[x] & madj[x]
    return ninv, tri, is_strong(n, fadj), H


def anneal(n, objective, iters, seed):
    rng = random.Random(seed)
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            adj[i] |= 1 << j
    adj[n - 1] |= 1; adj[0] &= ~(1 << (n - 1))          # strong: transitive + reversed extreme arc
    for _ in range(rng.randint(0, n)):
        i, j = sorted(rng.sample(range(n), 2))
        if adj[i] >> j & 1:
            adj[i] &= ~(1 << j); adj[j] |= 1 << i
    def score(a):
        arr = np.array(a, dtype=np.int64)
        if not is_strong(n, arr):
            return None
        ninv, tri, fs, H = stats(n, arr)
        val = -ninv if objective == 'I' else tri
        return val, ninv, tri, fs, H
    cur = score(adj)
    while cur is None:
        i, j = rng.sample(range(n), 2)
        if adj[i] >> j & 1:
            adj[i] &= ~(1 << j); adj[j] |= 1 << i
        cur = score(adj)
    best = (cur, list(adj))
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
        T = 1.0 * (1 - it / iters) + 0.05
        if sc[0] <= cur[0] or rng.random() < math.exp(-(sc[0] - cur[0]) / T):
            adj, cur = new, sc
            if cur[0] < best[0][0]:
                best = (cur, list(adj))
    return best


if __name__ == '__main__':
    NS = [int(a) for a in sys.argv[1:]] or [8, 9, 10, 11, 12]
    for n in NS:
        t0 = time.time()
        iters = 1200 if n <= 10 else 500
        for obj in ('I', 'C'):
            res = [anneal(n, obj, iters, s) for s in range(5 if n <= 10 else 3)]
            (val, ninv, tri, fs, H), adj = min(res, key=lambda r: r[0][0])
            what = 'max #inverted arcs' if obj == 'I' else 'min #3-cycles of M(T)'
            print(f'n = {n}: {what}: inverted = {ninv}, M 3-cycles = {tri}, F strongly connected = {fs}, H = {H}, adj = {adj}  [{time.time() - t0:.0f} s]', flush=True)
