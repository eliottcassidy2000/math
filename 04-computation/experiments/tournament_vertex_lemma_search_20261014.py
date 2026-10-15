#!/usr/bin/env python3
"""Adversarial search against the vertex lemma and its averaged form (HYP-9280), strong tournaments n = 8..12.

Objectives (simulated annealing over arc flips, strong tournaments only):
  (V) minimise over vertices w of  max over in-neighbours u of (2 p(u before w) - 1)   (vertex lemma: > 0)
  (D) minimise over vertices w of  D_in(w) = sum over in-neighbours u of (2 p(u before w) - 1)   (averaged form: > 0)
Session opus-2026-10-14-S36.
"""
import sys, time, math, random, os
import importlib.util
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('mc', os.path.join(HERE, 'tournament_majority_cycle_20261014.py'))
mc = importlib.util.module_from_spec(spec)
_argv = sys.argv; sys.argv = ['x']; spec.loader.exec_module(mc); sys.argv = _argv
pair_counts, is_strong = mc.pair_counts, mc.is_strong


@nb.njit(cache=True)
def margins(n, adj):
    H, C = pair_counts(n, adj)
    vmin = 10.0; dmin = 1e9
    for w in range(n):
        best = -2.0; dsum = 0.0
        for u in range(n):
            if u != w and (adj[u] >> w) & 1:
                m = (2.0 * C[u, w] - H) / H
                dsum += m
                if m > best:
                    best = m
        if best < vmin:
            vmin = best
        if dsum < dmin:
            dmin = dsum
    return vmin, dmin, H


def anneal(n, obj, iters, seed):
    rng = random.Random(seed)
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            adj[i] |= 1 << j
    adj[n - 1] |= 1; adj[0] &= ~(1 << (n - 1))
    for _ in range(rng.randint(0, n)):
        i, j = sorted(rng.sample(range(n), 2))
        if adj[i] >> j & 1:
            adj[i] &= ~(1 << j); adj[j] |= 1 << i
    def score(a):
        arr = np.array(a, dtype=np.int64)
        if not is_strong(n, arr):
            return None
        v, d, H = margins(n, arr)
        return (v if obj == 'V' else d), v, d, H
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
        T = 0.05 * (1 - it / iters) + 0.002
        if sc[0] <= cur[0] or rng.random() < math.exp(-(sc[0] - cur[0]) / T):
            adj, cur = new, sc
            if cur[0] < best[0][0]:
                best = (cur, list(adj))
    return best


if __name__ == '__main__':
    NS = [int(a) for a in sys.argv[1:]] or [8, 9, 10, 11, 12]
    for n in NS:
        t0 = time.time()
        iters = 1500 if n <= 10 else 600
        for obj in ('V', 'D'):
            res = [anneal(n, obj, iters, s) for s in range(5 if n <= 10 else 3)]
            (val, v, d, H), adj = min(res, key=lambda r: r[0][0])
            what = 'min vertex-lemma margin' if obj == 'V' else 'min averaged D_in'
            print(f'n = {n}: {what}: vertex margin = {v:.6f}, D_in = {d:.6f}, H = {H}, adj = {adj}  [{time.time() - t0:.0f} s]', flush=True)
