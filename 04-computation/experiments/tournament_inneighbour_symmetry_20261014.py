#!/usr/bin/env python3
"""Stochastic symmetry of the number of in-neighbours preceding a vertex (a sufficient condition for HYP-9280).

For a vertex w of in-degree d in a tournament T, let N = number of in-neighbours of w that precede w in a uniform
random Hamiltonian path.  Conjecture (SS): N stochastically dominates d - N, i.e. Pr[N >= d - k] >= Pr[N <= k] for all k.
(SS) gives E[N] >= d/2, hence some in-neighbour u has p(u before w) >= 1/2, and > 1/2 by parity (H odd): the vertex
lemma.  The vertex lemma at every vertex makes the forward-arc digraph F(T) have minimum in-degree >= 1, hence a directed
cycle, hence a directed 3-cycle of the majority tournament, hence a balanced pair in (1/2, 2/3]: the tournament 1/3-2/3
conjecture for strong tournaments (and, by the component reduction, for all non-transitive tournaments).
PROVED cases: k = 0 for every d (move w from the front to just after the last in-neighbour: injective, valid), which
gives the vertex lemma for d <= 2.
The distribution of N is computed exactly from prefix sets X (vertices before w):
  weight(X) = #(Ham paths of T[X] ending in I) * #(Ham paths of T[V - w - X] starting in O)   (empty factors = 1).
Exhaustive over labelled strong tournaments n <= NMAX.
Session opus-2026-10-14-S36.
"""
import sys, time, os
import importlib.util
import numpy as np
import numba as nb

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('mc', os.path.join(HERE, 'tournament_majority_cycle_20261014.py'))
mc = importlib.util.module_from_spec(spec)
_argv = sys.argv; sys.argv = ['x']; spec.loader.exec_module(mc); sys.argv = _argv
is_strong, fb_tables = mc.is_strong, mc.fb_tables


@nb.njit(cache=True)
def n_distribution(n, adj, w, F, B):
    full = (1 << n) - 1
    rest = full ^ (1 << w)
    I = 0; O = 0
    for u in range(n):
        if u == w:
            continue
        if (adj[u] >> w) & 1:
            I |= 1 << u
        else:
            O |= 1 << u
    d = 0
    m = I
    while m:
        m &= m - 1; d += 1
    dist = np.zeros(n, dtype=np.int64)
    # iterate over subsets X of rest
    X = rest
    while True:
        Y = rest ^ X
        # prefix factor
        if X == 0:
            a = 1
        else:
            a = 0
            for u in range(n):
                if (X >> u) & 1 and (I >> u) & 1:
                    a += F[X, u]
        if a:
            if Y == 0:
                b = 1
            else:
                b = 0
                for v in range(n):
                    if (Y >> v) & 1 and (O >> v) & 1:
                        b += B[Y, v]
            if b:
                k = 0
                mm = X & I
                while mm:
                    mm &= mm - 1; k += 1
                dist[k] += a * b
        if X == 0:
            break
        X = (X - 1) & rest
    return dist, d


@nb.njit(parallel=True, cache=True)
def sweep(n, pu, pv, nchunk):
    E = pu.shape[0]; total = 1 << E
    viol = np.zeros(nchunk, dtype=np.int64)
    tested = np.zeros(nchunk, dtype=np.int64)
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
            if not is_strong(n, adj):
                continue
            F, B = fb_tables(n, adj)
            for w in range(n):
                dist, d = n_distribution(n, adj, w, F, B)
                tested[c] += 1
                le = 0
                for k in range(d + 1):
                    le += dist[k]
                    ge = 0
                    for j in range(d - k, d + 1):
                        ge += dist[j]
                    if ge < le:
                        viol[c] += 1
                        break
    return viol, tested


if __name__ == '__main__':
    for n in range(3, NMAX + 1):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        pu = np.array([p[0] for p in pairs], dtype=np.int64); pv = np.array([p[1] for p in pairs], dtype=np.int64)
        t0 = time.time()
        v, t = sweep(n, pu, pv, 512)
        print(f'n = {n}: (strong T, vertex) pairs {int(t.sum())}; violations of Pr[N >= d - k] >= Pr[N <= k] (some k): {int(v.sum())}  [{time.time() - t0:.1f} s]', flush=True)
