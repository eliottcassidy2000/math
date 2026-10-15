#!/usr/bin/env python3
"""The pointwise complement pairing behind the stochastic symmetry (SS) of HYP-9280 is false: a random test.

Fix a strong tournament T and a vertex w with in-neighbours I and out-neighbours O, and write V' = V - w.
A Hamiltonian path of T is A w B, where A is a path on a set X subset of V' ending in I (or empty) and B is a path on
X^c = V' - X starting in O (or empty).  The weight of X is
    wt(X) = #{Hamiltonian paths of T[X] ending in I} * #{Hamiltonian paths of T[X^c] starting in O}
(an empty factor counts 1), so Pr[the set before w is X] = wt(X)/H(T) and N = |X cap I|.
The pointwise strengthening of (SS) would pair every X having fewer in-neighbours before w than after
(|X cap I| < |X^c cap I|) with its complement and ask wt(X^c) >= wt(X).  This script tests that inequality on
random strong tournaments with n in {5, 6, 7} (seed 5) and reports the number of tested triples (T, w, X) and of
violations.  (SS) itself held for all strong n <= 8 and fails at n = 9; this test shows (SS) cannot be proved
set-by-set.
Session opus-2026-10-14-S36.  Shipped after audit item E9.
"""
import itertools
import random


def count_paths(adj, verts, start_set=None, end_set=None):
    verts = list(verts)
    if not verts:
        return 1
    idx = {v: i for i, v in enumerate(verts)}
    m = len(verts)
    F = {}
    for v in verts:
        if start_set is None or v in start_set:
            F[(1 << idx[v], v)] = 1
    for S in range(1, 1 << m):
        for v in verts:
            c = F.get((S, v))
            if not c:
                continue
            for u in verts:
                if not S >> idx[u] & 1 and adj[v] >> u & 1:
                    key = (S | 1 << idx[u], u)
                    F[key] = F.get(key, 0) + c
    full = (1 << m) - 1
    return sum(F.get((full, v), 0) for v in verts if end_set is None or v in end_set)


def strongly_connected(n, adj):
    def reach(start, forward):
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for u in range(n):
                if u not in seen and ((adj[v] >> u & 1) if forward else (adj[u] >> v & 1)):
                    seen.add(u)
                    stack.append(u)
        return seen
    return len(reach(0, True)) == n and len(reach(0, False)) == n


if __name__ == '__main__':
    random.seed(5)
    viol = 0
    tested = 0
    ex = None
    ntourn = 0
    for trial in range(1500):
        n = random.choice([5, 6, 7])
        adj = [0] * n
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < 0.5:
                    adj[i] |= 1 << j
                else:
                    adj[j] |= 1 << i
        if not strongly_connected(n, adj):
            continue
        ntourn += 1
        w = random.randrange(n)
        I = {u for u in range(n) if u != w and adj[u] >> w & 1}
        O = {v for v in range(n) if v != w and adj[w] >> v & 1}
        rest = [v for v in range(n) if v != w]
        for r in range(len(rest) + 1):
            for X in itertools.combinations(rest, r):
                X = set(X)
                Xc = set(rest) - X
                if len(X & I) >= len(Xc & I):
                    continue
                wt = (count_paths(adj, X, end_set=I) if X else 1) * (count_paths(adj, Xc, start_set=O) if Xc else 1)
                wtc = (count_paths(adj, Xc, end_set=I) if Xc else 1) * (count_paths(adj, X, start_set=O) if X else 1)
                tested += 1
                if wtc < wt:
                    viol += 1
                    if ex is None:
                        ex = (n, adj, w, sorted(X), wt, wtc)
    print(f'random strong tournaments tested: {ntourn} (one random vertex w each)')
    print(f'triples (T, w, X) with |X cap I| < |X^c cap I|: {tested}; violations of wt(X^c) >= wt(X): {viol}')
    print(f'first violation (n, adj, w, X, wt(X), wt(X^c)): {ex}')
