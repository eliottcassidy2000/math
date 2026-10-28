#!/usr/bin/env python3
"""P1 of PREREGISTRATION.md (pushed in efdf0b0a5f before this script was written).  kind-pasteur-2026-10-27-S8.

(a) bigons B_k, (b) alternating C_{2k} + K_1, (c) dense core K_4 + K_2 + 2K_1, (d) positive controls: random oriented
trees (end-edge map) and random functional orientations (root-swap map) at n = 8, 9, with explicit equivariance checks.
"""
import random, itertools, json, sys
from math import comb
import networkx as nx
from carrier_core import analyse, automorphisms


def bigon(k):
    # s = 0, t = 1, path P: 0 -> 2 -> ... -> 1, path Q: 0 -> ... -> 1, solo = 2k
    p = [0] + list(range(2, 2 + k - 1)) + [1]
    q = [0] + list(range(2 + k - 1, 2 + 2 * (k - 1))) + [1]
    arcs = [(p[i], p[i + 1]) for i in range(k)] + [(q[i], q[i + 1]) for i in range(k)]
    return 2 * k + 1, arcs


def alternating(k):
    m = 2 * k
    arcs = [((i, (i + 1) % m) if i % 2 == 0 else ((i + 1) % m, i)) for i in range(m)]
    return m + 1, arcs


def end_edge_check(n, arcs, A):
    """Trees: Pi -> pairs of E via end edges of the connecting path; check bijective and equivariant (ordered)."""
    Gu = nx.Graph(); Gu.add_nodes_from(range(n)); Gu.add_edges_from(arcs)
    E = [frozenset(e) for e in Gu.edges()]
    def phi(u, w):
        path = nx.shortest_path(Gu, u, w)
        return (frozenset(path[:2]), frozenset(path[-2:]))
    Pi = [(a, b) for a, b in itertools.combinations(range(n), 2) if not Gu.has_edge(a, b)]
    images = [phi(a, b) for a, b in Pi]
    bij = len(set(frozenset(x) for x in images)) == len(Pi) == comb(len(E), 2) and all(x[0] != x[1] for x in images)
    eq = True
    for g in A:
        for (a, b), (ea, eb) in zip(Pi, images):
            ga, gb = g[a], g[b]
            ge = (frozenset(g[v] for v in ea), frozenset(g[v] for v in eb))
            if phi(ga, gb) != ge:
                eq = False; break
    return bij, eq


def root_swap_check(n, arcs, A):
    """Functional orientation: unique sink r, out-degree 1 elsewhere. pairs(V - r) -> Pi:
    {v, f(v)} -> {r, v} if f(v) != r, identity on non-edges; orientation (v -> f(v)) <-> (r -> v)."""
    f = dict(arcs)
    r = [v for v in range(n) if v not in f][0]
    X = [v for v in range(n) if v != r]
    E = {frozenset(a) for a in arcs}
    def phi(x, y):  # ordered pair of X -> ordered pair of V (non-edge)
        if f.get(x) == y and y != r:
            return (r, x)
        if f.get(y) == x and x != r:
            return (y, r)   # (x -> y) with y -> x = f(y) an arc: maps to (y <- r)?? orientation convention below
        return (x, y)
    imgs = {}
    ok = True
    for x, y in itertools.permutations(X, 2):
        u, w = phi(x, y)
        if frozenset((u, w)) in E or u == w:
            ok = False
        imgs[(x, y)] = (u, w)
    # bijectivity on unordered pairs
    unordered = {frozenset(v) for v in imgs.values()}
    bij = ok and len(unordered) == comb(n - 1, 2)
    # reversal compatibility: phi(y, x) must be the reverse of phi(x, y)
    rev = all(imgs[(y, x)] == (imgs[(x, y)][1], imgs[(x, y)][0]) for x, y in itertools.permutations(X, 2))
    eq = all(imgs[(g[x], g[y])] == (g[imgs[(x, y)][0]], g[imgs[(x, y)][1]]) for g in A for x, y in itertools.permutations(X, 2))
    return bij, rev, eq


def random_tree_arcs(n, rng):
    if n == 2:
        return [(0, 1)]
    pr = [rng.randrange(n) for _ in range(n - 2)]
    T = nx.from_prufer_sequence(pr)
    return [(a, b) if rng.random() < 0.5 else (b, a) for a, b in T.edges()]


def random_functional_arcs(n, rng):
    while True:
        r = rng.randrange(n)
        f = {}
        for v in range(n):
            if v == r:
                continue
            f[v] = rng.choice([w for w in range(n) if w != v])
        if any(f.get(f[v]) == v for v in f):
            continue
        return [(v, f[v]) for v in f]


def main():
    rng = random.Random(20261027)
    out = {'bigon': [], 'alternating': [], 'dense': [], 'trees': [], 'functional': []}
    for k in range(2, 6):
        n, arcs = bigon(k)
        r = analyse(n, arcs, 1)
        r['predicted_marked'] = 2 ** (comb(2 * k, 2) - 1); r['k'] = k; r['n'] = n
        out['bigon'].append(r)
        print('bigon', k, r, flush=True)
    for k in range(2, 6):
        n, arcs = alternating(k)
        r = analyse(n, arcs, 1); r['k'] = k; r['n'] = n
        out['alternating'].append(r)
        print('alternating', k, r, flush=True)
    # dense core: transitive K4 on 0..3, K2 4->5, isolated 6, 7
    arcs = [(a, b) for a, b in itertools.combinations(range(4), 2)] + [(4, 5)]
    r = analyse(8, arcs, 1); out['dense'].append(r)
    print('dense K4+K2+2K1', r, flush=True)
    for n in (8, 9):
        tally = {'n': n, 'samples': 0, 'end_edge_bijective': 0, 'end_edge_equivariant': 0, 'maxwell': 0, 'linear': 0}
        for _ in range(150):
            arcs = random_tree_arcs(n, rng)
            A = automorphisms(n, arcs)
            b, e = end_edge_check(n, arcs, A)
            r = analyse(n, arcs, 1, want_weak=False)
            tally['samples'] += 1; tally['end_edge_bijective'] += b; tally['end_edge_equivariant'] += e
            tally['maxwell'] += bool(r['maxwell']); tally['linear'] += bool(r['linear'])
        out['trees'].append(tally); print('trees', tally, flush=True)
        tally = {'n': n, 'samples': 0, 'root_swap_bijective': 0, 'reversal_compatible': 0, 'equivariant': 0, 'maxwell': 0, 'linear': 0}
        for _ in range(150):
            arcs = random_functional_arcs(n, rng)
            A = automorphisms(n, arcs)
            b, rv, e = root_swap_check(n, arcs, A)
            r = analyse(n, arcs, 1, want_weak=False)
            tally['samples'] += 1; tally['root_swap_bijective'] += b; tally['reversal_compatible'] += rv
            tally['equivariant'] += e; tally['maxwell'] += bool(r['maxwell']); tally['linear'] += bool(r['linear'])
        out['functional'].append(tally); print('functional', tally, flush=True)
    json.dump(out, open('families_and_controls.json', 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
