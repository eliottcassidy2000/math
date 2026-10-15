#!/usr/bin/env python3
"""S36 referee audit, part D2: brute-force (all 9! orders, pure Python) confirmation of the five n = 9 isomorphism classes
that the exhaustive class sweep (audit_tournaments_classes_n9.py) found violating (SS), and identification of the
session's two counterexamples among them (networkx digraph isomorphism)."""
import itertools
from fractions import Fraction as Fr
import networkx as nx

mine = [[2, 16, 3, 7, 269, 31, 63, 127, 239], [2, 256, 3, 7, 271, 159, 63, 95, 237], [2, 32, 3, 7, 15, 285, 63, 127, 223],
        [2, 32, 3, 7, 47, 397, 63, 95, 223], [2, 128, 3, 7, 15, 159, 63, 349, 127]]
session = {'SS-1': [502, 508, 504, 241, 480, 448, 384, 256, 8], 'SS-2': [250, 508, 497, 500, 480, 448, 384, 256, 1]}


def paths(adj):
    n = len(adj)
    return [p for p in itertools.permutations(range(n)) if all((adj[p[i]] >> p[i + 1]) & 1 for i in range(n - 1))]


def dg(adj):
    G = nx.DiGraph()
    n = len(adj)
    G.add_nodes_from(range(n))
    G.add_edges_from((x, y) for x in range(n) for y in range(n) if (adj[x] >> y) & 1)
    return G


for k, adj in enumerate(mine):
    n = len(adj)
    P = paths(adj)
    H = len(P)
    pos = [{v: i for i, v in enumerate(p)} for p in P]
    for w in range(n):
        I = [u for u in range(n) if u != w and (adj[u] >> w) & 1]
        d = len(I)
        dist = [0] * (d + 1)
        for q in pos:
            dist[sum(1 for u in I if q[u] < q[w])] += 1
        viol = [j for j in range(d + 1) if sum(dist[d - j:]) < sum(dist[:j + 1])]
        if viol:
            best = max(Fr(sum(1 for q in pos if q[u] < q[w]), H) for u in I)
            EN = Fr(sum(j * dist[j] for j in range(d + 1)), H)
            print(f'class {k}: H = {H}, vertex {w} (in-degree {d}): N-distribution {dist}, SS fails at k = {viol}, '
                  f'E[N] = {EN} > d/2: {EN > Fr(d, 2)}, best in-neighbour p = {best}')
for name, adj in session.items():
    hits = [k for k, m in enumerate(mine) if nx.is_isomorphic(dg(adj), dg(m))]
    print(f'session {name} is isomorphic to audit class {hits}')
