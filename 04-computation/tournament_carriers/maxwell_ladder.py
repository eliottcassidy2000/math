#!/usr/bin/env python3
"""P2 of PREREGISTRATION.md (pushed in efdf0b0a5f before this script was written).  kind-pasteur-2026-10-27-S8.

Maxwell ladder: oriented carriers with j n - C(j+1,2) edges on n vertices (complement = C(n-j,2) pairs).  Test, over
ALL oriented carriers (orientation classes up to Aut(G)) for j = 2 (n = 4, 5, 6) and j = 3 (n = 5, 6), whether linear
j-jectivity coincides with the equivariant Maxwell count |E^K| = j|V^K| - C(j+1,2) for all subgroups K <= Aut(D).
"""
import itertools, json, sys
from collections import Counter
import networkx as nx
from carrier_core import analyse


def graphs(n, m):
    return [sorted(tuple(sorted(e)) for e in g.edges()) for g in nx.graph_atlas_g()
            if g.number_of_nodes() == n and g.number_of_edges() == m]


def aut_graph(n, edges):
    E = set(edges); deg = Counter()
    for a, b in E:
        deg[a] += 1; deg[b] += 1
    res = []
    for p in itertools.permutations(range(n)):
        if any(deg[p[i]] != deg[i] for i in range(n)):
            continue
        if all(tuple(sorted((p[a], p[b]))) in E for a, b in E):
            res.append(p)
    return res


def orientation_reps(edges, autG):
    seen = set(); reps = []
    for mask in range(1 << len(edges)):
        if mask in seen:
            continue
        arcs = [(a, b) if not (mask >> i) & 1 else (b, a) for i, (a, b) in enumerate(edges)]
        orb = set()
        for p in autG:
            s = {(p[x], p[y]) for x, y in arcs}
            orb.add(sum(1 << i for i, (a, b) in enumerate(edges) if (b, a) in s))
        seen |= orb; reps.append(arcs)
    return reps


def main():
    cases = [(2, 4), (2, 5), (2, 6), (3, 5), (3, 6)]
    if len(sys.argv) > 1:
        cases = [tuple(map(int, c.split(','))) for c in sys.argv[1:]]
    out = {}
    for j, n in cases:
        m = j * n - j * (j + 1) // 2
        tally = Counter(); bad = []
        for edges in graphs(n, m):
            autG = aut_graph(n, edges)
            for arcs in orientation_reps(edges, autG):
                r = analyse(n, arcs, j, want_weak=True)
                key = (r['linear'], r['maxwell'], r['weak'])
                tally[key] += 1
                if r['linear'] != r['maxwell']:
                    bad.append({'arcs': arcs, **r})
        out[f'j={j},n={n}'] = {'edges': m, 'tally(linear,maxwell,weak)': {str(k): v for k, v in tally.items()},
                               'disagreements': bad}
        print(f'j={j} n={n} ({m} edges): (linear, maxwell, weak) -> {dict(tally)}; disagreements: {len(bad)}', flush=True)
    json.dump(out, open('maxwell_ladder.json', 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
