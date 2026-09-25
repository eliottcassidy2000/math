#!/usr/bin/env python3
"""Independent audit of decoder_minors / seam_threshold / thirtysix_multipartite:
 - G_N (edges {x,x+2},{x,2x} on 1..N) planar iff N<=15 (networkx LR test, own graph
   construction; verify the K33 subdivision paths at N=16 directly);
 - Q_N (square-sum graph) planarity flags N=14..16 and first nonplanar N;
 - L(K5): perfect matchings = 144 = 24 (distinct centers) + 120; complement of
   L(K5) is the Petersen graph; 24 regular labelled 5-tournaments.
 - number of rank-2 simple graphs on 6 labelled vertices and the odd-degree ones.
"""
import itertools
import collections
import networkx as nx
from math import isqrt


def G(N):
    g = nx.Graph()
    g.add_nodes_from(range(1, N + 1))
    for x in range(1, N + 1):
        if x + 2 <= N:
            g.add_edge(x, x + 2)
        if 2 * x <= N:
            g.add_edge(x, 2 * x)
    return g


def Q(N):
    g = nx.Graph()
    g.add_nodes_from(range(1, N + 1))
    for x in range(1, N + 1):
        for y in range(x + 1, N + 1):
            s = x + y
            if isqrt(s) ** 2 == s:
                g.add_edge(x, y)
    return g


def main():
    flags = [(N, nx.check_planarity(G(N))[0]) for N in range(1, 41)]
    print("G_N planar flags:", "".join("P" if f else "n" for _, f in flags))
    print("G_N last planar N:", max(N for N, f in flags if f), " first nonplanar:", min(N for N, f in flags if not f))
    g16 = G(16)
    paths = [[5, 3, 6], [5, 10], [5, 7, 14], [8, 6], [8, 10], [8, 16, 14], [12, 6], [12, 10], [12, 14]]
    ok = all(g16.has_edge(a, b) for p in paths for a, b in zip(p, p[1:]))
    interiors = [v for p in paths for v in p[1:-1]]
    print("G16 K33 subdivision edges valid:", ok, " interiors", interiors, "disjoint:", len(set(interiors)) == len(interiors),
          " interiors avoid branch vertices:", not (set(interiors) & {5, 8, 12, 6, 10, 14}))
    g15 = G(15)
    print("G15: V=%d E=%d" % (g15.number_of_nodes(), g15.number_of_edges()))
    qf = [(N, nx.check_planarity(Q(N))[0]) for N in range(1, 61)]
    print("Q_N planar flags N=1..60:", "".join("P" if f else "n" for _, f in qf))
    print("Q_N first nonplanar:", min(N for N, f in qf if not f))

    # L(K5) and Petersen
    K5 = nx.complete_graph(5)
    L = nx.line_graph(K5)
    comp = nx.complement(L)
    print("L(K5): V=%d E=%d regular deg %s ; complement E=%d, isomorphic to Petersen: %s"
          % (L.number_of_nodes(), L.number_of_edges(), set(dict(L.degree()).values()),
             comp.number_of_edges(), nx.is_isomorphic(comp, nx.petersen_graph())))
    # count perfect matchings of L(K5) by own recursion
    nodes = list(L.nodes())
    adj = {v: set(L.neighbors(v)) for v in nodes}

    def count_pm(rem):
        if not rem:
            return [[]]
        v = min(rem)
        out = []
        for w in adj[v]:
            if w in rem:
                for m in count_pm(rem - {v, w}):
                    out.append([(v, w)] + m)
        return out
    pms = count_pm(frozenset(nodes))
    centers = {}
    for m in pms:
        cs = []
        for (e1, e2) in m:
            c = set(e1) & set(e2)
            cs.append(next(iter(c)))
        prof = tuple(sorted(collections.Counter(cs).values()))
        centers[prof] = centers.get(prof, 0) + 1
    print("perfect matchings of L(K5):", len(pms), " center-multiplicity profiles:", centers)
    # regular 5-tournaments
    pairs = list(itertools.combinations(range(5), 2))
    reg = 0
    for mask in range(1 << 10):
        out = [0] * 5
        for i, (a, b) in enumerate(pairs):
            if (mask >> i) & 1:
                out[a] += 1
            else:
                out[b] += 1
        if all(o == 2 for o in out):
            reg += 1
    print("regular labelled 5-tournaments:", reg)

    # rank-2 graphs on 6 labelled vertices
    pairs6 = list(itertools.combinations(range(6), 2))
    def rank_f2(rows):
        rows = list(rows)
        rk = 0
        for col in range(6):
            piv = next((i for i in range(rk, 6) if (rows[i] >> col) & 1), None)
            if piv is None:
                continue
            rows[rk], rows[piv] = rows[piv], rows[rk]
            for i in range(6):
                if i != rk and (rows[i] >> col) & 1:
                    rows[i] ^= rows[rk]
            rk += 1
        return rk
    r2 = 0
    r2odd = {}
    for mask in range(1 << 15):
        rows = [0] * 6
        for i, (a, b) in enumerate(pairs6):
            if (mask >> i) & 1:
                rows[a] |= 1 << b
                rows[b] |= 1 << a
        if rank_f2(rows) == 2:
            r2 += 1
            degs = [bin(r).count("1") for r in rows]
            if all(d % 2 == 1 for d in degs):
                key = tuple(sorted(degs))
                r2odd[key] = r2odd.get(key, 0) + 1
    print("rank-2 labelled simple graphs on 6 vertices:", r2, " with all degrees odd:", r2odd)


if __name__ == "__main__":
    main()
