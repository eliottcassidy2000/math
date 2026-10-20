#!/usr/bin/env python3
"""HYP-9340 first open case, H = Z_3 of type (3,1^{n-3}). For each graph type X on n vertices:
 e(X) = [X even and X has an independent twin triple];
 b(X) = sum_{j>=1} (-1)^{j+1} #{Aut(X)-orbits of j-sets J of pairwise disjoint independent twin triples with the contracted
        sign s'_J trivial on Stab(J)}, s'_J(a) = sgn_X(a) * (-1)^{sum over adjacent twin-triple pairs swapped by a of
        (fix(a^2 on the triple) - 1)}.
(s' computed by contracting the marked triples) THM-4721-style argument gives sum_X b(X) = tau_T(Z_3) exactly; HYP-9340 for this H says sum_X (b - e) = 0.
Prints every type with b != e."""
import sys, subprocess, itertools
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
from collections import defaultdict
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
def sign_of(edges, g):
    s = 0
    for u, v in edges:
        if g[u] > g[v]: s ^= 1
    return -1 if s else 1
for n in range(5, NMAX + 1):
    out = subprocess.run(['geng', '-q', str(n)], capture_output=True, text=True).stdout.split()
    tot_b = tot_e = 0; defects = []
    for s6 in out:
        G = nx.from_graph6_bytes(s6.encode()); edges = [(min(u, v), max(u, v)) for u, v in G.edges()]
        nb = {v: frozenset(G[v]) for v in range(n)}
        # independent twin triples
        classes = defaultdict(list)
        for v in range(n): classes[nb[v]].append(v)
        triples = [frozenset(t) for C in classes.values() if len(C) >= 3 for t in itertools.combinations(C, 3)]
        if not triples: continue
        auts = [tuple(m[i] for i in range(n)) for m in GraphMatcher(G, G).isomorphisms_iter()]
        sg = {a: sign_of(edges, a) for a in auts}
        even = all(v == 1 for v in sg.values())
        e = 1 if even else 0
        # matchings of disjoint triples
        matchings = []
        def rec(start, cur, used):
            for i in range(start, len(triples)):
                t = triples[i]
                if t & used: continue
                matchings.append(frozenset(cur + [t])); rec(i + 1, cur + [t], used | t)
        rec(0, [], frozenset())
        adj = lambda t1, t2: any(G.has_edge(u, v) for u in t1 for v in t2)
        def act(a, J): return frozenset(frozenset(a[x] for x in t) for t in J)
        b = 0; seen = set()
        for J in matchings:
            if J in seen: continue
            orbit = {act(a, J) for a in auts}; seen |= orbit
            stab = [a for a in auts if act(a, J) == J]
            # contracted graph: marked triples -> their min vertex; s'(a) = sgn of the induced automorphism on it
            rep = {}
            for t in J:
                for x in t: rep[x] = min(t)
            cv = [v for v in range(n) if v not in rep or rep[v] == v]
            cedges = set()
            for u, v in edges:
                cu, cvv = rep.get(u, u), rep.get(v, v)
                if cu != cvv: cedges.add((min(cu, cvv), max(cu, cvv)))
            trivial = True
            for a in stab:
                alpha = {}
                for v in cv:
                    img = a[v]; alpha[v] = rep.get(img, img)
                sp = 0
                for u, v in cedges:
                    if alpha[u] > alpha[v]: sp ^= 1
                if sp: trivial = False; break
            if trivial: b += (-1) ** (len(J) + 1)
        tot_b += b; tot_e += e
        if b != e: defects.append((s6, len(edges), len(auts), even, b, [len(C) for C in classes.values() if len(C) >= 3]))
    print(f"n = {n}: sum b = {tot_b}, sum e = {tot_e}; types with b != e: {len(defects)}")
    for d in defects: print("   ", d)
