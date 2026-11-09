"""AUDIT2: the connected realisers of d = 37 (all on 6 vertices by the census; <= 7 vertices by the counting lemma),
with their automorphism groups and the reason they are even."""
import itertools, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a2_lib import *
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
reps = []
for n in range(2, 8):
    pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
    if n == 7: break          # census: no connected realiser at n = 7 (a2_census_n7.out, iso_ceven[37] = 0)
    for mask in range(1 << len(pairs)):
        g = from_edges(n, [pairs[e] for e in range(len(pairs)) if mask >> e & 1])
        G = nx_graph(g)
        if not nx.is_connected(G): continue
        if d_brute(g) != 37: continue
        if any(nx.is_isomorphic(G, nx_graph(h)) for h in reps): continue
        reps.append(g)
for g in reps:
    G = nx_graph(g)
    auts = list(GraphMatcher(G, G).isomorphisms_iter())
    E = edges(g)
    signs = [sum(1 for u, v in E if p[u] > p[v]) % 2 for p in auts]
    degs = [bin(a).count('1') for a in g[1]]
    univ = any(a == (1 << g[0]) - 1 - (1 << v) for v, a in enumerate(g[1]))
    twins_nonadj = [(u, v) for u in range(g[0]) for v in range(u + 1, g[0]) if g[1][u] == g[1][v]]
    print("n=%d edges=%s degs=%s |Aut|=%d odd auts=%d universal=%s non-adjacent twin pairs=%s"
          % (g[0], E, degs, len(auts), sum(signs), univ, twins_nonadj))
