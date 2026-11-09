"""AUDIT2: cross-check of a2_excess.py's list of connected graphs with excess <= 5, using the networkx graph atlas
(all 1253 graphs on <= 7 vertices; an independent enumeration source), and of the co-evenness of each."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a2_lib import *
import networkx as nx
from networkx.generators.atlas import graph_atlas_g
def excess(adj):
    seen = set()
    for N in adj:
        sub = N
        while True:
            if bin(sub).count('1') >= 2: seen.add(sub)
            if sub == 0: break
            sub = (sub - 1) & N
    return len(seen)
rows = []
for G in graph_atlas_g():
    n = G.number_of_nodes()
    if n < 3 or not nx.is_connected(G): continue
    g = from_edges(n, list(G.edges()))
    e = excess(g[1])
    if e <= 5:
        full = (1 << n) - 1
        comp = (n, [(full ^ a) & ~(1 << v) for v, a in enumerate(g[1])])
        co = run_check([(comp, -1)])[0][1] == 'EVEN'
        rows.append((e, n, G.number_of_edges(), co, sorted(d for _, d in G.degree())))
rows.sort()
for r in rows: print("e=%d order=%d edges=%d co-even=%s degrees=%s" % r)
print("total connected graphs (order >= 3) with excess <= 5:", len(rows))
# largest order with excess <= 5 in the atlas is 7 (P_7); the growth lemma excludes order >= 8
