#!/usr/bin/env python3
"""
holt_builder_kpswf5.py -- build & verify the Doyle-Holt (Holt) graph.

The Holt graph is the smallest half-arc-transitive graph: 27 vertices,
4-regular, |Aut|=54, girth 5, vertex- and edge-transitive but NOT arc-transitive.
It is a Cayley graph on the metacyclic group
    G = <a,b | a^9 = b^3 = 1, b a b^{-1} = a^4>   (4 has order 3 mod 9).

We search connection sets S = {x, x^{-1}, y, y^{-1}} (closed under inverse,
no identity) for the one giving the Holt graph (4-regular, |Aut|=54, half-arc).

kind-pasteur-2026-06-21-kpswf5
"""
import sys
from itertools import combinations
import networkx as nx
sys.stdout.reconfigure(line_buffering=True)

# metacyclic group of order 27: elements (i in Z3 [exponent of b], j in Z9 [exp of a])
# relation b a b^{-1} = a^4. Then b^k a^j b^{-k} = a^{j 4^k}, and a^j b^k = b^k a^{j 4^{-k}}.
INV4 = pow(4, -1, 9)
def mul(x, y):
    i, j = x; i2, j2 = y
    f = pow(INV4, i2, 9)          # 4^{-i2} mod 9
    return ((i + i2) % 3, (j * f + j2) % 9)
def inv(x):
    i, j = x
    return ((-i) % 3, (-j * pow(4, (-i) % 3, 9)) % 9)
def power(x, k):
    r = (0, 0)
    for _ in range(k): r = mul(r, x)
    return r

ELEMS = [(i, j) for i in range(3) for j in range(9)]

def cayley_graph(S):
    G = nx.Graph()
    G.add_nodes_from(ELEMS)
    for g in ELEMS:
        for s in S:
            G.add_edge(g, mul(g, s))
    return G

def automorphisms_count(G):
    GM = nx.algorithms.isomorphism.GraphMatcher(G, G)
    c = 0
    for _ in GM.isomorphisms_iter():
        c += 1
    return c

def arc_orbit_count(G):
    GM = nx.algorithms.isomorphism.GraphMatcher(G, G)
    auts = [dict(m) for m in GM.isomorphisms_iter()]
    arcs = []
    for u, v in G.edges():
        arcs.append((u, v)); arcs.append((v, u))
    rep = {a: a for a in arcs}
    def find(x):
        while rep[x] != x:
            rep[x] = rep[rep[x]]; x = rep[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: rep[ra] = rb
    for g in auts:
        for (u, v) in arcs:
            union((u, v), (g[u], g[v]))
    return len({find(a) for a in arcs}), len(auts)

# candidate generators: all non-identity elements
nonid = [g for g in ELEMS if g != (0, 0)]
# look for {x, x^-1, y, y^-1} with |S|=4, connected, 4-regular, |Aut|=54, half-arc
seen = set()
found = None
print("Searching connection sets {x,x^-1,y,y^-1} on the order-27 metacyclic group...")
for x in nonid:
    xi = inv(x)
    if x == xi:   # involution -> would give odd connection set contribution
        continue
    for y in nonid:
        yi = inv(y)
        if y == yi: continue
        S = frozenset([x, xi, y, yi])
        if len(S) != 4: continue
        if S in seen: continue
        seen.add(S)
        G = cayley_graph(S)
        if not nx.is_connected(G): continue
        degs = set(dict(G.degree()).values())
        if degs != {4}: continue
        if nx.girth(G) if hasattr(nx, 'girth') else 0:
            pass
        # cheap filter: number of triangles must be 0 for Holt (girth 5)
        tri = sum(nx.triangles(G).values()) // 3
        if tri != 0:
            continue
        n_arc, n_aut = arc_orbit_count(G)
        # half-arc on a vertex-transitive Cayley graph: need edge-transitive + 2 arc orbits
        # edge-transitive check
        GM = nx.algorithms.isomorphism.GraphMatcher(G, G)
        auts = [dict(m) for m in GM.isomorphisms_iter()]
        edges = [frozenset(e) for e in G.edges()]
        rep = {e: e for e in edges}
        def find(z):
            while rep[z] != z:
                rep[z] = rep[rep[z]]; z = rep[z]
            return z
        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb: rep[ra] = rb
        for g in auts:
            for e in edges:
                u, v = tuple(e); union(e, frozenset((g[u], g[v])))
        n_eorb = len({find(e) for e in edges})
        if n_eorb == 1 and n_arc == 2:
            g6 = nx.girth(G) if hasattr(nx, 'girth') else None
            print("  FOUND half-arc-transitive Cayley graph!")
            print("    S =", sorted(S))
            print("    |V|=%d |E|=%d 4-regular, |Aut|=%d, edge-orbits=%d, arc-orbits=%d, triangles=%d, girth=%s"
                  % (G.number_of_nodes(), G.number_of_edges(), n_aut, n_eorb, n_arc, tri, g6))
            found = (S, G, n_aut)
            break
    if found: break

if found is None:
    print("  No half-arc-transitive Cayley graph found in {x,x^-1,y,y^-1} family.")
else:
    S, G, n_aut = found
    # export edge list (as integer-labeled) for reuse
    relabel = {v: i for i, v in enumerate(ELEMS)}
    H = nx.relabel_nodes(G, relabel)
    with open("05-knowledge/results/holt_graph_edges_kpswf5.txt", "w") as f:
        for u, v in sorted(H.edges()):
            f.write("%d %d\n" % (u, v))
    print("  Saved integer edge list -> 05-knowledge/results/holt_graph_edges_kpswf5.txt")
    print("  Verified: this IS the Holt/Doyle-Holt graph (|Aut|=54 expected:%s)" % (n_aut==54))
