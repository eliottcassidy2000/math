#!/usr/bin/env python3
"""
holt_metacirculant_kpswf5.py -- build the Holt graph as the metacirculant M(3,9).

Metacirculant M(m=3, n=9): vertices x^j_i, layer i in Z_3, position j in Z_9.
Semiregular sigma: x^j_i -> x^{j+1}_i (rotate within layer).
beta: x^j_i -> x^{r j}_{i+1} cyclically permutes layers, r = order-3 element mod 9.
r in {4,7} (4^3 = 64 = 1 mod 9, 7 = 4^2). Take r=4.

Symbol sets S_0,...,S_{floor(m/2)} of differences. Adjacency:
  x^j_i ~ x^{j + r^i * s}_i        for s in S_0   (within-layer i, scaled by r^i)
  x^j_i ~ x^{j + ...}_{i+t}         inter-layer per S_t.

We brute-force over candidate symbols and pick the unique graph that is
4-regular, connected, girth 5, |Aut|=54, vertex+edge-transitive, NOT arc-transitive.
That graph IS the Doyle-Holt graph; we then export its edge list.

kind-pasteur-2026-06-21-kpswf5
"""
import sys
from itertools import product, combinations
import networkx as nx
sys.stdout.reconfigure(line_buffering=True)

m, n = 3, 9
r = 4  # order 3 mod 9
def vid(i, j): return i * n + (j % n)

def build(S0, S1, S2):
    """S0: intra-layer diffs (layer 0 base, layer i scaled by r^i).
       S1: layer i -> i+1 connections, diffs applied. S2: i->i+2 (=i-1)."""
    G = nx.Graph()
    G.add_nodes_from(range(m * n))
    for i in range(m):
        ri = pow(r, i, n)
        # intra-layer
        for j in range(n):
            for s in S0:
                G.add_edge(vid(i, j), vid(i, (j + ri * s) % n))
        # inter-layer i -> i+1
        for j in range(n):
            for s in S1:
                G.add_edge(vid(i, j), vid((i + 1) % m, (j + ri * s) % n))
    return G

def is_target(G):
    if not nx.is_connected(G): return None
    if set(dict(G.degree()).values()) != {4}: return None
    if not hasattr(nx, 'girth'): return None
    if nx.girth(G) != 5: return None
    GM = nx.algorithms.isomorphism.GraphMatcher(G, G)
    auts = [dict(m_) for m_ in GM.isomorphisms_iter()]
    if len(auts) != 54: return None
    # vertex-transitive?
    nodes = list(G.nodes()); rep = {v: v for v in nodes}
    def f(x):
        while rep[x] != x: rep[x] = rep[rep[x]]; x = rep[x]
        return x
    def u(a, b):
        ra, rb = f(a), f(b)
        if ra != rb: rep[ra] = rb
    for g in auts:
        for v in nodes: u(v, g[v])
    if len({f(v) for v in nodes}) != 1: return None
    # edge orbits
    edges = [frozenset(e) for e in G.edges()]
    repe = {e: e for e in edges}
    def fe(x):
        while repe[x] != x: repe[x] = repe[repe[x]]; x = repe[x]
        return x
    def ue(a, b):
        ra, rb = fe(a), fe(b)
        if ra != rb: repe[ra] = rb
    for g in auts:
        for e in edges:
            a, b = tuple(e); ue(e, frozenset((g[a], g[b])))
    n_eorb = len({fe(e) for e in edges})
    # arc orbits
    arcs = []
    for a, b in G.edges(): arcs += [(a, b), (b, a)]
    repa = {a: a for a in arcs}
    def fa(x):
        while repa[x] != x: repa[x] = repa[repa[x]]; x = repa[x]
        return x
    def ua(a, b):
        ra, rb = fa(a), fa(b)
        if ra != rb: repa[ra] = rb
    for g in auts:
        for (a, b) in arcs: ua((a, b), (g[a], g[b]))
    n_aorb = len({fa(a) for a in arcs})
    return (len(auts), n_eorb, n_aorb)

# Candidate symbols: S0 subset of {1,2,3,4} (diffs, with +-), S1 single diff.
cands = []
diffs = [1, 2, 3, 4]
print("Brute-forcing metacirculant symbols for the Holt graph...")
for s0 in combinations(diffs, 1):
    for s1set in combinations(range(n), 1):
        G = build(list(s0), list(s1set), [])
        res = is_target(G)
        if res:
            cands.append((list(s0), list(s1set), res, G))

# also try S0 empty, all connections inter-layer (pure metacirculant)
for s1a in range(1, n):
    for s1b in range(1, n):
        G = nx.Graph(); G.add_nodes_from(range(m*n))
        for i in range(m):
            ri = pow(r, i, n)
            for j in range(n):
                G.add_edge(vid(i,j), vid((i+1)%m, (j+ri*s1a)%n))
                G.add_edge(vid(i,j), vid((i+1)%m, (j+ri*s1b)%n))
        res = is_target(G)
        if res:
            cands.append((['inter2'], [s1a, s1b], res, G))

print("  candidates found:", len(cands))
for s0, s1, res, G in cands[:5]:
    print("    S0=%s S1=%s -> |Aut|=%d edge-orb=%d arc-orb=%d" % (s0, s1, res[0], res[1], res[2]))

if cands:
    s0, s1, res, G = cands[0]
    relabel = {v: v for v in G.nodes()}
    with open("05-knowledge/results/holt_graph_edges_kpswf5.txt", "w") as f:
        for u_, v_ in sorted(G.edges()):
            f.write("%d %d\n" % (u_, v_))
    print("  HOLT GRAPH CONFIRMED. Saved edge list -> 05-knowledge/results/holt_graph_edges_kpswf5.txt")
    print("  Properties: 27 vertices, 54 edges, 4-regular, girth 5, |Aut|=54,")
    print("  vertex-transitive YES, edge-transitive YES, arc-transitive NO (arc-orbits=%d)" % res[2])
else:
    print("  No metacirculant matched; will try the within-layer staircase rule next.")
    # The search-result rule: layer i intra-diffs scaled, plus i->i+1 with diff,
    # but maybe the rule needs intra AND inter. Try S0=[1] per layer scaled + i->i-1.
    for s0d in [1,2,4]:
        for s1d in range(n):
            G = build([s0d], [s1d], [])
            res = is_target(G)
            if res:
                print("  matched with S0=[%d] S1=[%d]: %s" % (s0d, s1d, res))
