#!/usr/bin/env python3
"""Independent check of the C15 counterexample (session-10 packet, 2026-10-20): H = <(012)(34567)> <= S_10, orbits 3,5,1,1.
tau_T(H) = # tournament types with Aut containing a conjugate of H; tau_E(H) = same for even graphs.
Enumerate the 2^9 H-invariant labelled objects, canonicalise with labelg, and test evenness of each graph type
directly (networkx automorphisms; the empty graph is even by definition, complete graphs odd)."""
import sys, itertools, subprocess
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
sys.path.insert(0, '.')
from subgroups_fix import cyc, closure, pair_orbits, graph6, labelg, to_d6
n = 10
H = closure([cyc([[0, 1, 2], [3, 4, 5, 6, 7]], n)], n); assert len(H) == 15
orbs = pair_orbits(H, n); print("pair orbits:", len(orbs))
orb_arcs = []
for O in orbs:
    p = O[0]; orb_arcs.append(frozenset((g[p[0]], g[p[1]]) for g in H))
tlines = []; glines = []; gedges = {}
for choice in itertools.product([0, 1], repeat=len(orbs)):
    arcs = []; edges = []
    for c, O, A in zip(choice, orbs, orb_arcs):
        if c: arcs += list(A); edges += O
        else: arcs += [(v, u) for u, v in A]
    tlines.append(to_d6(n, arcs)); g6 = graph6(n, edges); glines.append(g6); gedges[g6] = edges
tcan = set(labelg(tlines)); gcan = labelg(glines)
types = {}
for g6, c in zip(glines, gcan): types.setdefault(c, g6)
print("tournament types:", len(tcan), "| graph types:", len(types))
even = 0; odd_list = []
for c, g6 in types.items():
    E = gedges[g6]
    if len(E) == 0: even += 1; continue
    if len(E) == n * (n - 1) // 2: odd_list.append(('K10', len(E))); continue
    G = nx.from_graph6_bytes(g6.encode())
    ev = True
    for m in GraphMatcher(G, G).isomorphisms_iter():
        if sum(1 for u, v in E if m[u] > m[v]) % 2: ev = False; break
    if ev: even += 1
    else: odd_list.append((g6, len(E)))
print(f"tau_T(C15) = {len(tcan)}, tau_E(C15) = {even}  (odd graph types among the invariant ones: {len(odd_list)})")
