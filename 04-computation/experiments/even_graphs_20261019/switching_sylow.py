#!/usr/bin/env python3
"""Switching-class half of HYP-9172, numerical probe for odd n (THM-4531 P1: switching classes of tournaments on odd n
correspond S_n-equivariantly to tournaments with all scores = (n-1)/2 mod 2). For odd-order P <= S_n:
 tau_C(P)  = # types of mod-4-balanced tournaments with Aut containing a conjugate of P,
 tau_EE(P) = # types of even Euler graphs (all degrees even, no odd automorphism) with Aut containing a conjugate of P.
Both by enumerating the P-invariant labelled objects (unions of pair orbits) and canonicalising with labelg."""
import sys, subprocess, itertools
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
sys.path.insert(0, '.')
from subgroups_fix import cyc, mul, closure, pair_orbits, graph6, labelg, to_d6, subgroup_reps
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
for n in range(5, NMAX + 1, 2):
    out = subprocess.run(['geng', '-q', str(n)], capture_output=True, text=True).stdout.split()
    canon = labelg(out); even_of = {}
    for s, c in zip(out, canon):
        G = nx.from_graph6_bytes(s.encode()); edges = [(min(u, v), max(u, v)) for u, v in G.edges()]
        ev = True
        for m in GraphMatcher(G, G).isomorphisms_iter():
            t = 0
            for u, v in edges:
                if m[u] > m[v]: t ^= 1
            if t: ev = False; break
        even_of[c] = ev
    print(f"\n=== n = {n} (scores must be = {(n-1)//2 % 2} mod 2) ===")
    for name, H in subgroup_reps(n).items():
        orbs = pair_orbits(H, n)
        orb_arcs = []
        for O in orbs:
            p = O[0]; orb_arcs.append(frozenset((g[p[0]], g[p[1]]) for g in H))
        tlines = []; glines = []
        for choice in itertools.product([0, 1], repeat=len(orbs)):
            arcs = []; edges = []
            for c, O, A in zip(choice, orbs, orb_arcs):
                if c: arcs += list(A); edges += O
                else: arcs += [(v, u) for u, v in A]
            outdeg = [0] * n
            for u, v in arcs: outdeg[u] += 1
            if all(d % 2 == (n - 1) // 2 % 2 for d in outdeg): tlines.append(to_d6(n, arcs))
            deg = [0] * n
            for u, v in edges: deg[u] += 1; deg[v] += 1
            if all(d % 2 == 0 for d in deg): glines.append(graph6(n, edges))
        tc = set(labelg(tlines)) if tlines else set()
        gc = set(labelg(glines)) if glines else set()
        ee = sum(1 for c in gc if even_of[c])
        flag = "" if len(tc) == ee else "  <-- DIFFER"
        print(f"  P = {name} (|P|={len(H)}): tau_C = {len(tc)}, tau_EE = {ee}, Euler types {len(gc)}{flag}")
