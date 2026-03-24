#!/usr/bin/env python3
"""
even_odd_shadow_s20du.py — Even/odd graph deep dive through the metagraph lens
kind-pasteur-2026-03-23-S20du

The merged metagraph G_n/Z_2 applies to BOTH tournaments and even graphs
(since they are Burnside-equivalent). This script explores:

1. The complete even/odd partition of all graphs on n vertices
2. What happens to even graphs under edge flip (the graph metagraph)
3. The "shadow" of the tournament metagraph onto graph space
4. How odd graphs relate to tournaments (the "inverse" relationship)
5. The even graph metagraph vs the tournament metagraph
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  EVEN/ODD GRAPHS AND THE SHADOW OF THE METAGRAPH")
print("  kind-pasteur-2026-03-23-S20du")
print("=" * 80)

for n in [3, 4, 5]:
    t0 = time.time()
    ALL_EDGES = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(ALL_EDGES)
    perms = list(permutations(range(n)))

    def graph_adj(bits):
        A = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(ALL_EDGES):
            if bits & (1 << k): A[i][j] = 1; A[j][i] = 1
        return A

    def canon_graph(bits):
        best = None
        for p in perms:
            nb = 0
            for k, (i,j) in enumerate(ALL_EDGES):
                pi, pj = min(p[i],p[j]), max(p[i],p[j])
                nk = ALL_EDGES.index((pi, pj))
                if (bits >> k) & 1: nb |= (1 << nk)
            if best is None or nb < best: best = nb
        return best

    def graph_auts(A):
        return [p for p in perms if all(A[p[i]][p[j]] == A[i][j]
                for i in range(n) for j in range(n))]

    def is_odd_graph(bits):
        A = graph_adj(bits)
        edge_list = [(i,j) for k, (i,j) in enumerate(ALL_EDGES) if bits & (1 << k)]
        if not edge_list: return False
        auts = graph_auts(A)
        num_e = len(edge_list)
        for orient in range(1 << num_e):
            oriented = []
            for k, (i,j) in enumerate(edge_list):
                if orient & (1 << k): oriented.append((j,i))
                else: oriented.append((i,j))
            for sigma in auts:
                if sigma == list(range(n)): continue
                rev = 0
                for u, v in oriented:
                    su, sv = sigma[u], sigma[v]
                    e = (min(su,sv), max(su,sv))
                    idx = edge_list.index(e)
                    ou, ov = oriented[idx]
                    if (su == ov and sv == ou): rev += 1
                if rev % 2 == 1: return True
        return False

    # Build all graph iso classes with even/odd classification
    graph_classes = {}
    for bits in range(1 << m):
        cn = canon_graph(bits)
        if cn not in graph_classes:
            A = graph_adj(bits)
            deg = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
            num_e = bin(bits).count('1')
            aut = len(graph_auts(A))
            odd = is_odd_graph(bits)
            comp_bits = ((1 << m) - 1) ^ bits
            comp_cn = canon_graph(comp_bits)
            graph_classes[cn] = {
                'bits': bits, 'deg': deg, 'edges': num_e,
                'aut': aut, 'odd': odd, 'even': not odd,
                'comp_cn': comp_cn
            }

    # Build graph metagraph (single edge flip)
    graph_edges = set()
    for cn, info in graph_classes.items():
        bits = info['bits']
        for k in range(m):
            bits2 = bits ^ (1 << k)
            cn2 = canon_graph(bits2)
            if cn2 != cn:
                graph_edges.add((min(cn, cn2), max(cn, cn2)))

    even_classes = {cn: info for cn, info in graph_classes.items() if info['even']}
    odd_classes = {cn: info for cn, info in graph_classes.items() if info['odd']}

    print(f"\n{'#'*60}")
    print(f"  n = {n}")
    print(f"{'#'*60}")
    print(f"  Total graphs: {len(graph_classes)}, Even: {len(even_classes)}, Odd: {len(odd_classes)}")
    print(f"  Graph metagraph edges: {len(graph_edges)}")

    # Classify metagraph edges by even/odd endpoint types
    ee_edges = 0; eo_edges = 0; oo_edges = 0
    for a, b in graph_edges:
        a_even = graph_classes[a]['even']
        b_even = graph_classes[b]['even']
        if a_even and b_even: ee_edges += 1
        elif a_even != b_even: eo_edges += 1
        else: oo_edges += 1

    print(f"\n  GRAPH METAGRAPH EDGE TYPES:")
    print(f"    Even-Even: {ee_edges}")
    print(f"    Even-Odd: {eo_edges}")
    print(f"    Odd-Odd: {oo_edges}")

    # Self-loops in graph metagraph (edge flip preserving class)
    graph_sl = 0
    for cn, info in graph_classes.items():
        bits = info['bits']
        for k in range(m):
            bits2 = bits ^ (1 << k)
            if canon_graph(bits2) == cn:
                graph_sl += 1
    # Divide by orbit size
    print(f"  Self-loops (total raw): {graph_sl}")

    # Compare with tournament metagraph
    def make_tourn(bits):
        A = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(ALL_EDGES):
            if bits & (1 << k): A[i][j] = 1
            else: A[j][i] = 1
        return A

    def canon_tourn(A):
        best = None
        for p in perms:
            s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best: best = s
        return best

    tourn_classes = {}
    for bits in range(1 << m):
        A = make_tourn(bits)
        cn = canon_tourn(A)
        if cn not in tourn_classes:
            tourn_classes[cn] = {'bits': bits}

    tourn_edges = set()
    for cn, info in tourn_classes.items():
        bits = info['bits']
        for k in range(m):
            bits2 = bits ^ (1 << k)
            A2 = make_tourn(bits2)
            cn2 = canon_tourn(A2)
            if cn2 != cn:
                tourn_edges.add((min(cn, cn2), max(cn, cn2)))

    print(f"\n  TOURNAMENT METAGRAPH: V={len(tourn_classes)}, E={len(tourn_edges)}")
    print(f"  GRAPH METAGRAPH: V={len(graph_classes)}, E={len(graph_edges)}")

    # The SHADOW: how does the even subgraph relate?
    print(f"\n  EVEN GRAPH SUB-METAGRAPH:")
    even_internal = sum(1 for a,b in graph_edges
                       if graph_classes[a]['even'] and graph_classes[b]['even'])
    print(f"    V={len(even_classes)}, E(internal)={even_internal}")
    print(f"    Tournament metagraph: V={len(tourn_classes)}, E={len(tourn_edges)}")
    print(f"    Same V? {len(even_classes) == len(tourn_classes)}")

    # Complement structure
    print(f"\n  COMPLEMENT STRUCTURE:")
    even_comp_even = sum(1 for cn in even_classes
                        if graph_classes[graph_classes[cn]['comp_cn']]['even'])
    even_comp_odd = len(even_classes) - even_comp_even
    odd_comp_even = sum(1 for cn in odd_classes
                       if graph_classes[graph_classes[cn]['comp_cn']]['even'])
    odd_comp_odd = len(odd_classes) - odd_comp_even
    print(f"    Even -> Even complement: {even_comp_even}")
    print(f"    Even -> Odd complement: {even_comp_odd}")
    print(f"    Odd -> Even complement: {odd_comp_even}")
    print(f"    Odd -> Odd complement: {odd_comp_odd}")

    # Detail: list all even graphs with their properties
    print(f"\n  EVEN GRAPH DETAILS:")
    print(f"  {'edges':>5} {'deg':>20} {'Aut':>4} {'comp_type':>10}")
    for cn in sorted(even_classes.keys(), key=lambda c: even_classes[c]['edges']):
        info = even_classes[cn]
        comp_type = 'even' if graph_classes[info['comp_cn']]['even'] else 'odd'
        print(f"  {info['edges']:5d} {str(info['deg']):>20} {info['aut']:4d} {comp_type:>10}")

    print(f"\n  ODD GRAPH DETAILS:")
    print(f"  {'edges':>5} {'deg':>20} {'Aut':>4} {'comp_type':>10}")
    for cn in sorted(odd_classes.keys(), key=lambda c: odd_classes[c]['edges']):
        info = odd_classes[cn]
        comp_type = 'even' if graph_classes[info['comp_cn']]['even'] else 'odd'
        print(f"  {info['edges']:5d} {str(info['deg']):>20} {info['aut']:4d} {comp_type:>10}")

    print(f"\n  Time: {time.time()-t0:.1f}s")

print(f"\n{'='*60}")
print("SYNTHESIS")
print(f"{'='*60}")
print("""
THE SHADOW OF THE METAGRAPH:

The graph metagraph (ALL graphs on n vertices, edges = single edge flips)
contains BOTH even and odd graphs. The even/odd partition creates
THREE types of metagraph edges:

  Even-Even (EE): both endpoints are even graphs
  Even-Odd (EO): one even, one odd
  Odd-Odd (OO): both endpoints are odd graphs

The EVEN SUB-METAGRAPH (EE edges only) is the "shadow" of the
tournament metagraph onto graph space. It has the SAME number of
vertices as the tournament metagraph (by Royle et al.) but
potentially DIFFERENT edge structure.

The ODD graphs are the "anti-tournaments" — they are what's LEFT
after removing the tournament-equivalent objects from graph space.
The EO edges are the BOUNDARY between the tournament world and
its complement.

KEY QUESTION: Is the even sub-metagraph ISOMORPHIC to the tournament
metagraph? Or just equinumerous in vertices?
""")

print("DONE.")
print("=" * 80)
