#!/usr/bin/env python3
"""
even_graph_invariants_s315.py — Deep invariants of E_n
opus-2026-03-24-S315

KEY QUESTION: Why is E_n PERFECT (χ = ω)?

APPROACH: Check if E_n is a comparability graph, or interval graph,
or has some ordering property that forces perfectness.

Also compute even graph analogues of tournament invariants.
"""

import sys
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def build_even_graph_metagraph(n):
    edges = [(i,j) for i in range(n) for j in range(i+1, n)]
    m_full = len(edges)
    edge_idx = {e: i for i, e in enumerate(edges)}
    base_edges = [(i, i+1) for i in range(n-1)]
    non_base = [e for e in edges if e not in base_edges]
    m = len(non_base)
    
    fund_cycles = []
    for e in non_base:
        i, j = e
        cycle_bits = 0
        for k in range(i, j): cycle_bits |= (1 << edge_idx[(k, k+1)])
        cycle_bits |= (1 << edge_idx[(i, j)])
        fund_cycles.append(cycle_bits)
    
    def tiling_to_even(tile_bits):
        result = 0
        for k in range(m):
            if tile_bits & (1 << k): result ^= fund_cycles[k]
        return result
    
    def graph_edge_count(g): return bin(g).count('1')
    
    all_perms = list(permutations(range(n)))
    def graph_canon(g_bits):
        best = None
        for p in all_perms:
            nb = 0
            for k, (i,j) in enumerate(edges):
                pi, pj = min(p[i], p[j]), max(p[i], p[j])
                nk = edge_idx[(pi, pj)]
                if g_bits & (1 << k): nb |= (1 << nk)
            if best is None or nb < best: best = nb
        return best
    
    class_map = {}; cid = 0; class_info = {}
    tiling_class = {}
    
    for tb in range(1 << m):
        eg = tiling_to_even(tb)
        cn = graph_canon(eg)
        if cn not in class_map:
            class_map[cn] = cid
            ec = graph_edge_count(eg)
            # Degree sequence
            deg = [0] * n
            for k, (i,j) in enumerate(edges):
                if eg & (1 << k): deg[i] += 1; deg[j] += 1
            ds = tuple(sorted(deg))
            # Triangles
            tri = 0
            for i in range(n):
                for j in range(i+1, n):
                    for k in range(j+1, n):
                        if (eg & (1<<edge_idx[(i,j)])) and (eg & (1<<edge_idx[(i,k)])) and (eg & (1<<edge_idx[(j,k)])):
                            tri += 1
            # Components
            adj_g = defaultdict(set)
            for k, (i,j) in enumerate(edges):
                if eg & (1 << k): adj_g[i].add(j); adj_g[j].add(i)
            visited = set(); comps = 0
            for v in range(n):
                if v not in visited:
                    comps += 1; stack = [v]
                    while stack:
                        u = stack.pop()
                        if u in visited: continue
                        visited.add(u)
                        for w in adj_g[u]:
                            if w not in visited: stack.append(w)
            
            class_info[cid] = {'canon': cn, 'edges': ec, 'deg_seq': ds, 'tri': tri, 'comp': comps, 'fiber': 0}
            cid += 1
        
        class_info[class_map[cn]]['fiber'] += 1
        tiling_class[tb] = class_map[cn]
    
    V = cid
    meta_edges = set()
    for tb in range(1 << m):
        ca = tiling_class[tb]
        for tile in range(m):
            cb = tiling_class[tb ^ (1 << tile)]
            if ca != cb: meta_edges.add((min(ca, cb), max(ca, cb)))
    
    adj_list = defaultdict(set)
    for (a, b) in meta_edges: adj_list[a].add(b); adj_list[b].add(a)
    
    return V, class_info, meta_edges, adj_list, tiling_class, m, edge_idx, edges

print("=" * 72)
print("  EVEN GRAPH INVARIANTS AND PERFECTNESS")
print("  opus-2026-03-24-S315")
print("=" * 72)

for n in [5, 6]:
    print(f"\n{'#'*72}")
    print(f"  n = {n}")
    print(f"{'#'*72}")
    
    V, ci, edges, adj, tc, m, eidx, all_edges = build_even_graph_metagraph(n)
    E = len(edges)
    print(f"  V(E_{n}) = {V}, E(E_{n}) = {E}")
    
    # 1. IS E_n A COMPARABILITY GRAPH?
    # A comparability graph = the edges are a transitive relation.
    # Check: does the edge set form a partial order (after orientation)?
    # For a perfect graph, the complement is also perfect.
    # Check complement: omega(complement) = alpha(E_n)
    
    comp_edges = set()
    for a in range(V):
        for b in range(a+1, V):
            if (a, b) not in edges:
                comp_edges.add((a, b))
    print(f"\n  COMPLEMENT: V={V}, E(complement)={len(comp_edges)}")
    
    # Find alpha(E_n) = omega(complement)
    comp_adj = defaultdict(set)
    for (a, b) in comp_edges: comp_adj[a].add(b); comp_adj[b].add(a)
    
    # Greedy max clique on complement = max independent set of E_n
    best_clique = []
    for start in range(V):
        clique = [start]
        candidates = sorted(comp_adj[start], key=lambda v: -len(comp_adj[v]))
        for c in candidates:
            if all((min(c, x), max(c, x)) in comp_edges for x in clique):
                clique.append(c)
        if len(clique) > len(best_clique):
            best_clique = clique
    
    alpha = len(best_clique)
    print(f"  α(E_{n}) = {alpha} (independent set)")
    print(f"  Independent set: {best_clique}")
    print(f"  Their edge counts: {[ci[v]['edges'] for v in best_clique]}")
    
    # 2. CHECK: Is the EDGE COUNT a valid coloring?
    # If edges(G1) = edges(G2) and they're adjacent, this fails.
    edge_count_coloring = True
    for (a, b) in edges:
        if ci[a]['edges'] == ci[b]['edges']:
            edge_count_coloring = False
            break
    print(f"\n  Edge count as coloring? {edge_count_coloring}")
    
    # Edge count distribution
    ec_dist = Counter(ci[v]['edges'] for v in range(V))
    print(f"  Edge count distribution: {sorted(ec_dist.items())}")
    n_edge_counts = len(ec_dist)
    print(f"  Distinct edge counts: {n_edge_counts}")
    
    # 3. Check if EDGE COUNT determines an ordering that makes E_n a comparability graph
    # Two classes are adjacent → do their edge counts differ?
    same_ec_edges = sum(1 for (a,b) in edges if ci[a]['edges'] == ci[b]['edges'])
    print(f"  Same-edge-count edges: {same_ec_edges}/{E} ({100*same_ec_edges/E:.1f}%)")
    
    # 4. EDGE COUNT MOD k coloring
    for k in range(2, 8):
        valid = True
        for (a, b) in edges:
            if ci[a]['edges'] % k == ci[b]['edges'] % k:
                valid = False; break
        if valid:
            print(f"  Edge count mod {k}: VALID coloring ({len(set(ci[v]['edges'] % k for v in range(V)))} colors)")
    
    # 5. Is the complement of E_n triangle-free?
    comp_tri = 0
    for a in range(V):
        for b in comp_adj[a]:
            if b <= a: continue
            for c in comp_adj[a]:
                if c <= b: continue
                if (min(b,c), max(b,c)) in comp_edges:
                    comp_tri += 1
    print(f"\n  Complement triangles: {comp_tri}")
    if comp_tri == 0:
        print(f"  *** Complement is TRIANGLE-FREE! ***")
        print(f"  This means alpha(E_n) <= 2 and E_n has no independent set of size 3!")
    
    # 6. Is E_n chordal? (every cycle of length >= 4 has a chord)
    # A chordal graph is perfect. Check by looking for chordless 4-cycles.
    has_chordless_4 = False
    for a in range(V):
        for b in adj[a]:
            if b <= a: continue
            for c in adj[b]:
                if c <= a or c == a: continue
                if (min(a,c), max(a,c)) in edges: continue  # a-c is a chord
                for d in adj[c]:
                    if d <= a or d == b: continue
                    if (min(a,d), max(a,d)) not in edges: continue  # need a-d edge
                    if (min(b,d), max(b,d)) not in edges: continue  # need this for 4-cycle
                    # Wait, need to check: a-b, b-c, c-d, d-a is a 4-cycle
                    # with no chord a-c or b-d
                    if (min(a,c), max(a,c)) not in edges and (min(b,d), max(b,d)) not in edges:
                        has_chordless_4 = True
                        break
            if has_chordless_4: break
        if has_chordless_4: break
    
    print(f"  Chordal? {'No (has chordless 4-cycle)' if has_chordless_4 else 'YES!'}")
    
    # 7. Is there a NATURAL partial order making E_n a comparability graph?
    # Try: order by edge count (most natural)
    # In a comparability graph, u-v edge iff u < v in some partial order
    # If we orient edges by edge count (lower → higher), check transitivity
    print(f"\n  COMPARABILITY TEST (ordered by edge count):")
    transitive_violations = 0
    for a in range(V):
        for b in adj[a]:
            if ci[a]['edges'] >= ci[b]['edges']: continue  # skip reverse direction
            # a → b (a has fewer edges). Check: for every c adjacent to b with more edges,
            # is a also adjacent to c?
            for c in adj[b]:
                if ci[b]['edges'] >= ci[c]['edges']: continue
                if (min(a,c), max(a,c)) not in edges:
                    transitive_violations += 1
    print(f"    Transitivity violations: {transitive_violations}")

print("\nDONE.")
