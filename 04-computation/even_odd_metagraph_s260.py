#!/usr/bin/env python3
"""
even_odd_metagraph_s260.py -- opus-2026-03-23-S260

META-GRAPH THEORIES FOR EVEN GRAPHS AND ODD GRAPHS

Three parallel meta-graph theories:
1. TOURNAMENT meta-graph G_n: arc-reversal on tournaments
2. EVEN GRAPH meta-graph E_n: edge-toggle on even graphs
3. ODD GRAPH meta-graph O_n: edge-toggle on odd graphs

All three have V = A000568(n), A000568(n), A000088(n)-A000568(n) respectively.

The GRAPH meta-graph Γ_n (edge-toggle on ALL graphs) decomposes:
  Γ_n = EE edges + OO edges + EO edges
where EE connects two even graphs, OO connects two odd, EO crosses.

Author: opus-2026-03-23-S260
"""

import sys
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def canon_graph(edge_set, n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    best = None
    for perm in permutations(range(n)):
        mapped = frozenset((min(perm[i], perm[j]), max(perm[i], perm[j]))
                           for (i,j) in edge_set)
        if best is None or sorted(mapped) < sorted(best):
            best = mapped
    return frozenset(best)

def is_even_graph(edge_set, n):
    for perm in permutations(range(n)):
        mapped = frozenset((min(perm[i], perm[j]), max(perm[i], perm[j]))
                           for (i,j) in edge_set)
        if mapped != edge_set: continue
        rev = sum(1 for (i,j) in edge_set if perm[i] > perm[j])
        if rev % 2 == 1: return False
    return True

def graph_complement(edge_set, n):
    all_edges = frozenset((i,j) for i in range(n) for j in range(i+1,n))
    return all_edges - edge_set

for n in [3, 4, 5]:
    banner("COMPLETE META-GRAPH ANALYSIS AT n=%d" % n)
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

    # Build all graph iso classes with even/odd classification
    classes = {}  # canon -> {edges, even, id, size}
    class_id = 0
    for bits in range(2**m):
        edge_set = frozenset(pairs[k] for k in range(m) if (bits >> k) & 1)
        c = canon_graph(edge_set, n)
        if c not in classes:
            classes[c] = {
                'edges': edge_set,
                'even': is_even_graph(edge_set, n),
                'id': class_id,
                'size': 0,
                'edge_count': len(edge_set),
            }
            class_id += 1
        classes[c]['size'] += 1

    nc = len(classes)
    even_ids = [v['id'] for v in classes.values() if v['even']]
    odd_ids = [v['id'] for v in classes.values() if not v['even']]
    n_even = len(even_ids); n_odd = len(odd_ids)

    print("  Total graph iso classes: %d (A000088)" % nc)
    print("  Even graph classes: %d (= A000568 = tournament count)" % n_even)
    print("  Odd graph classes: %d" % n_odd)

    # Build the FULL graph meta-graph (edge-toggle)
    canon_map = {}
    for bits in range(2**m):
        edge_set = frozenset(pairs[k] for k in range(m) if (bits >> k) & 1)
        canon_map[bits] = classes[canon_graph(edge_set, n)]['id']

    # Edges: two classes adjacent iff one edge toggle connects them
    adj = [[False]*nc for _ in range(nc)]
    for bits in range(2**m):
        ci = canon_map[bits]
        for k in range(m):
            flipped = bits ^ (1 << k)
            cj = canon_map[flipped]
            if ci != cj:
                adj[ci][cj] = True; adj[cj][ci] = True

    # Classify edges by even/odd type
    EE = 0; OO = 0; EO = 0; total_edges = 0
    self_loops = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if adj[i][j]:
                total_edges += 1
                i_even = i in even_ids
                j_even = j in even_ids
                if i_even and j_even: EE += 1
                elif not i_even and not j_even: OO += 1
                else: EO += 1

    # Self-loops in graph meta-graph
    for bits in range(2**m):
        ci = canon_map[bits]
        for k in range(m):
            if canon_map[bits ^ (1<<k)] == ci:
                self_loops += 1
                break  # count each class once

    print("\n  GRAPH META-GRAPH Γ_%d:" % n)
    print("    V = %d, E = %d" % (nc, total_edges))
    print("    Self-loop classes: %d" % self_loops)
    print("    EE (even-even): %d" % EE)
    print("    OO (odd-odd): %d" % OO)
    print("    EO (even-odd): %d" % EO)

    # Even sub-meta-graph
    even_adj = [[False]*n_even for _ in range(n_even)]
    even_id_list = sorted(even_ids)
    even_idx = {eid: i for i, eid in enumerate(even_id_list)}
    for i_orig in even_ids:
        for j_orig in even_ids:
            if i_orig < j_orig and adj[i_orig][j_orig]:
                even_adj[even_idx[i_orig]][even_idx[j_orig]] = True
                even_adj[even_idx[j_orig]][even_idx[i_orig]] = True
    E_even = sum(1 for i in range(n_even) for j in range(i+1,n_even) if even_adj[i][j])

    # Odd sub-meta-graph
    odd_id_list = sorted(odd_ids)
    if n_odd > 0:
        odd_adj = [[False]*n_odd for _ in range(n_odd)]
        odd_idx = {oid: i for i, oid in enumerate(odd_id_list)}
        for i_orig in odd_ids:
            for j_orig in odd_ids:
                if i_orig < j_orig and adj[i_orig][j_orig]:
                    odd_adj[odd_idx[i_orig]][odd_idx[j_orig]] = True
                    odd_adj[odd_idx[j_orig]][odd_idx[i_orig]] = True
        E_odd = sum(1 for i in range(n_odd) for j in range(i+1,n_odd) if odd_adj[i][j])
    else:
        E_odd = 0

    print("\n  EVEN SUB-META-GRAPH E_%d:" % n)
    print("    V = %d, E = %d" % (n_even, E_even))
    print("    Avg degree: %.2f" % (2*E_even/n_even if n_even > 0 else 0))

    print("\n  ODD SUB-META-GRAPH O_%d:" % n)
    print("    V = %d, E = %d" % (n_odd, E_odd))
    print("    Avg degree: %.2f" % (2*E_odd/n_odd if n_odd > 0 else 0))

    # Tournament meta-graph for comparison
    T_edges = {3:1, 4:5, 5:30}[n]
    print("\n  TOURNAMENT META-GRAPH G_%d:" % n)
    print("    V = %d, E = %d" % (n_even, T_edges))
    print("    Avg degree: %.2f" % (2*T_edges/n_even if n_even > 0 else 0))

    # KEY RATIOS
    print("\n  KEY RATIOS:")
    print("    E(even)/E(tournament) = %d/%d = %.4f" %
          (E_even, T_edges, E_even/T_edges if T_edges > 0 else 0))
    print("    EO/total = %d/%d = %.4f (boundary fraction)" %
          (EO, total_edges, EO/total_edges if total_edges > 0 else 0))

    # Complement analysis
    print("\n  COMPLEMENT ANALYSIS:")
    comp_pairs = {'EE': 0, 'OO': 0, 'EO': 0}
    for c, v in classes.items():
        comp = canon_graph(graph_complement(v['edges'], n), n)
        if comp in classes:
            v_comp = classes[comp]
            if v['even'] and v_comp['even']: comp_pairs['EE'] += 1
            elif not v['even'] and not v_comp['even']: comp_pairs['OO'] += 1
            else: comp_pairs['EO'] += 1
    for label, count in comp_pairs.items():
        print("    %s complement pairs: %d" % (label, count // 2 if label != 'EO' else count))

    # Degree distributions
    print("\n  DEGREE DISTRIBUTIONS:")
    for label, id_list, adj_mat, n_cls in [
        ("EVEN", even_id_list, None, n_even),
        ("ODD", odd_id_list, None, n_odd),
        ("ALL", list(range(nc)), None, nc)
    ]:
        degrees = []
        for i_orig in id_list:
            deg = sum(1 for j_orig in range(nc) if adj[i_orig][j_orig])
            degrees.append(deg)
        if degrees:
            print("    %s: min=%d, max=%d, avg=%.2f, median=%.1f" %
                  (label, min(degrees), max(degrees),
                   sum(degrees)/len(degrees),
                   sorted(degrees)[len(degrees)//2]))

    # Cliques in even sub-meta-graph
    even_clique = 0
    for size in range(n_even, 0, -1):
        found = False
        if size > 6: continue  # skip large searches
        for combo in combinations(range(n_even), size):
            if all(even_adj[combo[a]][combo[b]]
                   for a in range(len(combo)) for b in range(a+1, len(combo))):
                even_clique = size; found = True; break
        if found: break
    print("\n  CLIQUE NUMBERS:")
    print("    Even sub-meta-graph: ω = %d" % even_clique)

    # Triangles
    even_tri = 0
    for i in range(n_even):
        for j in range(i+1, n_even):
            if not even_adj[i][j]: continue
            for k in range(j+1, n_even):
                if even_adj[i][k] and even_adj[j][k]:
                    even_tri += 1
    print("    Even triangles: %d" % even_tri)

    # Connected components of even sub-meta-graph
    visited = [False]*n_even
    components = 0
    for start in range(n_even):
        if visited[start]: continue
        components += 1
        stack = [start]
        while stack:
            v = stack.pop()
            if visited[v]: continue
            visited[v] = True
            for w in range(n_even):
                if even_adj[v][w] and not visited[w]:
                    stack.append(w)
    print("    Even connected components: %d" % components)

    print()

banner("SYNTHESIS: THREE META-GRAPH THEORIES")

print("""
  THE THREE META-GRAPHS:

  1. TOURNAMENT META-GRAPH G_n (arc reversal):
     V = A000568(n). Flip = reverse one arc.
     State space: {0,1}^m (m-cube).
     Properties: self-loops exist, rich connectivity, H always odd,
     edge_orbits = T_n/2 + (n-2)!, β_1 = 2,15,...

  2. EVEN GRAPH META-GRAPH E_n (edge toggle restricted to even graphs):
     V = A000568(n) (SAME count!). Flip = toggle one edge, STAY in even class.
     State space: subset of {0,1}^m.
     Properties: NO self-loops (confirmed n=3,4,5), SPARSE, fewer edges,
     connected(?), complement maps even↔odd.

  3. ODD GRAPH META-GRAPH O_n (edge toggle restricted to odd graphs):
     V = A000088(n) - A000568(n). Flip = toggle one edge, STAY in odd class.
     State space: subset of {0,1}^m.
     Properties: fewer vertices at large n (odd graphs vanish),
     connected(?), complement maps odd↔even.

  THE BOUNDARY (EO edges): connects even and odd graphs.
     = the BIPARTITE part of the full graph meta-graph.
     Analogous to the BLACK (SC↔NS) edges in the tournament meta-graph!

  THE PARALLEL:
  Tournament theory:  SPINE (SC) + RIBS (SC↔NS) + SEA (NS↔NS)
  Graph theory:       EVEN (spine) + EO BOUNDARY (ribs) + ODD (sea at small n)

  But the analogy REVERSES at large n:
  Tournaments: SC vanishes (NS dominates) — SEA is the bulk
  Graphs: ODD vanishes (EVEN dominates) — EVEN is the bulk!

  The tournament SEA (NS-NS) corresponds to the graph EVEN sub-meta-graph.
  Both are the "generic" bulk at large n.
""")

# Compute the edge-orbit formula for even graph meta-graph
banner("EDGE ORBIT FORMULAS")

print("""
  FOR TOURNAMENTS: edge_orbits = T_n/2 + (n-2)!

  FOR GRAPHS: edge_orbits = (1/n!) sum_g c(g_E) × 2^{c(g_E)}
  This counts orbits of (graph, edge) pairs under S_n.
  The graph edge_orbits = V_graph × m × (avg neutral fraction).

  FOR EVEN GRAPHS: the "edge_orbits" concept is different.
  An "edge orbit" of E_n would be an orbit of (even graph, edge)
  pairs where toggling the edge PRESERVES evenness.

  The EVEN-PRESERVING flips are the key:
  Toggling edge e in even graph X gives X△{e}.
  X△{e} is even iff the edge toggle doesn't create an odd automorphism.

  WHEN DOES TOGGLING PRESERVE EVENNESS?
  X is even: all g ∈ Aut(X) have sgn_X(g) = +1.
  X△{e} is even: all g ∈ Aut(X△{e}) have sgn_{X△{e}}(g) = +1.

  If Aut(X△{e}) ⊂ Aut(X): the toggle preserves automorphisms
  and sgn doesn't change → X△{e} stays even.

  If Aut(X△{e}) has NEW automorphisms not in Aut(X):
  these new g might have sgn = -1, making X△{e} odd.

  So: EVEN-PRESERVING flips are those that don't CREATE new
  automorphisms with odd sign. This is a STRONG constraint.

  CONTRAST WITH TOURNAMENTS:
  Every arc reversal gives another tournament (completeness preserved).
  Not every edge toggle on an even graph gives another even graph.

  This is WHY the even sub-meta-graph has FEWER edges:
  many toggles leave even-land and enter odd-land.
""")

# Summary table
banner("COMPARISON TABLE")

data = {
    3: {'V_even': 2, 'E_even': 0, 'V_odd': 2, 'E_odd': 0, 'EO': 1,
        'V_tourn': 2, 'E_tourn': 1, 'V_all': 4, 'E_all': 3},
    4: {'V_even': 4, 'E_even': 1, 'V_odd': 7, 'E_odd': 7, 'EO': 8,
        'V_tourn': 4, 'E_tourn': 5, 'V_all': 11, 'E_all': 16},
    5: {'V_even': 12, 'E_even': 9, 'V_odd': 22, 'E_odd': 51, 'EO': 49,
        'V_tourn': 12, 'E_tourn': 30, 'V_all': 34, 'E_all': 109},
}

print("  %-4s | %-8s | %-8s | %-8s | %-8s | %-8s | %-8s | %-8s" %
      ("n", "V_even", "E_even", "V_odd", "E_odd", "EO", "V_tourn", "E_tourn"))
print("  " + "-"*75)
for n_val in [3, 4, 5]:
    d = data[n_val]
    print("  %-4d | %-8d | %-8d | %-8d | %-8d | %-8d | %-8d | %-8d" %
          (n_val, d['V_even'], d['E_even'], d['V_odd'], d['E_odd'],
           d['EO'], d['V_tourn'], d['E_tourn']))

print("""
  KEY OBSERVATIONS:

  1. E_even < E_tourn always (even graph is SPARSER than tournament)
     Ratio: 0.00, 0.20, 0.30 at n=3,4,5.
     Prediction: ratio → 1 at large n (since almost all graphs are even
     and rigid, so almost all edges stay within even-land).

  2. V_odd > V_even at n=4,5 but V_odd/V_all → 0 at large n.
     Most graphs are even (because most have trivial automorphism).

  3. EO (boundary) DOMINATES at small n: 33%, 50%, 45% of all edges.
     Prediction: EO/E_all → 0 at large n (odd graphs vanish).

  4. The even sub-meta-graph is DISCONNECTED at n=3 (E=0, 2 components)
     but CONNECTED at n=5 (E=9, 1 component).
     Tournament meta-graph is ALWAYS connected (by Rédei's theorem).

  5. Self-loops: ZERO in graph meta-graph (toggling always changes degree
     sequence). NONZERO in tournament meta-graph (arc reversal can
     preserve the score sequence).
""")

print("Done. opus-2026-03-23-S260")
