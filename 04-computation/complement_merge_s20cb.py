#!/usr/bin/env python3
"""
complement_merge_s20cb.py -- kind-pasteur-2026-03-22-S20cb

THE COMPLEMENT-MERGED META-GRAPH G_n / Z_2.

Merge each complement pair (C, C^op) into a single node.
SC classes stay as individual nodes (they're already self-inverse).

The merged graph has:
  Vertices: SC_count + NS_count/2 = (A000568 + SC_count) / 2
  Edges: decompose into blue (SC-preserving) and black (SC-changing)

In the merged graph:
  - Blue edges between SC nodes: preserved
  - Blue edges between merged NS nodes: also preserved
  - Black edges between SC and NS: SC node <-> merged NS node
  - Black edges between different NS pairs: merged <-> merged

The user's key insight: after merging, the distinction between
SC and merged-NS nodes becomes LESS important, because both are
"self-inverse" in the complement sense.

Author: kind-pasteur-2026-03-22-S20cb
"""
import sys
import numpy as np
from math import comb, factorial
from collections import defaultdict
from itertools import permutations
sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best: best = form
    return best

def is_sc(A, n):
    A_comp = np.zeros_like(A)
    for i in range(n):
        for j in range(n):
            if i != j: A_comp[i][j] = 1 - A[i][j]
    for perm in permutations(range(n)):
        if all(A[perm[i]][perm[j]] == A_comp[i][j] for i in range(n) for j in range(n) if i != j):
            return True
    return False

print("=" * 70)
print("  THE COMPLEMENT-MERGED META-GRAPH G_n / Z_2")
print("=" * 70)

for n in [3, 4, 5, 6]:
    print(f"\n{'='*70}")
    print(f"  n = {n}")
    print(f"{'='*70}\n")

    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)

    # Build iso classes
    canon_map = defaultdict(list)
    H_map = {}
    for bits in range(2**m):
        A = np.zeros((n,n), dtype=np.int8)
        for k, (i,j) in enumerate(pairs):
            if (bits >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        H = count_hp(A, n)
        H_map[bits] = H
        cf = canonical(A, n)
        canon_map[cf].append(bits)

    classes = []
    cf_to_id = {}
    for cf, members in sorted(canon_map.items(), key=lambda x: H_map[x[1][0]]):
        cid = len(classes)
        cf_to_id[cf] = cid
        A = np.zeros((n,n), dtype=np.int8)
        for k, (i,j) in enumerate(pairs):
            if (members[0] >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        sc = is_sc(A, n) if n <= 6 else None
        classes.append({'id': cid, 'H': H_map[members[0]], 'size': len(members),
                       'sc': sc, 'members': set(members)})

    bits_to_class = {}
    for c in classes:
        for b in c['members']:
            bits_to_class[b] = c['id']

    N = len(classes)

    # Find complement pairs
    complement_map = {}  # class_id -> complement class_id
    for c in classes:
        T = list(c['members'])[0]
        # Complement: flip all bits
        T_comp = T ^ ((1 << m) - 1)
        comp_class = bits_to_class[T_comp]
        complement_map[c['id']] = comp_class

    # Identify SC classes and NS pairs
    sc_classes = [c['id'] for c in classes if c['sc']]
    ns_pairs = []
    ns_seen = set()
    for c in classes:
        if c['sc']: continue
        if c['id'] in ns_seen: continue
        comp = complement_map[c['id']]
        ns_pairs.append((c['id'], comp))
        ns_seen.add(c['id'])
        ns_seen.add(comp)

    n_merged = len(sc_classes) + len(ns_pairs)
    print(f"  Original: {N} iso classes ({len(sc_classes)} SC + {2*len(ns_pairs)} NS)")
    print(f"  Merged: {n_merged} nodes ({len(sc_classes)} SC + {len(ns_pairs)} merged NS pairs)")
    print(f"  Formula: ({N} + {len(sc_classes)}) / 2 = {(N + len(sc_classes)) / 2}")

    # Build the merged graph
    # Node IDs: SC classes keep their ID. NS pairs get a new merged ID.
    merged_ids = {}
    for sc in sc_classes:
        merged_ids[sc] = sc  # SC keeps its ID
    for (a, b) in ns_pairs:
        merged_id = min(a, b)  # use the smaller ID for the pair
        merged_ids[a] = merged_id
        merged_ids[b] = merged_id

    # Build adjacency for G_n (original)
    adj = defaultdict(set)
    edge_colors = {}  # (min_id, max_id) -> 'blue' or 'black'

    for c in classes:
        T = list(c['members'])[0]
        for k in range(m):
            nb = T ^ (1 << k)
            nb_class = bits_to_class[nb]
            if nb_class != c['id']:
                edge = (min(c['id'], nb_class), max(c['id'], nb_class))
                adj[c['id']].add(nb_class)
                # Color
                if c['sc'] == classes[nb_class]['sc']:
                    edge_colors[edge] = 'blue'
                else:
                    edge_colors[edge] = 'black'

    n_blue = sum(1 for v in edge_colors.values() if v == 'blue')
    n_black = sum(1 for v in edge_colors.values() if v == 'black')
    n_edges = len(edge_colors)

    print(f"\n  Original G_n: {n_edges} edges = {n_blue} blue + {n_black} black")

    # Build merged adjacency
    merged_adj = defaultdict(set)
    merged_edge_colors = {}

    for (a, b), color in edge_colors.items():
        ma = merged_ids[a]
        mb = merged_ids[b]
        if ma != mb:
            edge = (min(ma, mb), max(ma, mb))
            merged_adj[ma].add(mb)
            merged_adj[mb].add(ma)
            # In merged graph: what's the color?
            # The edge connects two merged nodes.
            merged_edge_colors[edge] = color

    # Self-loops in merged graph (edges within a complement pair)
    merged_self_loops = 0
    for (a, b), color in edge_colors.items():
        ma = merged_ids[a]
        mb = merged_ids[b]
        if ma == mb:
            merged_self_loops += 1

    n_merged_edges = len(merged_edge_colors)
    n_merged_blue = sum(1 for v in merged_edge_colors.values() if v == 'blue')
    n_merged_black = sum(1 for v in merged_edge_colors.values() if v == 'black')

    print(f"  Merged G_n/Z_2: {n_merged_edges} edges = {n_merged_blue} blue + {n_merged_black} black")
    print(f"  Collapsed edges (self-loops in merged): {merged_self_loops}")
    print(f"  Check: {n_merged_edges} + {merged_self_loops} = {n_merged_edges + merged_self_loops} vs original {n_edges}")

    # ================================================================
    # THE BLUE/BLACK DECOMPOSITION
    # ================================================================
    print(f"\n  BLUE/BLACK DECOMPOSITION:")
    print(f"    Original: {n_blue} blue + {n_black} black = {n_edges}")
    print(f"    Merged:   {n_merged_blue} blue + {n_merged_black} black = {n_merged_edges}")
    print(f"    Collapsed: {merged_self_loops}")

    # Degree analysis in merged graph
    merged_nodes = set(merged_ids.values())
    print(f"\n  MERGED DEGREE SEQUENCE:")
    degrees_merged = {}
    for node in sorted(merged_nodes):
        d = len(merged_adj[node])
        is_sc_node = node in sc_classes
        node_type = "SC" if is_sc_node else "NS-pair"
        # Find the H values this merged node represents
        if is_sc_node:
            h = classes[node]['H']
        else:
            h = classes[node]['H']  # both in pair have same H
        degrees_merged[node] = d
        if n <= 5:
            print(f"    node {node:>3d} ({node_type:>7s}, H={h:>3d}): degree {d}")

    # ================================================================
    # THE KEY: IS THE MERGED GRAPH SIMPLER?
    # ================================================================
    print(f"\n  SIMPLICITY COMPARISON:")
    print(f"    Original: V={N}, E={n_edges}, density={2*n_edges/(N*(N-1)):.4f}")
    print(f"    Merged:   V={n_merged}, E={n_merged_edges}, density={2*n_merged_edges/(n_merged*(n_merged-1)):.4f}" if n_merged > 1 else "")

    # Blue fraction
    print(f"    Blue fraction original: {n_blue}/{n_edges} = {n_blue/n_edges:.3f}")
    if n_merged_edges > 0:
        print(f"    Blue fraction merged:   {n_merged_blue}/{n_merged_edges} = {n_merged_blue/n_merged_edges:.3f}")

print(f"\n{'='*70}")
print(f"  SYNTHESIS: THE COMPLEMENT-MERGED PICTURE")
print(f"{'='*70}\n")

print(f"""  THE MERGED GRAPH G_n/Z_2 has:
    Vertices = (A000568(n) + SC_count(n)) / 2
    Edges = blue_edges + black_edges - collapsed_edges

  After merging complement pairs:
  1. SC nodes and merged-NS nodes are BOTH "self-inverse"
  2. Blue edges connect same-type nodes (SC-SC or merged-merged)
  3. Black edges connect different-type nodes (SC-merged or merged-merged cross)
  4. Some original edges COLLAPSE into self-loops (edges within a complement pair)

  The blue/black split is cleaner in the merged graph because
  each node is now self-inverse, removing the asymmetry between
  SC and NS classes.

  EDGE COUNT = BLUE + BLACK + COLLAPSED:
  This decomposition separates the edge count into three parts,
  each potentially having a cleaner formula than the total.
""")
