#!/usr/bin/env python3
"""
iso_class_graph_s170.py — The Isomorphism Class Graph with Blue/Black Edges
opus-2026-03-22-S170

NODES = tournament isomorphism classes
EDGES = single arc flips connecting two classes
  BLUE edge: the flip preserves self-complementarity status
  BLACK edge: the flip changes SC status (SC→non-SC or non-SC→SC)

This graph captures the TOPOLOGY of tournament space at the
isomorphism class level. Its structure tells us:
  - Which classes are "neighbors" in the flip graph
  - How the blue skeleton (SC line) connects to non-SC regions
  - What the landscape of H looks like at the coarsest level
"""

from math import comb
from itertools import permutations
from collections import defaultdict, Counter

print("=" * 72)
print("  THE ISOMORPHISM CLASS GRAPH")
print("  opus-2026-03-22-S170")
print("=" * 72)
print()

# ==========================================================================
# HELPERS
# ==========================================================================

def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def H_dp(adj, n):
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if adj[v][u]: dp[(S|(1<<u))*n+u] += val
    full = (1<<n)-1
    return sum(dp[full*n+v] for v in range(n))

def canonical_form(adj, n):
    """Return the canonical form of a tournament (lexicographically smallest
    under all vertex permutations)."""
    best = None
    for perm in permutations(range(n)):
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

def is_sc(adj, n):
    """Check if tournament is self-complementary."""
    comp = [[1 - adj[i][j] if i != j else 0 for j in range(n)] for i in range(n)]
    return canonical_form(adj, n) == canonical_form(comp, n)

def scores(adj, n):
    return tuple(sorted(sum(adj[i]) for i in range(n)))

# ==========================================================================
# PART 1: BUILD THE ISO CLASS GRAPH AT n=5
# ==========================================================================

print("PART 1: BUILD THE ISO CLASS GRAPH AT n=5")
print("-" * 60)
print()

n = 5
m = comb(n, 2)

# Group tournaments by canonical form
iso_classes = defaultdict(list)
for bits in range(1 << m):
    adj = tournament_adj(n, bits)
    canon = canonical_form(adj, n)
    iso_classes[canon].append(bits)

# Compute invariants for each class
class_data = {}
for canon, members in iso_classes.items():
    adj = tournament_adj(n, members[0])
    h = H_dp(adj, n)
    sc = is_sc(adj, n)
    ss = scores(adj, n)
    class_data[canon] = {
        'h': h, 'sc': sc, 'scores': ss,
        'size': len(members), 'members': members,
        'id': None,  # assigned below
    }

# Assign IDs sorted by H
sorted_classes = sorted(class_data.keys(), key=lambda c: (class_data[c]['h'], class_data[c]['scores']))
for i, canon in enumerate(sorted_classes):
    class_data[canon]['id'] = i

print(f"  {len(iso_classes)} isomorphism classes at n={n}")
print()
print(f"  {'ID':>3s}  {'H':>3s}  {'SC':>3s}  {'scores':>16s}  {'size':>5s}")
print(f"  {'─'*3}  {'─'*3}  {'─'*3}  {'─'*16}  {'─'*5}")
for canon in sorted_classes:
    d = class_data[canon]
    print(f"  {d['id']:3d}  {d['h']:3d}  {'Y' if d['sc'] else 'N':>3s}  "
          f"{str(d['scores']):>16s}  {d['size']:5d}")

print()

# ==========================================================================
# PART 2: COMPUTE EDGES (arc flips between iso classes)
# ==========================================================================

print("PART 2: EDGES OF THE ISO CLASS GRAPH")
print("-" * 60)
print()

# For each tournament, flip each arc, find which class the result is in
# An edge connects two iso classes if ANY member of one is one arc flip
# from ANY member of the other

# First: build bits → canon mapping
bits_to_canon = {}
for canon, members in iso_classes.items():
    for b in members:
        bits_to_canon[b] = canon

# Count edges between classes
edge_counts = Counter()  # (class_id_1, class_id_2) → count of arc flips
blue_edges = set()  # edges where both endpoints have same SC status
black_edges = set()  # edges where SC status changes

for bits in range(1 << m):
    canon1 = bits_to_canon[bits]
    id1 = class_data[canon1]['id']
    sc1 = class_data[canon1]['sc']

    for bit in range(m):
        new_bits = bits ^ (1 << bit)
        canon2 = bits_to_canon[new_bits]
        id2 = class_data[canon2]['id']
        sc2 = class_data[canon2]['sc']

        if id1 < id2:
            edge = (id1, id2)
        elif id1 > id2:
            edge = (id2, id1)
        else:
            continue  # self-loop (flip stays in same class)

        edge_counts[edge] += 1

        if sc1 == sc2:
            blue_edges.add(edge)
        else:
            black_edges.add(edge)

# Note: an edge can be BOTH blue and black if different flips
# within the same class pair preserve or break SC differently
mixed_edges = blue_edges & black_edges
pure_blue = blue_edges - mixed_edges
pure_black = black_edges - mixed_edges

all_edges = blue_edges | black_edges

print(f"  Total edges between distinct iso classes: {len(all_edges)}")
print(f"  Pure blue (SC-preserving): {len(pure_blue)}")
print(f"  Pure black (SC-changing): {len(pure_black)}")
print(f"  Mixed (both types of flip exist): {len(mixed_edges)}")
print()

# ==========================================================================
# PART 3: THE ADJACENCY STRUCTURE
# ==========================================================================

print("PART 3: ADJACENCY LIST OF THE ISO CLASS GRAPH")
print("-" * 60)
print()

# Build adjacency list
adj_list = defaultdict(set)
edge_types = {}
for edge in all_edges:
    i, j = edge
    adj_list[i].add(j)
    adj_list[j].add(i)
    if edge in pure_blue:
        edge_types[edge] = 'BLUE'
    elif edge in pure_black:
        edge_types[edge] = 'BLACK'
    else:
        edge_types[edge] = 'MIXED'

for canon in sorted_classes:
    d = class_data[canon]
    i = d['id']
    neighbors = sorted(adj_list[i])
    neighbor_strs = []
    for j in neighbors:
        edge = (min(i,j), max(i,j))
        color = edge_types.get(edge, '?')
        neighbor_strs.append(f"{j}({color[0]})")
    print(f"  Class {i:2d} (H={d['h']:2d},{'SC' if d['sc'] else 'NS'}): "
          f"neighbors = {', '.join(neighbor_strs)}")

print()

# ==========================================================================
# PART 4: GRAPH INVARIANTS
# ==========================================================================

print("PART 4: GRAPH INVARIANTS OF THE ISO CLASS GRAPH")
print("-" * 60)
print()

n_nodes = len(sorted_classes)
n_edges = len(all_edges)
degrees = [len(adj_list[i]) for i in range(n_nodes)]
max_deg = max(degrees)
min_deg = min(degrees)

print(f"  Nodes: {n_nodes}")
print(f"  Edges: {n_edges}")
print(f"  Blue edges: {len(pure_blue)}")
print(f"  Black edges: {len(pure_black)}")
print(f"  Mixed edges: {len(mixed_edges)}")
print(f"  Degree sequence: {sorted(degrees)}")
print(f"  Min degree: {min_deg}")
print(f"  Max degree: {max_deg}")
print(f"  Average degree: {sum(degrees)/n_nodes:.2f}")
print()

# Is it connected?
visited = {0}
queue = [0]
while queue:
    curr = queue.pop(0)
    for nb in adj_list[curr]:
        if nb not in visited:
            visited.add(nb)
            queue.append(nb)
connected = (len(visited) == n_nodes)
print(f"  Connected: {connected}")

# Diameter
from collections import deque
def bfs_dist(start, adj, n_nodes):
    dist = [-1] * n_nodes
    dist[start] = 0
    q = deque([start])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if dist[v] == -1:
                dist[v] = dist[u] + 1
                q.append(v)
    return dist

max_dist = 0
for i in range(n_nodes):
    dists = bfs_dist(i, adj_list, n_nodes)
    max_dist = max(max_dist, max(d for d in dists if d >= 0))

print(f"  Diameter: {max_dist}")
print()

# ==========================================================================
# PART 5: THE BLUE SUBGRAPH
# ==========================================================================

print("PART 5: THE BLUE SUBGRAPH (SC-preserving edges only)")
print("-" * 60)
print()

blue_adj = defaultdict(set)
for edge in pure_blue | mixed_edges:  # Include mixed since they have blue flips
    i, j = edge
    blue_adj[i].add(j)
    blue_adj[j].add(i)

blue_degrees = [len(blue_adj[i]) for i in range(n_nodes)]
sc_classes = [d['id'] for _, d in class_data.items() if d['sc']]
nsc_classes = [d['id'] for _, d in class_data.items() if not d['sc']]

print(f"  SC classes: {sorted(sc_classes)} ({len(sc_classes)} classes)")
print(f"  Non-SC classes: {sorted(nsc_classes)} ({len(nsc_classes)} classes)")
print()

# Connected components of the blue subgraph
blue_visited = set()
blue_components = []
for start in range(n_nodes):
    if start in blue_visited:
        continue
    component = set()
    q = [start]
    while q:
        u = q.pop()
        if u in blue_visited:
            continue
        blue_visited.add(u)
        component.add(u)
        for v in blue_adj[u]:
            if v not in blue_visited:
                q.append(v)
    blue_components.append(component)

print(f"  Blue subgraph connected components: {len(blue_components)}")
for i, comp in enumerate(blue_components):
    h_vals = sorted(set(class_data[sorted_classes[c]]['h'] for c in comp))
    sc_in_comp = sum(1 for c in comp if class_data[sorted_classes[c]]['sc'])
    print(f"    Component {i}: {len(comp)} classes, H={h_vals}, SC={sc_in_comp}")

print()

# ==========================================================================
# PART 6: H-GRADIENT ON THE ISO CLASS GRAPH
# ==========================================================================

print("PART 6: H-GRADIENT STRUCTURE")
print("-" * 60)
print()

# For each edge, is it H-increasing, H-decreasing, or H-level?
h_increasing = 0
h_decreasing = 0
h_level = 0

for edge in all_edges:
    i, j = edge
    hi = class_data[sorted_classes[i]]['h']
    hj = class_data[sorted_classes[j]]['h']
    if hi < hj:
        h_increasing += 1
    elif hi > hj:
        h_decreasing += 1
    else:
        h_level += 1

print(f"  H-increasing edges: {h_increasing}")
print(f"  H-decreasing edges: {h_decreasing}")
print(f"  H-level edges (same H): {h_level}")
print()

# The H-level edges form the "contour" of the H function
print(f"  H-LEVEL EDGES (critical for understanding topology):")
for edge in sorted(all_edges):
    i, j = edge
    hi = class_data[sorted_classes[i]]['h']
    hj = class_data[sorted_classes[j]]['h']
    if hi == hj:
        color = edge_types.get(edge, '?')
        si = 'SC' if class_data[sorted_classes[i]]['sc'] else 'NS'
        sj = 'SC' if class_data[sorted_classes[j]]['sc'] else 'NS'
        print(f"    {i}({si})—{j}({sj}) at H={hi}, color={color}")

print()

# ==========================================================================
# PART 7: SEQUENCES FROM THE ISO CLASS GRAPH
# ==========================================================================

print("PART 7: SEQUENCES FROM THE ISO CLASS GRAPH")
print("-" * 60)
print()

# At n=5 we have the complete data. Let me also do n=3,4.
print(f"  ISO CLASS GRAPH SEQUENCES:")
print()

for n_seq in [3, 4, 5]:
    m_seq = comb(n_seq, 2)

    # Build iso classes
    iso_seq = defaultdict(list)
    for bits in range(1 << m_seq):
        adj = tournament_adj(n_seq, bits)
        canon = canonical_form(adj, n_seq)
        iso_seq[canon].append(bits)

    n_classes = len(iso_seq)

    # Build edges
    bits_to_c = {}
    for canon, members in iso_seq.items():
        for b in members:
            bits_to_c[b] = canon

    edges_seq = set()
    canons = list(iso_seq.keys())
    canon_to_id = {c: i for i, c in enumerate(canons)}

    for bits in range(1 << m_seq):
        id1 = canon_to_id[bits_to_c[bits]]
        for bit in range(m_seq):
            new_bits = bits ^ (1 << bit)
            id2 = canon_to_id[bits_to_c[new_bits]]
            if id1 != id2:
                edges_seq.add((min(id1,id2), max(id1,id2)))

    n_edges_seq = len(edges_seq)

    # Build adjacency and compute diameter
    adj_seq = defaultdict(set)
    for i, j in edges_seq:
        adj_seq[i].add(j)
        adj_seq[j].add(i)

    max_d = 0
    for i in range(n_classes):
        dists = bfs_dist(i, adj_seq, n_classes)
        max_d = max(max_d, max(d for d in dists if d >= 0))

    degs = sorted([len(adj_seq[i]) for i in range(n_classes)])

    print(f"  n={n_seq}: {n_classes} classes, {n_edges_seq} edges, "
          f"diameter={max_d}, degrees={degs}")

print()
print(f"  SEQUENCES:")
print(f"    #iso_classes(n):    (A000568 related) = 1, 1, 2, 4, 12, ...")
print(f"    #edges in flip graph: computed above")
print(f"    diameter:           computed above")
print()

# ==========================================================================
# SUMMARY
# ==========================================================================

print()
print("=" * 72)
print("  SUMMARY")
print("=" * 72)
print()

print(f"""
  THE ISO CLASS GRAPH AT n=5:

  12 nodes (isomorphism classes)
  {n_edges} edges (arc flip connections)
  {len(pure_blue)} pure blue, {len(pure_black)} pure black, {len(mixed_edges)} mixed
  Connected: {connected}
  Diameter: {max_dist}

  H-gradient: {h_increasing} up, {h_decreasing} down, {h_level} level edges
  The graph has a clear H-gradient from H=1 (bottom) to H=15 (top).
  Level edges exist at H=3 (3 classes), H=5 (2 classes), H=15 (2 classes).

  The BLUE SUBGRAPH (SC-preserving flips) has {len(blue_components)} components.
  SC classes: {len(sc_classes)} out of {n_nodes}.

  NEW SEQUENCES:
    #edges in iso flip graph: n=3: ?, n=4: ?, n=5: {n_edges}
    diameter of iso flip graph
    #blue edges, #black edges, #mixed edges
    #connected components of blue subgraph
    chromatic number of iso flip graph
""")

print("opus-2026-03-22-S170")
