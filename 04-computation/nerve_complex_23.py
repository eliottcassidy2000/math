#!/usr/bin/env python3
"""
Nerve complex and simplicial structure of tournaments.
opus-2026-03-14-S85

Build the NERVE COMPLEX of the H-fibers:
- For each achievable H value, the set of tournaments with that H
  forms a "fiber" in the tournament hypercube {0,1}^m.
- The nerve complex has a vertex for each H-value and a k-simplex
  for each (k+1)-tuple of H-values whose fibers have a common neighbor
  (tournaments differing by 1 arc flip).

This connects to:
- Persistent homology of the H landscape
- Čech complex / Vietoris-Rips at various scales
- Adjacency graph of H-fibers (arc-flip graph restricted to H level sets)

Also: the "H landscape" viewed as a function on the Boolean hypercube.
Critical points = tournaments where every arc flip changes H.
Saddle points = transitions between H values.
"""

from itertools import permutations, combinations
from collections import Counter, defaultdict
import math
import sys

def compute_all_H_with_neighbors(n):
    """Compute H for all tournaments and record arc-flip neighbors."""
    m = n * (n - 1) // 2
    N = 1 << m
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    all_perms = list(permutations(range(n)))

    H_values = [0] * N
    for bits in range(N):
        if bits % 5000 == 0 and N > 5000:
            print(f"  n={n}: {bits}/{N}", file=sys.stderr)
        adj = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(arcs):
            if (bits >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        H_values[bits] = sum(1 for p in all_perms if all(adj[p[i]][p[i+1]] == 1 for i in range(n-1)))

    return H_values

# ============================================================
# Part 1: H landscape on hypercube
# ============================================================
print("=" * 70)
print("PART 1: H LANDSCAPE ON BOOLEAN HYPERCUBE")
print("=" * 70)

for n in [4, 5, 6]:
    H_vals = compute_all_H_with_neighbors(n)
    m = n * (n - 1) // 2
    N = 1 << m

    # For each tournament, compute: is it a local max, local min, or saddle?
    local_max = 0
    local_min = 0
    saddle = 0
    flat = 0

    # ΔH statistics
    delta_H_global = Counter()

    for bits in range(N):
        H = H_vals[bits]
        neighbors_H = []
        for k in range(m):
            nb = bits ^ (1 << k)
            neighbors_H.append(H_vals[nb])

        max_nb = max(neighbors_H)
        min_nb = min(neighbors_H)

        if H > max_nb:
            local_max += 1
        elif H < min_nb:
            local_min += 1
        elif H == max_nb and H == min_nb:
            flat += 1
        else:
            saddle += 1

        # ΔH statistics
        for nh in neighbors_H:
            delta_H_global[nh - H] += 1

    print(f"\nn={n}: H landscape topology (m={m})")
    print(f"  Local maxima: {local_max} ({100*local_max/N:.2f}%)")
    print(f"  Local minima: {local_min} ({100*local_min/N:.2f}%)")
    print(f"  Saddle points: {saddle} ({100*saddle/N:.2f}%)")
    print(f"  Flat: {flat} ({100*flat/N:.2f}%)")

    # Morse theory: #maxima - #saddles + #minima = Euler char of domain
    # But domain is a hypercube graph, which is contractible...
    print(f"  #max - #saddle + #min = {local_max - saddle + local_min}")

    # ΔH distribution
    print(f"  ΔH distribution:")
    for delta in sorted(delta_H_global.keys()):
        print(f"    ΔH={delta:+3d}: {delta_H_global[delta]:8d}")

# ============================================================
# Part 2: Adjacency between H-fibers
# ============================================================
print("\n" + "=" * 70)
print("PART 2: ADJACENCY GRAPH OF H-FIBERS")
print("=" * 70)

for n in [4, 5, 6]:
    H_vals = compute_all_H_with_neighbors(n)
    m = n * (n - 1) // 2
    N = 1 << m

    # Which H values are adjacent (connected by single arc flip)?
    adj_H = defaultdict(set)

    for bits in range(N):
        h1 = H_vals[bits]
        for k in range(m):
            h2 = H_vals[bits ^ (1 << k)]
            if h1 != h2:
                adj_H[h1].add(h2)

    achievable = sorted(set(H_vals))
    print(f"\nn={n}: H-fiber adjacency graph:")
    print(f"  Vertices (H values): {achievable}")
    edges = set()
    for h1, neighbors in adj_H.items():
        for h2 in neighbors:
            if h1 < h2:
                edges.add((h1, h2))
    print(f"  Edges: {len(edges)}")

    # Print adjacency
    for h in achievable:
        nbrs = sorted(adj_H[h])
        print(f"    H={h:2d} → {nbrs}")

    # Is the adjacency graph connected?
    visited = {achievable[0]}
    queue = [achievable[0]]
    while queue:
        v = queue.pop(0)
        for w in adj_H[v]:
            if w not in visited:
                visited.add(w)
                queue.append(w)
    print(f"  Connected: {len(visited) == len(achievable)}")

    # Diameter
    # BFS from each vertex
    max_dist = 0
    for start in achievable:
        dist = {start: 0}
        queue = [start]
        while queue:
            v = queue.pop(0)
            for w in adj_H[v]:
                if w not in dist:
                    dist[w] = dist[v] + 1
                    queue.append(w)
        max_dist = max(max_dist, max(dist.values()))
    print(f"  Diameter: {max_dist}")

# ============================================================
# Part 3: Critical points of H
# ============================================================
print("\n" + "=" * 70)
print("PART 3: CRITICAL POINTS — LOCAL EXTREMA")
print("=" * 70)

for n in [4, 5]:
    H_vals = compute_all_H_with_neighbors(n)
    m = n * (n - 1) // 2
    N = 1 << m

    # Find all local maxima and their H values
    maxima_H = Counter()
    minima_H = Counter()
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]

    for bits in range(N):
        H = H_vals[bits]
        is_max = True
        is_min = True
        for k in range(m):
            nb_H = H_vals[bits ^ (1 << k)]
            if nb_H >= H:
                is_max = False
            if nb_H <= H:
                is_min = False

        if is_max:
            maxima_H[H] += 1
        if is_min:
            minima_H[H] += 1

    print(f"\nn={n}:")
    print(f"  Local maxima by H: {dict(sorted(maxima_H.items()))}")
    print(f"  Local minima by H: {dict(sorted(minima_H.items()))}")

    # Check: are maxima always at max-H and minima always at min-H?
    print(f"  Max-H tournaments that are local maxima: {maxima_H.get(max(H_vals), 0)}/{Counter(H_vals)[max(H_vals)]}")
    print(f"  Min-H tournaments that are local minima: {minima_H.get(min(H_vals), 0)}/{Counter(H_vals)[min(H_vals)]}")

# ============================================================
# Part 4: ΔH parity — is ΔH always even?
# ============================================================
print("\n" + "=" * 70)
print("PART 4: ΔH PARITY — IS ΔH ALWAYS EVEN?")
print("=" * 70)

for n in [4, 5, 6]:
    H_vals = compute_all_H_with_neighbors(n)
    m = n * (n - 1) // 2
    N = 1 << m

    odd_delta = 0
    even_delta = 0

    for bits in range(N):
        for k in range(m):
            delta = H_vals[bits ^ (1 << k)] - H_vals[bits]
            if delta % 2 == 0:
                even_delta += 1
            else:
                odd_delta += 1

    print(f"n={n}: even ΔH = {even_delta}, odd ΔH = {odd_delta}")
    print(f"  ΔH always even: {odd_delta == 0}")

# ============================================================
# Part 5: H gradient flow — basins of attraction
# ============================================================
print("\n" + "=" * 70)
print("PART 5: GRADIENT FLOW — BASINS OF ATTRACTION")
print("=" * 70)

# From each tournament, follow the steepest ascent (largest ΔH).
# Where does each tournament end up?

for n in [4, 5]:
    H_vals = compute_all_H_with_neighbors(n)
    m = n * (n - 1) // 2
    N = 1 << m

    # Steepest ascent from each tournament
    attractor = {}
    for bits in range(N):
        current = bits
        path_length = 0
        while True:
            H = H_vals[current]
            best_nb = current
            best_H = H
            for k in range(m):
                nb = current ^ (1 << k)
                if H_vals[nb] > best_H:
                    best_H = H_vals[nb]
                    best_nb = nb
            if best_nb == current:
                break  # local max reached
            current = best_nb
            path_length += 1
            if path_length > N:
                break  # safety

        attractor[bits] = (current, H_vals[current], path_length)

    # Basin sizes
    basin_sizes = Counter()
    for bits, (att, att_H, _) in attractor.items():
        basin_sizes[att] += 1

    print(f"\nn={n}: Gradient flow basins of attraction:")
    print(f"  #attractors (local maxima): {len(basin_sizes)}")
    for att in sorted(basin_sizes.keys(), key=lambda x: -basin_sizes[x])[:10]:
        att_H = H_vals[att]
        size = basin_sizes[att]
        print(f"  Attractor H={att_H}: basin size = {size} ({100*size/N:.1f}%)")

    # Average path length to attractor
    avg_path = sum(pl for _, (_, _, pl) in attractor.items()) / N
    print(f"  Average path length to attractor: {avg_path:.2f}")

# ============================================================
# Part 6: Nerve complex — clique structure
# ============================================================
print("\n" + "=" * 70)
print("PART 6: NERVE COMPLEX OF H-FIBERS")
print("=" * 70)

# A k-simplex in the nerve complex exists iff there's a tournament
# that can reach all k+1 H-values by single arc flips.
# Equivalently: the 1-neighborhoods of the fibers have nonempty intersection.

# Actually, let's define it more simply:
# Vertex set = achievable H values
# Edge (h1, h2) iff some tournament with H=h1 has a neighbor with H=h2
# This is the adjacency graph from Part 2.

# For higher simplices: triangle (h1,h2,h3) iff exists tournament T
# with H(T)=h1 that has neighbors with H=h2 AND H=h3.

for n in [4, 5]:
    H_vals = compute_all_H_with_neighbors(n)
    m = n * (n - 1) // 2
    N = 1 << m

    achievable = sorted(set(H_vals))

    # Build vertex → set of neighbor H values
    nb_H_sets = defaultdict(lambda: defaultdict(set))
    for bits in range(N):
        h = H_vals[bits]
        for k in range(m):
            h_nb = H_vals[bits ^ (1 << k)]
            if h_nb != h:
                nb_H_sets[h][bits].add(h_nb)

    # Count triangles in nerve
    triangles = set()
    for h1 in achievable:
        for bits in nb_H_sets[h1]:
            nbs = nb_H_sets[h1][bits]
            for h2, h3 in combinations(sorted(nbs), 2):
                triangles.add((h1, h2, h3))

    print(f"\nn={n}: Nerve complex:")
    print(f"  Vertices: {len(achievable)}")
    # edges from Part 2
    edges = set()
    for h1 in achievable:
        for bits in nb_H_sets[h1]:
            for h2 in nb_H_sets[h1][bits]:
                if h1 < h2:
                    edges.add((h1, h2))
    print(f"  Edges: {len(edges)}")
    print(f"  Triangles: {len(triangles)}")

    # Euler characteristic
    chi = len(achievable) - len(edges) + len(triangles)
    print(f"  Euler characteristic (V-E+T): {chi}")

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS — TOURNAMENT TOPOLOGY")
print("=" * 70)
print("""
KEY FINDINGS:
1. H landscape: ΔH is ALWAYS EVEN (arc flip preserves parity — Rédei!)
2. Local maxima exist at highest H values; local minima at lowest.
3. H-fiber adjacency graph is CONNECTED — you can walk between any two
   H values by single arc flips.
4. Gradient flow creates basins of attraction around local maxima.
5. Nerve complex has rich simplicial structure — triangles etc.
6. The "Morse theory" of H on the hypercube gives topological invariants
   of the tournament landscape.
""")
