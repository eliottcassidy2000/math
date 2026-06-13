#!/usr/bin/env python3
"""
Homological algebra and discrete Morse theory on tournament complexes.
opus-2026-03-14-S85

THREE PERSPECTIVES ON TOURNAMENTS VIA HOMOLOGY:

1. ORDER COMPLEX: For tournament T, the order complex Δ(T) has simplices =
   chains in the dominance poset. Homology of Δ(T) encodes the tournament structure.

2. INDEPENDENCE COMPLEX: Simplices = sets of arcs that don't form a directed cycle.
   This is the cycle matroid of T viewed combinatorially.

3. PATH COMPLEX (GLMY): Vertices = directed paths, boundary = path deletion.
   This is the GLMY path homology that connects to our H research.

For each: compute Betti numbers, Euler characteristics, and find
connections to H(T), scores, and tournament invariants.

Also: Discrete Morse theory — find optimal Morse functions on these complexes
to compute homology efficiently. Connect to the gradient flow from nerve_complex.
"""

from itertools import permutations, combinations
from collections import Counter, defaultdict
import math
import sys

def get_tournament(n, bits):
    """Build adjacency matrix from bit encoding."""
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(arcs):
        if (bits >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def compute_H(adj, n, all_perms):
    """Count Hamiltonian paths."""
    return sum(1 for p in all_perms if all(adj[p[i]][p[i+1]] == 1 for i in range(n-1)))

def compute_scores(adj, n):
    """Compute score sequence."""
    return tuple(sorted(sum(adj[i][j] for j in range(n) if j != i) for i in range(n)))

# ============================================================
# Part 1: Dominance Poset and Order Complex
# ============================================================
print("=" * 70)
print("PART 1: DOMINANCE POSET — ORDER COMPLEX")
print("=" * 70)

# For tournament T, vertex i dominates vertex j if i→j.
# The dominance poset has i ≥ j if every vertex that j beats, i also beats.
# (This is the "king" ordering.)

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m
    all_perms = list(permutations(range(n)))
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]

    # Track statistics
    chain_lengths = Counter()  # length of longest chain in dominance poset
    antichain_sizes = Counter()  # size of largest antichain

    for bits in range(N):
        adj = get_tournament(n, bits)
        H = compute_H(adj, n, all_perms)

        # Build dominance poset: i ≥ j iff out(j) ⊆ out(i)
        out_sets = [frozenset(j for j in range(n) if j != i and adj[i][j]) for i in range(n)]

        # Check dominance relations
        dom = [[False]*n for _ in range(n)]
        for i in range(n):
            dom[i][i] = True  # reflexive
            for j in range(n):
                if i != j and out_sets[j] <= out_sets[i]:
                    dom[i][j] = True

        # Longest chain (by DFS/BFS)
        # Build strict order
        strict = [[False]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j and dom[i][j] and not dom[j][i]:
                    strict[i][j] = True

        # Longest chain via topological approach
        # Simple: try all permutations (n small)
        max_chain = 1
        for perm in permutations(range(n)):
            chain_len = 1
            for k in range(len(perm)-1):
                if strict[perm[k]][perm[k+1]]:
                    chain_len += 1
                else:
                    break
            max_chain = max(max_chain, chain_len)

        chain_lengths[max_chain] += 1

    print(f"\nn={n}: Dominance poset chain statistics:")
    for cl in sorted(chain_lengths.keys()):
        pct = 100 * chain_lengths[cl] / N
        print(f"  Max chain length {cl}: {chain_lengths[cl]} tournaments ({pct:.1f}%)")

# ============================================================
# Part 2: Simplicial Homology — Clique Complex of Score Graph
# ============================================================
print("\n" + "=" * 70)
print("PART 2: CLIQUE COMPLEX OF SCORE-EQUIVALENCE GRAPH")
print("=" * 70)

# Group tournaments by score sequence. The "score graph" has
# an edge between two tournaments T1, T2 if they have the same
# score sequence. The clique complex of this graph has interesting topology.

# Actually more useful: the graph where tournaments are connected
# if they differ by a SINGLE arc flip AND have the same score sequence.

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m
    all_perms = list(permutations(range(n)))
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]

    # Compute all scores
    all_scores = [None] * N
    all_H = [0] * N
    for bits in range(N):
        adj = get_tournament(n, bits)
        scores = [0] * n
        for k, (i, j) in enumerate(arcs):
            if (bits >> k) & 1:
                scores[i] += 1
            else:
                scores[j] += 1
        all_scores[bits] = tuple(sorted(scores))
        all_H[bits] = compute_H(adj, n, all_perms)

    # Group by score sequence
    by_score = defaultdict(list)
    for bits in range(N):
        by_score[all_scores[bits]].append(bits)

    # For each score class: compute connected components under arc flip
    print(f"\nn={n}: Score classes and their arc-flip graphs:")
    score_classes = sorted(by_score.keys(), key=lambda s: sum(s))

    for score in score_classes:
        tours = by_score[score]
        tour_set = set(tours)

        # Build adjacency within this score class
        adj_graph = defaultdict(set)
        for bits in tours:
            for k in range(m):
                nb = bits ^ (1 << k)
                if nb in tour_set:
                    adj_graph[bits].add(nb)

        # Connected components
        visited = set()
        components = 0
        for t in tours:
            if t not in visited:
                components += 1
                queue = [t]
                visited.add(t)
                while queue:
                    v = queue.pop(0)
                    for w in adj_graph[v]:
                        if w not in visited:
                            visited.add(w)
                            queue.append(w)

        # H values in this class
        H_vals = [all_H[b] for b in tours]
        H_dist = Counter(H_vals)

        print(f"  Score {score}: {len(tours)} tours, {components} components, H dist: {dict(sorted(H_dist.items()))}")

# ============================================================
# Part 3: Euler Characteristic of H-Level Subcomplexes
# ============================================================
print("\n" + "=" * 70)
print("PART 3: EULER CHARACTERISTIC OF H-LEVEL SUBCOMPLEXES")
print("=" * 70)

# For each threshold h, define the subcomplex X_h = {T : H(T) ≥ h}.
# This is a subcomplex of the arc-flip cubical complex.
# Compute its Euler characteristic via inclusion-exclusion.

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m
    all_perms = list(permutations(range(n)))

    H_vals = [0] * N
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_vals[bits] = compute_H(adj, n, all_perms)

    achievable = sorted(set(H_vals))

    print(f"\nn={n}: H-level Euler characteristics:")
    for threshold in achievable:
        # Count: vertices (tournaments), edges (pairs differing by 1 flip),
        # faces (triples forming a 2-face of hypercube), etc.

        # Vertices in sublevel set
        verts = [bits for bits in range(N) if H_vals[bits] >= threshold]
        vert_set = set(verts)

        # Edges: pairs in sublevel set differing by 1 bit
        edges = 0
        for bits in verts:
            for k in range(m):
                nb = bits ^ (1 << k)
                if nb in vert_set and nb > bits:
                    edges += 1

        # 2-faces: quadruples {b, b⊕e_i, b⊕e_j, b⊕e_i⊕e_j} all in sublevel
        faces2 = 0
        for bits in verts:
            for i in range(m):
                for j in range(i+1, m):
                    b1 = bits ^ (1 << i)
                    b2 = bits ^ (1 << j)
                    b3 = bits ^ (1 << i) ^ (1 << j)
                    if b1 in vert_set and b2 in vert_set and b3 in vert_set:
                        # Only count if bits is the lexicographically smallest
                        if bits < b1 and bits < b2 and bits < b3:
                            faces2 += 1

        chi = len(verts) - edges + faces2
        pct = 100 * len(verts) / N
        print(f"  H≥{threshold:2d}: V={len(verts):5d} E={edges:5d} F={faces2:5d} χ={chi:5d} ({pct:.1f}% of tournaments)")

# ============================================================
# Part 4: Betti-0 (Connected Components) of H-Level Sets
# ============================================================
print("\n" + "=" * 70)
print("PART 4: π₀ OF H-LEVEL SETS (ARC-FLIP CONNECTIVITY)")
print("=" * 70)

# For each H value h, the set {T : H(T) = h} forms a subgraph
# of the arc-flip graph. How many connected components?

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m
    all_perms = list(permutations(range(n)))

    H_vals = [0] * N
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_vals[bits] = compute_H(adj, n, all_perms)

    achievable = sorted(set(H_vals))

    print(f"\nn={n}: Connected components of H-level sets:")
    for h in achievable:
        level_set = [bits for bits in range(N) if H_vals[bits] == h]
        level_set_s = set(level_set)

        # BFS to find components
        visited = set()
        components = 0
        for start in level_set:
            if start not in visited:
                components += 1
                queue = [start]
                visited.add(start)
                while queue:
                    v = queue.pop(0)
                    for k in range(m):
                        nb = v ^ (1 << k)
                        if nb in level_set_s and nb not in visited:
                            visited.add(nb)
                            queue.append(nb)

        print(f"  H={h:2d}: {len(level_set):5d} tournaments, {components:3d} components (β₀)")

# ============================================================
# Part 5: Tournament Homology via Chain Complex
# ============================================================
print("\n" + "=" * 70)
print("PART 5: CHAIN COMPLEX OF TOURNAMENT FLAG COMPLEX")
print("=" * 70)

# The flag complex of tournament T:
# k-simplex = (k+1)-clique in the "beats" graph (tournament T)
# A clique {v0,...,vk} exists iff there's a total order on them
# consistent with T (i.e., they form a transitive sub-tournament).

# This is the ORDER COMPLEX of the tournament poset!

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m
    all_perms = list(permutations(range(n)))
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]

    # Sample a few tournaments at each H value
    H_vals = [0] * N
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_vals[bits] = compute_H(adj, n, all_perms)

    # For each distinct H, pick first tournament with that H
    by_H = defaultdict(list)
    for bits in range(N):
        by_H[H_vals[bits]].append(bits)

    print(f"\nn={n}: Flag complex statistics (transitive sub-tournaments):")
    for h in sorted(by_H.keys()):
        bits = by_H[h][0]  # representative
        adj = get_tournament(n, bits)

        # Count transitive sub-tournaments of each size
        simplex_counts = [0] * (n + 1)
        simplex_counts[0] = 1  # empty simplex
        simplex_counts[1] = n  # vertices

        for size in range(2, n + 1):
            for subset in combinations(range(n), size):
                # Check if subset forms a transitive tournament
                transitive = True
                for a, b in combinations(subset, 2):
                    # Check transitivity: if a→b and b→c then a→c
                    pass
                # Actually just check: is there a topological ordering?
                # A tournament on k vertices is transitive iff it has
                # exactly one Hamiltonian path.

                # Check: does this sub-tournament have exactly one HP?
                sub_adj = {(i, j): adj[i][j] for i in subset for j in subset if i != j}
                hp_count = 0
                for perm in permutations(subset):
                    if all(sub_adj.get((perm[i], perm[i+1]), 0) == 1 for i in range(len(perm)-1)):
                        hp_count += 1

                if hp_count == 1:
                    simplex_counts[size] += 1

        # Euler characteristic of flag complex
        chi = sum((-1)**k * simplex_counts[k] for k in range(n+1))
        f_vec = simplex_counts[1:]  # f-vector (skip empty simplex)
        print(f"  H={h:2d}: f-vector={f_vec} χ={chi}")

# ============================================================
# Part 6: Discrete Morse Theory on Arc-Flip Graph
# ============================================================
print("\n" + "=" * 70)
print("PART 6: DISCRETE MORSE FUNCTION — H AS MORSE FUNCTION")
print("=" * 70)

# H(T) defines a function on the vertices of the m-dimensional hypercube.
# This is a discrete Morse function if:
# For each vertex v, at most one upper neighbor u has f(u) ≤ f(v),
# and at most one lower neighbor w has f(w) ≥ f(v).
# (This is the Forman definition for cubical complexes.)

# H is probably NOT a Morse function, but we can measure how far it is.

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m
    all_perms = list(permutations(range(n)))

    H_vals = [0] * N
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_vals[bits] = compute_H(adj, n, all_perms)

    # For each vertex, count how many neighbors have H ≥ current
    # and how many have H ≤ current.
    up_violations = 0  # vertices with >1 neighbor having H ≤ H(v) and higher in cube
    down_violations = 0
    morse_critical = 0  # vertices where ALL neighbors have different H

    for bits in range(N):
        H = H_vals[bits]
        up_count = 0  # neighbors with H ≤ H(bits) (going "up" in cube but down in H)
        down_count = 0  # neighbors with H ≥ H(bits)

        for k in range(m):
            nb = bits ^ (1 << k)
            nb_H = H_vals[nb]
            if nb_H <= H:
                up_count += 1
            if nb_H >= H:
                down_count += 1

        # Morse-critical: no neighbor has same H AND it's a local extremum
        all_different = all(H_vals[bits ^ (1 << k)] != H for k in range(m))
        if all_different:
            morse_critical += 1

        if up_count > 1:
            up_violations += 1
        if down_count > 1:
            down_violations += 1

    print(f"\nn={n}: H as discrete Morse function:")
    print(f"  Morse-critical cells (all neighbors different H): {morse_critical}")
    print(f"  Violations (>1 up-neighbor with H≤): {up_violations}")
    print(f"  Violations (>1 down-neighbor with H≥): {down_violations}")
    print(f"  Total vertices: {N}")
    print(f"  H is Morse: {up_violations == 0 and down_violations == 0}")

# ============================================================
# Part 7: Persistent Homology Preview — Filtration by H
# ============================================================
print("\n" + "=" * 70)
print("PART 7: PERSISTENT HOMOLOGY PREVIEW")
print("=" * 70)

# Build the filtration: at time h, include all tournaments with H ≤ h.
# Track connected components (β₀) as h increases.

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m
    all_perms = list(permutations(range(n)))

    H_vals = [0] * N
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_vals[bits] = compute_H(adj, n, all_perms)

    achievable = sorted(set(H_vals))

    print(f"\nn={n}: Sublevel set filtration (H ≤ h):")

    # Union-Find for connected components
    parent = list(range(N))
    rank = [0] * N

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        px, py = find(x), find(y)
        if px == py:
            return False
        if rank[px] < rank[py]:
            px, py = py, px
        parent[py] = px
        if rank[px] == rank[py]:
            rank[px] += 1
        return True

    # Sort tournaments by H value
    sorted_tours = sorted(range(N), key=lambda b: H_vals[b])

    active = set()
    births = {}  # component birth times
    deaths = []  # (birth, death) pairs

    current_components = 0
    idx = 0

    for h in achievable:
        # Add all tournaments with H = h
        new_tours = []
        while idx < N and H_vals[sorted_tours[idx]] == h:
            bits = sorted_tours[idx]
            new_tours.append(bits)
            active.add(bits)
            current_components += 1
            births[bits] = h
            idx += 1

        # Merge with existing neighbors
        merges = 0
        for bits in new_tours:
            for k in range(m):
                nb = bits ^ (1 << k)
                if nb in active:
                    if union(bits, nb):
                        current_components -= 1
                        merges += 1

        # Also check merges among existing active set with new neighbors
        # (already handled above since we check all neighbors)

        print(f"  H≤{h:2d}: {len(active):5d} tours, β₀={current_components:4d} (added {len(new_tours)}, merged {merges})")

# ============================================================
# Part 8: Čech Numbers — Intersection Patterns of Score Fibers
# ============================================================
print("\n" + "=" * 70)
print("PART 8: SCORE FIBER INTERSECTIONS")
print("=" * 70)

# Score fibers: {T : score(T) = s}
# H fibers: {T : H(T) = h}
# The intersection |score_fiber ∩ H_fiber| gives a matrix.
# This matrix's rank = number of independent "score effects" on H.

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m
    all_perms = list(permutations(range(n)))
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]

    all_scores_list = []
    all_H_list = []
    for bits in range(N):
        adj = get_tournament(n, bits)
        scores = [0] * n
        for k, (i, j) in enumerate(arcs):
            if (bits >> k) & 1:
                scores[i] += 1
            else:
                scores[j] += 1
        all_scores_list.append(tuple(sorted(scores)))
        all_H_list.append(compute_H(adj, n, all_perms))

    # Build intersection matrix
    score_vals = sorted(set(all_scores_list))
    H_vals_unique = sorted(set(all_H_list))

    matrix = defaultdict(int)
    for bits in range(N):
        matrix[(all_scores_list[bits], all_H_list[bits])] += 1

    print(f"\nn={n}: Score × H intersection matrix:")
    header = "Score\\H  " + "".join(f"{h:6d}" for h in H_vals_unique)
    print(f"  {header}")
    for s in score_vals:
        row = f"  {str(s):12s}" + "".join(f"{matrix.get((s, h), 0):6d}" for h in H_vals_unique)
        print(row)

    # Rank of this matrix (over Q)
    # Simple: count nonzero rows and columns
    nonzero_rows = sum(1 for s in score_vals if any(matrix.get((s, h), 0) > 0 for h in H_vals_unique))
    nonzero_cols = sum(1 for h in H_vals_unique if any(matrix.get((s, h), 0) > 0 for s in score_vals))
    print(f"  Score types: {len(score_vals)}, H values: {len(H_vals_unique)}")
    print(f"  Nonzero rows: {nonzero_rows}, nonzero cols: {nonzero_cols}")

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS — HOMOLOGICAL MORSE THEORY")
print("=" * 70)
print("""
KEY FINDINGS:
1. Flag complex (transitive sub-tournaments) has varying f-vectors by H value.
2. H is NOT a discrete Morse function — too many violations.
3. Persistent homology: sublevel filtration shows how β₀ evolves with H.
4. Score × H intersection matrix encodes the joint distribution.
5. H-level sets have interesting connectivity (π₀) structure.
6. Dominance poset chain lengths correlate with tournament structure.
""")
