"""
persistent_H_landscape.py -- kind-pasteur-2026-03-14-S74
Persistent homology of the H-landscape on the tournament hypercube.

IDEA: The sublevel sets L_k = {T : H(T) <= k} form a filtration of the
tournament hypercube. As k increases, the topology changes (births and deaths
of connected components, loops, voids). The barcode/persistence diagram
encodes this topological information.

ALSO: The SUPERLEVEL sets U_k = {T : H(T) >= k} filtration, tracking
how the "peaks" merge as we lower the threshold.

CONNECTION TO H=7 GAP: The filtration jumps from L_5 to L_9 (no tournaments
with H=7). This creates a "topological gap" — a large birth-death interval.

COMPUTATIONAL APPROACH:
For n=5 (1024 tournaments on 10-dim hypercube):
1. Build the adjacency graph of the hypercube (tournaments connected by single arc flip)
2. Compute connected components of L_k for each k
3. Track births and deaths = persistence diagram
4. Compute Betti numbers of sublevel sets

Also: MORSE THEORY on the landscape
- Critical points = local maxima/minima/saddle points (already computed in S71)
- Morse inequalities: #critical points of index k >= beta_k of sublevel sets
"""

import numpy as np
from collections import Counter, defaultdict, deque
import sys

sys.stdout.reconfigure(encoding='utf-8')

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def compute_H_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def union_find_create(n):
    return list(range(n))

def union_find_find(parent, i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i

def union_find_union(parent, rank, i, j):
    ri, rj = union_find_find(parent, i), union_find_find(parent, j)
    if ri == rj:
        return False
    if rank[ri] < rank[rj]:
        ri, rj = rj, ri
    parent[rj] = ri
    if rank[ri] == rank[rj]:
        rank[ri] += 1
    return True

def main():
    print("=" * 70)
    print("PERSISTENT HOMOLOGY OF THE H-LANDSCAPE")
    print("kind-pasteur-2026-03-14-S74")
    print("=" * 70)

    for n in [4, 5]:
        m = n * (n - 1) // 2
        N = 2**m
        print(f"\n{'='*70}")
        print(f"n = {n}, m = {m}, N = {N} tournaments")
        print(f"{'='*70}")

        # Compute all H values
        H_values = np.zeros(N, dtype=int)
        for bits in range(N):
            A = bits_to_adj(bits, n)
            H_values[bits] = compute_H_dp(A, n)

        H_set = sorted(set(H_values))
        print(f"\n  H values: {H_set}")
        print(f"  H distribution: {dict(Counter(H_values))}")

        # ========================================
        # SUBLEVEL SET PERSISTENCE (H_0)
        # ========================================
        print(f"\n  --- SUBLEVEL SET PERSISTENCE ---")
        print(f"  Track connected components of L_k = {{T : H(T) <= k}}")

        # Sort tournaments by H value (ascending)
        sorted_tournaments = sorted(range(N), key=lambda b: H_values[b])

        # Union-Find for connected components
        parent = union_find_create(N)
        rank_uf = [0] * N
        active = set()
        births = []  # (birth_H, component_id)
        deaths = []  # (birth_H, death_H)

        component_birth = {}  # root -> birth_H

        prev_H = -1
        for idx, bits in enumerate(sorted_tournaments):
            H = H_values[bits]

            # Add this tournament to active set
            active.add(bits)
            component_birth[bits] = H  # born at H

            # Check neighbors (single arc flips)
            for arc_bit in range(m):
                nbr = bits ^ (1 << arc_bit)
                if nbr in active:
                    root_bits = union_find_find(parent, bits)
                    root_nbr = union_find_find(parent, nbr)
                    if root_bits != root_nbr:
                        # Merge: the YOUNGER component dies
                        birth_bits = component_birth.get(root_bits, H)
                        birth_nbr = component_birth.get(root_nbr, H)

                        if birth_bits <= birth_nbr:
                            # nbr's component is younger, it dies
                            deaths.append((birth_nbr, H))
                            union_find_union(parent, rank_uf, root_bits, root_nbr)
                            new_root = union_find_find(parent, bits)
                            component_birth[new_root] = birth_bits
                        else:
                            deaths.append((birth_bits, H))
                            union_find_union(parent, rank_uf, root_nbr, root_bits)
                            new_root = union_find_find(parent, bits)
                            component_birth[new_root] = birth_nbr

        # Components still alive
        alive_roots = set()
        for bits in active:
            alive_roots.add(union_find_find(parent, bits))
        alive_births = [component_birth[r] for r in alive_roots]

        print(f"\n  Persistence pairs (birth, death):")
        persistence = sorted(deaths, key=lambda x: x[1] - x[0], reverse=True)
        for b, d in persistence[:15]:
            lifetime = d - b
            print(f"    born at H={b}, dies at H={d}, lifetime={lifetime}")

        print(f"\n  Still alive at H=max: {len(alive_roots)} component(s) born at H={sorted(set(alive_births))}")

        # Persistence diagram statistics
        lifetimes = [d - b for b, d in deaths]
        if lifetimes:
            print(f"\n  Lifetime statistics:")
            print(f"    Max lifetime: {max(lifetimes)}")
            print(f"    Mean lifetime: {np.mean(lifetimes):.2f}")
            print(f"    Median lifetime: {np.median(lifetimes):.1f}")
            print(f"    # pairs with lifetime > 0: {sum(1 for l in lifetimes if l > 0)}")
            print(f"    # pairs with lifetime = 0: {sum(1 for l in lifetimes if l == 0)}")

        # Betti_0 at each H level
        print(f"\n  Connected components (Betti_0) at each H level:")
        parent2 = union_find_create(N)
        rank_uf2 = [0] * N
        active2 = set()
        components_at_H = {}

        for H_level in H_set:
            # Add all tournaments with H == H_level
            new_bits = [b for b in range(N) if H_values[b] == H_level]
            for bits in new_bits:
                active2.add(bits)

            # Connect to neighbors
            for bits in new_bits:
                for arc_bit in range(m):
                    nbr = bits ^ (1 << arc_bit)
                    if nbr in active2:
                        union_find_union(parent2, rank_uf2, bits, nbr)

            # Count components
            roots = set(union_find_find(parent2, b) for b in active2)
            components_at_H[H_level] = len(roots)

        for H_level in H_set:
            count = sum(1 for b in range(N) if H_values[b] == H_level)
            print(f"    L_{H_level}: {sum(1 for b in range(N) if H_values[b] <= H_level)} tournaments, "
                  f"beta_0 = {components_at_H[H_level]} components")

        # ========================================
        # THE H=7 GAP IN TOPOLOGY
        # ========================================
        if 7 not in set(H_values) and n == 5:
            print(f"\n  --- THE H=7 GAP IN PERSISTENT TOPOLOGY ---")
            # L_5 and L_9 differ by the H=9 tournaments
            # How many connected components does L_5 have?
            # Do the H=9 tournaments connect previously-separate components?
            print(f"    L_5 has {components_at_H[5]} components")
            print(f"    L_9 has {components_at_H[9]} components")
            print(f"    The H=9 tournaments ({sum(1 for b in range(N) if H_values[b] == 9)} of them) connect")
            print(f"    {components_at_H[5] - components_at_H[9]} pairs of components from L_5")
            print(f"    This is the TOPOLOGICAL SIGNATURE of the H=7 gap!")

        # ========================================
        # SUPERLEVEL SET PERSISTENCE (descending)
        # ========================================
        print(f"\n  --- SUPERLEVEL SET PERSISTENCE ---")
        print(f"  Track connected components of U_k = {{T : H(T) >= k}}")

        sorted_desc = sorted(range(N), key=lambda b: -H_values[b])

        parent3 = union_find_create(N)
        rank_uf3 = [0] * N
        active3 = set()
        component_birth3 = {}
        deaths3 = []

        for bits in sorted_desc:
            H = H_values[bits]
            active3.add(bits)
            component_birth3[bits] = H

            for arc_bit in range(m):
                nbr = bits ^ (1 << arc_bit)
                if nbr in active3:
                    root_bits = union_find_find(parent3, bits)
                    root_nbr = union_find_find(parent3, nbr)
                    if root_bits != root_nbr:
                        birth_bits = component_birth3.get(root_bits, H)
                        birth_nbr = component_birth3.get(root_nbr, H)
                        # In superlevel: higher birth = older
                        if birth_bits >= birth_nbr:
                            deaths3.append((birth_nbr, H))
                            union_find_union(parent3, rank_uf3, root_bits, root_nbr)
                            new_root = union_find_find(parent3, bits)
                            component_birth3[new_root] = birth_bits
                        else:
                            deaths3.append((birth_bits, H))
                            union_find_union(parent3, rank_uf3, root_nbr, root_bits)
                            new_root = union_find_find(parent3, bits)
                            component_birth3[new_root] = birth_nbr

        persistence3 = sorted(deaths3, key=lambda x: x[0] - x[1], reverse=True)
        print(f"\n  Superlevel persistence pairs (birth_H, death_H) [descending]:")
        for b, d in persistence3[:10]:
            lifetime = b - d
            print(f"    born at H={b}, dies at H={d}, lifetime={lifetime}")

    # ========================================
    # EULER CHARACTERISTIC OF SUBLEVEL SETS
    # ========================================
    print(f"\n{'='*70}")
    print("EULER CHARACTERISTIC OF SUBLEVEL SETS")
    print("  chi(L_k) = #vertices - #edges + #triangles - ...")
    print("  where the simplicial complex is the clique complex of the")
    print("  subgraph of the hypercube restricted to L_k")
    print(f"{'='*70}")

    n = 5
    m = n * (n - 1) // 2
    N = 2**m
    H_values_5 = np.zeros(N, dtype=int)
    for bits in range(N):
        A = bits_to_adj(bits, n)
        H_values_5[bits] = compute_H_dp(A, n)

    H_set_5 = sorted(set(H_values_5))

    for H_level in H_set_5:
        # Vertices in L_k
        verts = [b for b in range(N) if H_values_5[b] <= H_level]
        vert_set = set(verts)

        # Edges in L_k (arc-flip neighbors both in L_k)
        edges = 0
        for b in verts:
            for arc_bit in range(m):
                nbr = b ^ (1 << arc_bit)
                if nbr in vert_set and nbr > b:
                    edges += 1

        # Euler characteristic approximation (just vertices - edges)
        chi_approx = len(verts) - edges
        print(f"  L_{H_level}: {len(verts)} vertices, {edges} edges, chi_approx = {chi_approx}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
