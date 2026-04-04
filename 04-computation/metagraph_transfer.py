"""
metagraph_transfer.py — opus-2026-04-04-S7

THE TRANSFER MATRIX REPRESENTATION.

Key insight: the multilinear polynomial H(t) can be computed via a
TRANSFER MATRIX that processes tiles strip-by-strip along the staircase.

Strip k (anti-diagonal r+c=k) has (k-1) tiles at n.
The transfer matrix M_k maps states at strip k-1 to states at strip k.
H(t) = trace of the product M_1 · M_2 · ... · M_{n-2}.

This gives a COMPACT representation of the metagraph as a product of
small transfer matrices.

For the METAGRAPH (quotient by S_n), the transfer matrix acts on
isomorphism-class states, not individual tile states.
"""
import numpy as np
from collections import defaultdict

def make_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

def tiling_to_adj(n, tiles, bits):
    T = np.zeros((n, n), dtype=np.int8)
    for k in range(n-1):
        T[k+1][k] = 1
    for i, (x, y) in enumerate(tiles):
        a, b = x-1, y-1
        if bits & (1 << i):
            T[b][a] = 1
        else:
            T[a][b] = 1
    return T

def count_hp(T, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    full = (1 << n) - 1
    total = 0
    for mask in range(1, 1 << n):
        for v in range(n):
            cnt = dp.get((mask, v), 0)
            if cnt == 0: continue
            if mask == full: total += cnt; continue
            for u in range(n):
                if mask & (1 << u): continue
                if T[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + cnt
    return total

def strip_decomposition(n):
    """Decompose tiles into strips (anti-diagonals of the staircase).

    Tile (x,y) has strip index = x - 1 (in 1-indexed).
    Strip k (for k=2,...,n-1) contains tiles (k+1, y) for y=1,...,k-1.
    Wait, tile (x,y) with x-y≥2. Strip = x-1 in 0-indexed = x in 1-indexed minus 1.

    Actually: tile (x,y) in 1-indexed has row r=x-y-1, col c=y.
    Anti-diagonal: r+c = x-y-1+y = x-1. So strip index = x-1.

    Strip x-1 contains tiles (x, y) for y=1,...,x-2. Count = x-2.
    Strips: x-1=2,...,n-1, i.e., strips 2 through n-1.
    """
    tiles = make_tiles(n)
    strips = defaultdict(list)
    for i, (x, y) in enumerate(tiles):
        strip = x - 1  # 0-indexed strip
        strips[strip].append((i, (x, y)))

    print(f"Strip decomposition for n={n}:")
    for s in sorted(strips.keys()):
        tile_list = [(x,y) for _, (x,y) in strips[s]]
        print(f"  Strip {s}: {len(strips[s])} tiles — {tile_list}")

    return strips

def compute_strip_transfer(n):
    """Compute the CONDITIONAL H given partial tiling up to each strip.

    Process strips left to right (strip 2, 3, ..., n-1).
    At each strip, the "state" is the partial tournament on the vertices
    seen so far, and the "transfer" extends it to include the new vertex.

    For the METAGRAPH representation, the state would be the isomorphism
    class of the partial tournament.
    """
    tiles = make_tiles(n)
    m = len(tiles)
    strips = strip_decomposition(n)

    # The strips add vertices in order: strip k adds vertex k+1 (1-indexed)
    # Strip 2: tiles involving vertex 3 = {(3,1)} — vertex 3 connected to 1
    # Strip 3: tiles involving vertex 4 = {(4,1),(4,2)} — vertex 4 connected to 1,2
    # Strip k: tiles involving vertex k+1 = {(k+1,1),...,(k+1,k-1)} — vertex k+1 connected to 1,...,k-1

    # Wait, strip x-1 has tiles (x,y) for y=1..x-2. These connect vertex x to
    # vertices y (for y < x-1, non-consecutive). So adding strip x-1 adds the
    # arcs from vertex x to all vertices 1,...,x-2 (non-consecutive pairs).

    # But vertex x already has arcs to x-1 (base path). So after strip x-1,
    # vertex x is fully connected to all vertices 1,...,x-1.

    # The PARTIAL tournament after processing strips 2,...,k includes all arcs
    # among vertices 1,...,k+1.

    # This is a VERTEX-ADDITION process:
    # Start with vertices {1,2} (base path 2→1, one arc, no tiles)
    # Add vertex 3: strip 2 adds tile (3,1) connecting 3 to 1
    # Add vertex 4: strip 3 adds tiles (4,1),(4,2) connecting 4 to 1,2
    # ...
    # Add vertex n: strip n-1 adds tiles (n,1),...,(n,n-2)

    # At each step, the NEW vertex k+1 has:
    # - Base path arc to k (fixed: (k+1)→k)
    # - Tile arcs to 1,...,k-1 (binary choices: k-1 tiles)

    # The "state" after adding vertex k+1 is the tournament on {1,...,k+1}.
    # The transfer matrix maps states at k to states at k+1.

    print(f"\nVERTEX-ADDITION PROCESS for n={n}:")

    # Compute at each step k: how many distinct tournaments, their H values
    for k in range(2, n+1):
        sub_tiles = [(x,y) for x,y in tiles if x <= k and y >= 1 and x-y >= 2]
        m_k = len(sub_tiles)
        # Number of tournaments on vertices {1,...,k} with base path k→k-1→...→1
        num_tours = 1 << m_k

        # Compute H for each
        h_values = defaultdict(int)
        for bits in range(num_tours):
            T = np.zeros((k, k), dtype=np.int8)
            for j in range(k-1):
                T[j+1][j] = 1
            for idx, (x, y) in enumerate(sub_tiles):
                a, b = x-1, y-1
                if bits & (1 << idx):
                    T[b][a] = 1
                else:
                    T[a][b] = 1
            h = count_hp(T, k)
            h_values[h] += 1

        print(f"  k={k}: {m_k} tiles, {num_tours} tournaments, "
              f"{len(h_values)} distinct H, H range [{min(h_values)}, {max(h_values)}]")

    # THE KEY: compute the TRANSFER MATRIX from step k to k+1.
    # At step k: state = tournament on {1,...,k} (encoded by tiling bits for tiles ≤ k)
    # At step k+1: state = tournament on {1,...,k+1} (adds new tile bits for strip k)
    #
    # The transfer matrix T_k has:
    #   T_k[new_state, old_state] = 1 if new_state extends old_state
    #   (or 0 otherwise)
    #
    # But the number of states grows exponentially. For the METAGRAPH
    # representation, we want to work with ISOMORPHISM CLASSES instead.

    # For now, let me just show the transfer sizes.
    print(f"\nTransfer matrix dimensions:")
    for k in range(2, n):
        sub_tiles_k = [(x,y) for x,y in tiles if x <= k and x-y >= 2]
        sub_tiles_k1 = [(x,y) for x,y in tiles if x <= k+1 and x-y >= 2]
        new_tiles = len(sub_tiles_k1) - len(sub_tiles_k)
        print(f"  Step {k}→{k+1}: state dim 2^{len(sub_tiles_k)} = {1<<len(sub_tiles_k)}, "
              f"new tiles = {new_tiles}, "
              f"transfer dim = 2^{len(sub_tiles_k1)} × 2^{len(sub_tiles_k)}")

def vertex_addition_metagraph(n):
    """Build the metagraph using vertex-addition transfer matrices.

    At each step k→k+1, the transfer maps (partial iso classes) to (new iso classes).
    This gives a RECURSIVE construction of the metagraph.
    """
    tiles = make_tiles(n)

    print(f"\n{'='*60}")
    print(f"VERTEX-ADDITION METAGRAPH n={n}")
    print(f"{'='*60}")

    # At each k, compute isomorphism classes of k-tournaments with base path
    from itertools import permutations

    prev_classes = None
    for k in range(3, n+1):
        sub_tiles = [(x,y) for x,y in tiles if x <= k and x-y >= 2]
        m_k = len(sub_tiles)

        # Canonical form for k-tournaments
        classes = defaultdict(list)
        for bits in range(1 << m_k):
            T = np.zeros((k, k), dtype=np.int8)
            for j in range(k-1):
                T[j+1][j] = 1
            for idx, (x, y) in enumerate(sub_tiles):
                a, b = x-1, y-1
                if bits & (1 << idx):
                    T[b][a] = 1
                else:
                    T[a][b] = 1

            # Canonical form (brute force, OK for k ≤ 7)
            best = None
            for perm in permutations(range(k)):
                enc = tuple(T[perm[i]][perm[j]] for i in range(k) for j in range(i+1,k))
                if best is None or enc < best:
                    best = enc
            classes[best].append(bits)

        h_per_class = {}
        for canon, members in classes.items():
            bits0 = members[0]
            T = np.zeros((k, k), dtype=np.int8)
            for j in range(k-1):
                T[j+1][j] = 1
            for idx, (x, y) in enumerate(sub_tiles):
                a, b = x-1, y-1
                if bits0 & (1 << idx):
                    T[b][a] = 1
                else:
                    T[a][b] = 1
            h_per_class[canon] = count_hp(T, k)

        nc = len(classes)
        h_vals = sorted(set(h_per_class.values()))

        print(f"\n  k={k}: {nc} isomorphism classes, H values: {h_vals}")

        # If we have previous classes, compute the TRANSITION
        if prev_classes is not None:
            # How many k-classes does each (k-1)-class split into?
            prev_sub_tiles = [(x,y) for x,y in tiles if x <= k-1 and x-y >= 2]
            m_prev = len(prev_sub_tiles)

            # For each (k-1)-class, find which k-classes its members extend to
            # Each (k-1)-tiling extends by choosing the new strip tiles
            new_strip = [(x,y) for x,y in tiles if x == k and x-y >= 2]
            n_new = len(new_strip)

            expansion = defaultdict(set)
            for prev_canon, prev_members in prev_classes.items():
                # For each member, extend with all strip choices
                for prev_bits in prev_members[:1]:  # just one representative
                    for ext in range(1 << n_new):
                        # Build extended tiling
                        T = np.zeros((k, k), dtype=np.int8)
                        for j in range(k-1):
                            T[j+1][j] = 1
                        for idx, (x, y) in enumerate(prev_sub_tiles):
                            a, b = x-1, y-1
                            if prev_bits & (1 << idx):
                                T[b][a] = 1
                            else:
                                T[a][b] = 1
                        for idx, (x, y) in enumerate(new_strip):
                            a, b = x-1, y-1
                            if ext & (1 << idx):
                                T[b][a] = 1
                            else:
                                T[a][b] = 1

                        best = None
                        for perm in permutations(range(k)):
                            enc = tuple(T[perm[i]][perm[j]] for i in range(k) for j in range(i+1,k))
                            if best is None or enc < best:
                                best = enc
                        expansion[prev_canon].add(best)

            print(f"  Transition {k-1}→{k}:")
            for prev_canon in sorted(expansion.keys(), key=lambda c: prev_h[c]):
                targets = expansion[prev_canon]
                target_h = sorted(set(h_per_class[t] for t in targets))
                print(f"    H={prev_h[prev_canon]} → {len(targets)} classes, "
                      f"H values: {target_h}")

        prev_classes = classes
        prev_h = h_per_class

if __name__ == "__main__":
    strip_decomposition(5)
    strip_decomposition(7)
    compute_strip_transfer(6)
    vertex_addition_metagraph(6)
