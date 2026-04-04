"""
recursive_n_n2.py — opus-2026-04-04-S4

THE n → n+2 RECURSIVE STRUCTURE OF THE MULTILINEAR POLYNOMIAL.

Key insight: A tournament on n+2 vertices with base path fixed decomposes as:
  - Overlap: δ_{n-2} = tiles of a tournament on n inner vertices {2,...,n+1}
  - Bottom wiring: tiles (x, 1) connecting vertex 1 to inner vertices
  - Top wiring: tiles (n+2, y) connecting vertex n+2 to inner vertices
  - Apex: tile (n+2, 1) connecting the two endpoints

CLAIM: H_{n+2}(t_inner, 0_boundary) = H_n(t_inner)

When boundary tiles = 0 (forward/transitive):
  - Vertex n+2 beats all others (source)
  - Vertex 1 loses to all others (sink)
  - Every HP goes: (n+2) → [HP of inner tournament] → 1
  - So the HP count = HP count of inner tournament on {2,...,n+1}

This gives the n → n+2 recursion:
  H_{n+2}(t) = H_n(t_inner) + Σ_{S touching boundary} c_S · Π t_i
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

def tiling_to_tournament(n, tiles, bits):
    T = np.zeros((n, n), dtype=int)
    for k in range(n-1):
        T[k+1][k] = 1  # base path k+1 → k
    for i, (x, y) in enumerate(tiles):
        a, b = x-1, y-1  # 0-indexed
        if bits & (1 << i):
            T[b][a] = 1  # backward
        else:
            T[a][b] = 1  # forward
    return T

def count_hp(T, n):
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u): continue
                if T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])

def verify_recursion(n_inner, n_outer=None):
    """
    Verify: H_{n_outer}(t_inner, 0_boundary) = H_{n_inner}(t_inner)

    n_outer = n_inner + 2 by default.
    Inner vertices of n_outer are {2,...,n_outer-1} (1-indexed) = {1,...,n_outer-2} (0-indexed).
    These correspond to an n_inner-vertex tournament.
    """
    if n_outer is None:
        n_outer = n_inner + 2

    tiles_inner = make_tiles(n_inner)
    tiles_outer = make_tiles(n_outer)
    m_inner = len(tiles_inner)
    m_outer = len(tiles_outer)

    print(f"\n{'='*60}")
    print(f"n → n+2 RECURSION: {n_inner} → {n_outer}")
    print(f"Inner tiles: {m_inner}, Outer tiles: {m_outer}")
    print(f"Boundary tiles: {m_outer - m_inner}")
    print(f"{'='*60}")

    # Identify which outer tiles are "inner" (overlap with n_inner's tiles)
    # Inner vertices in n_outer (1-indexed): {2, ..., n_outer-1}
    # These have the same relative structure as vertices {1, ..., n_inner} in n_inner

    # Mapping: inner tile (x, y) in n_inner corresponds to outer tile (x+1, y+1) in n_outer
    # Because: in n_inner, vertices are 1..n_inner, tiles (x,y) with x-y≥2
    # In n_outer, inner vertices are 2..n_outer-1, tiles (x+1,y+1) where x,y from n_inner

    inner_to_outer = {}  # inner tile index → outer tile index
    outer_tile_idx = {t: i for i, t in enumerate(tiles_outer)}

    for i_inner, (x, y) in enumerate(tiles_inner):
        outer_tile = (x + 1, y + 1)  # shift by 1 for the n_outer labeling
        if outer_tile in outer_tile_idx:
            inner_to_outer[i_inner] = outer_tile_idx[outer_tile]

    print(f"\nTile mapping (inner → outer):")
    for i_inner in range(min(5, m_inner)):
        i_outer = inner_to_outer[i_inner]
        print(f"  inner tile {i_inner} {tiles_inner[i_inner]} → "
              f"outer tile {i_outer} {tiles_outer[i_outer]}")

    # Identify boundary tiles (those in outer but not in inner)
    inner_outer_set = set(inner_to_outer.values())
    boundary_indices = [i for i in range(m_outer) if i not in inner_outer_set]
    print(f"\nBoundary tiles ({len(boundary_indices)}):")
    for i in boundary_indices:
        x, y = tiles_outer[i]
        if y == 1:
            role = "BOTTOM" if x < n_outer else "APEX"
        elif x == n_outer:
            role = "TOP"
        else:
            role = "???"
        print(f"  tile {i} {tiles_outer[i]} — {role}")

    # VERIFY: for each inner tiling, H_outer(inner_bits, 0) == H_inner(inner_bits)
    print(f"\nVerifying H_{n_outer}(inner, 0) = H_{n_inner}(inner)...")
    mismatches = 0
    for inner_bits in range(1 << m_inner):
        # Compute H_inner
        T_inner = tiling_to_tournament(n_inner, tiles_inner, inner_bits)
        h_inner = count_hp(T_inner, n_inner)

        # Construct outer tiling: inner bits mapped, boundary = 0
        outer_bits = 0
        for i_inner in range(m_inner):
            if inner_bits & (1 << i_inner):
                i_outer = inner_to_outer[i_inner]
                outer_bits |= (1 << i_outer)
        # boundary bits remain 0

        T_outer = tiling_to_tournament(n_outer, tiles_outer, outer_bits)
        h_outer = count_hp(T_outer, n_outer)

        if h_inner != h_outer:
            mismatches += 1
            if mismatches <= 5:
                print(f"  MISMATCH: inner_bits={inner_bits:0{m_inner}b}, "
                      f"H_{n_inner}={h_inner}, H_{n_outer}={h_outer}")

    total = 1 << m_inner
    print(f"\n  Checked {total} inner tilings: {total - mismatches} match, {mismatches} mismatch")
    if mismatches == 0:
        print(f"  *** H_{n_outer}(inner, 0) = H_{n_inner}(inner) VERIFIED! ***")

    return mismatches == 0

def analyze_boundary_contributions(n_inner, n_outer=None):
    """
    Decompose H_{n_outer} into inner + boundary contributions.

    H_{n_outer}(t) = H_{n_inner}(t_inner) + Σ_{boundary S} c_S Π t_i

    The boundary corrections come from flipping tiles involving vertices 1 or n_outer.
    """
    if n_outer is None:
        n_outer = n_inner + 2

    tiles_inner = make_tiles(n_inner)
    tiles_outer = make_tiles(n_outer)
    m_inner = len(tiles_inner)
    m_outer = len(tiles_outer)

    # Build tile mapping
    outer_tile_idx = {t: i for i, t in enumerate(tiles_outer)}
    inner_to_outer = {}
    for i, (x, y) in enumerate(tiles_inner):
        inner_to_outer[i] = outer_tile_idx[(x+1, y+1)]

    inner_outer_set = set(inner_to_outer.values())
    boundary_indices = sorted(i for i in range(m_outer) if i not in inner_outer_set)

    # Compute full H_outer for all tilings
    H_outer = np.zeros(1 << m_outer, dtype=np.int64)
    for bits in range(1 << m_outer):
        T = tiling_to_tournament(n_outer, tiles_outer, bits)
        H_outer[bits] = count_hp(T, n_outer)

    # Möbius inversion for the outer polynomial
    coeffs_outer = {}
    for S in range(1 << m_outer):
        val = 0
        T = S
        while True:
            sign = (-1) ** bin(S ^ T).count('1')
            val += sign * int(H_outer[T])
            if T == 0: break
            T = (T - 1) & S
        if val != 0:
            coeffs_outer[S] = val

    # Classify coefficients: purely inner vs touching boundary
    pure_inner = {}  # Only involves inner tiles
    boundary_touching = {}  # Involves at least one boundary tile

    for S, c in coeffs_outer.items():
        tiles_in_S = [i for i in range(m_outer) if S & (1 << i)]
        if all(i in inner_outer_set for i in tiles_in_S):
            pure_inner[S] = c
        else:
            boundary_touching[S] = c

    print(f"\n{'='*60}")
    print(f"BOUNDARY CONTRIBUTION ANALYSIS: {n_inner} → {n_outer}")
    print(f"{'='*60}")
    print(f"Total nonzero coefficients: {len(coeffs_outer)}")
    print(f"  Pure inner: {len(pure_inner)}")
    print(f"  Boundary-touching: {len(boundary_touching)}")

    # Verify: pure inner coefficients match H_{n_inner} coefficients
    H_inner_arr = np.zeros(1 << m_inner, dtype=np.int64)
    for bits in range(1 << m_inner):
        T = tiling_to_tournament(n_inner, tiles_inner, bits)
        H_inner_arr[bits] = count_hp(T, n_inner)

    coeffs_inner = {}
    for S in range(1 << m_inner):
        val = 0
        T = S
        while True:
            sign = (-1) ** bin(S ^ T).count('1')
            val += sign * int(H_inner_arr[T])
            if T == 0: break
            T = (T - 1) & S
        if val != 0:
            coeffs_inner[S] = val

    # Map inner coefficients to outer tile indices
    inner_coeffs_mapped = {}
    for S_inner, c in coeffs_inner.items():
        S_outer = 0
        for i_inner in range(m_inner):
            if S_inner & (1 << i_inner):
                S_outer |= (1 << inner_to_outer[i_inner])
        inner_coeffs_mapped[S_outer] = c

    # Compare
    inner_match = 0
    inner_mismatch = 0
    for S_outer, c_inner in inner_coeffs_mapped.items():
        c_outer = pure_inner.get(S_outer, 0)
        if c_inner == c_outer:
            inner_match += 1
        else:
            inner_mismatch += 1
            if inner_mismatch <= 3:
                ts = [tiles_outer[i] for i in range(m_outer) if S_outer & (1 << i)]
                print(f"  INNER MISMATCH: {ts}: inner={c_inner}, outer={c_outer}")

    print(f"\n  Pure inner coefficients match H_{n_inner}: "
          f"{inner_match}/{len(inner_coeffs_mapped)}")

    # Analyze boundary-touching coefficients
    print(f"\n  Boundary-touching coefficients by order:")
    by_order = defaultdict(list)
    for S, c in boundary_touching.items():
        order = bin(S).count('1')
        n_boundary = sum(1 for i in range(m_outer) if (S & (1 << i)) and i in boundary_indices)
        n_inner_in_S = sum(1 for i in range(m_outer) if (S & (1 << i)) and i in inner_outer_set)
        by_order[order].append((S, c, n_boundary, n_inner_in_S))

    for order in sorted(by_order.keys()):
        entries = by_order[order]
        bnd_counts = [e[2] for e in entries]
        inn_counts = [e[3] for e in entries]
        print(f"    Order {order}: {len(entries)} coefficients, "
              f"boundary tiles per coeff: {min(bnd_counts)}-{max(bnd_counts)}, "
              f"inner tiles: {min(inn_counts)}-{max(inn_counts)}")

    # Show the pure boundary coefficients (only boundary tiles)
    print(f"\n  Pure boundary coefficients (no inner tiles involved):")
    for S, c in sorted(boundary_touching.items(), key=lambda x: bin(x[0]).count('1')):
        tiles_in_S = [i for i in range(m_outer) if S & (1 << i)]
        if all(i in boundary_indices for i in tiles_in_S):
            ts = [tiles_outer[i] for i in tiles_in_S]
            print(f"    {ts}: c = {c}")

    return pure_inner, boundary_touching

if __name__ == "__main__":
    # n=3 → n=5
    verify_recursion(3, 5)

    # n=4 → n=6
    verify_recursion(4, 6)

    # n=5 → n=7
    verify_recursion(5, 7)

    # Boundary analysis
    analyze_boundary_contributions(3, 5)
    analyze_boundary_contributions(4, 6)
    analyze_boundary_contributions(5, 7)
