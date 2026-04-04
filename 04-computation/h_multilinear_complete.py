"""
h_multilinear_complete.py — opus-2026-04-04-S1

Complete extraction of H(T) as multilinear polynomial in tile variables.
H(t_1,...,t_m) = sum_S c_S prod_{i in S} t_i
where c_S = sum_{T subset S} (-1)^{|S|-|T|} H(1_T)

This gives the FULL structure — not just linear/quadratic but all orders.
Key questions:
  - What is the maximum degree (order of interaction)?
  - How do coefficients depend on tile geometry (skip, position)?
  - Why does degree cap at 4 for n=5,6 but jump to 6 at n=7?
  - What is the apex effect correction at n>=7?
"""
import numpy as np
from itertools import combinations
from collections import defaultdict
import sys

def make_tiles(n):
    """Return list of tiles (x,y) in the standard explorer order."""
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

def tiling_to_tournament(n, tiles, tiling_bits):
    """Convert a tiling (bit vector) to tournament adjacency matrix.
    Bit 0 = forward arc a->b (base path direction), bit 1 = backward arc b->a."""
    T = np.zeros((n, n), dtype=int)
    # Base path: n -> n-1 -> ... -> 1 (vertices labeled 1..n, using 0-indexed: n-1 -> n-2 -> ... -> 0)
    for k in range(n-1):
        T[k+1][k] = 1  # k+1 -> k (base path arc)

    # Tiles
    for i, (x, y) in enumerate(tiles):
        # Convert to 0-indexed
        a, b = x-1, y-1  # a > b+1 (non-consecutive)
        if tiling_bits & (1 << i):
            T[b][a] = 1  # backward: b -> a
        else:
            T[a][b] = 1  # forward: a -> b
    return T

def count_hamiltonian_paths(T, n):
    """Count directed Hamiltonian paths via DP bitmask."""
    # dp[mask][v] = number of paths visiting exactly the vertices in mask, ending at v
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])

def compute_all_H(n):
    """Compute H for all 2^m tilings."""
    tiles = make_tiles(n)
    m = len(tiles)
    H = np.zeros(1 << m, dtype=np.int64)
    for bits in range(1 << m):
        T = tiling_to_tournament(n, tiles, bits)
        H[bits] = count_hamiltonian_paths(T, n)
    return H, tiles, m

def mobius_inversion(H_values, m):
    """Extract multilinear coefficients via Möbius inversion on Boolean lattice.
    c_S = sum_{T subset S} (-1)^{|S|-|T|} H(1_T)
    """
    coeffs = {}
    for S in range(1 << m):
        val = 0
        # Iterate over all subsets T of S
        T = S
        while True:
            sign = (-1) ** (bin(S ^ T).count('1'))  # |S \ T|
            val += sign * H_values[T]
            if T == 0:
                break
            T = (T - 1) & S
        if val != 0:
            coeffs[S] = val
    return coeffs

def subset_to_tiles(S, tiles):
    """Convert bitmask S to list of tile positions."""
    return [tiles[i] for i in range(len(tiles)) if S & (1 << i)]

def analyze_coefficients(coeffs, tiles, m, n):
    """Analyze the multilinear coefficient spectrum."""
    print(f"\n{'='*60}")
    print(f"n={n}, m={m} tiles, 2^m={1<<m} tilings")
    print(f"{'='*60}")

    # Group by order (number of tiles in the interaction)
    by_order = defaultdict(list)
    for S, c in coeffs.items():
        order = bin(S).count('1')
        tile_set = subset_to_tiles(S, tiles)
        skips = [x-y for x,y in tile_set]
        by_order[order].append((S, c, tile_set, skips))

    print(f"\nTotal nonzero coefficients: {len(coeffs)}")
    print(f"Maximum order (degree): {max(by_order.keys()) if by_order else 0}")

    for order in sorted(by_order.keys()):
        entries = by_order[order]
        values = [c for _, c, _, _ in entries]
        print(f"\n--- Order {order}: {len(entries)} nonzero coefficients ---")
        print(f"  Values: min={min(values)}, max={max(values)}, distinct={len(set(values))}")

        if order == 0:
            print(f"  c_empty = {values[0]} (= H at all-zero tiling = H(transitive))")
        elif order == 1:
            print(f"  Linear coefficients (should be 2^(skip-1) by THM-284):")
            for S, c, tile_set, skips in sorted(entries, key=lambda x: x[2]):
                (x,y) = tile_set[0]
                skip = x - y
                expected = 2**(skip-1)
                match = "✓" if c == expected else f"✗ (expected {expected})"
                print(f"    tile ({x},{y}) skip={skip}: c={c} {match}")
        elif order == 2:
            print(f"  Quadratic coefficients:")
            # Group by relationship between tiles
            same_vertex = []
            disjoint = []
            for S, c, tile_set, skips in entries:
                (x1,y1), (x2,y2) = tile_set
                shared = set()
                if x1 == x2: shared.add(f"upper={x1}")
                if y1 == y2: shared.add(f"lower={y1}")
                if x1 == y2: shared.add(f"x1=y2={x1}")
                if y1 == x2: shared.add(f"y1=x2={y1}")
                if shared:
                    same_vertex.append((c, tile_set, skips, shared))
                else:
                    disjoint.append((c, tile_set, skips))

            if same_vertex:
                print(f"  Shared vertex ({len(same_vertex)} pairs):")
                for c, ts, sk, sh in sorted(same_vertex, key=lambda x: x[0]):
                    print(f"    {ts} skips={sk} shared={sh}: c₂={c}")
            if disjoint:
                print(f"  Disjoint ({len(disjoint)} pairs):")
                for c, ts, sk in sorted(disjoint, key=lambda x: x[0]):
                    print(f"    {ts} skips={sk}: c₂={c}")
        else:
            # Higher orders — just show summary and a few examples
            print(f"  Sample coefficients:")
            for S, c, tile_set, skips in sorted(entries, key=lambda x: abs(x[1]), reverse=True)[:5]:
                print(f"    {tile_set} skips={skips}: c={c}")

    # Special analysis: degree structure
    max_order = max(by_order.keys()) if by_order else 0
    print(f"\n{'='*60}")
    print(f"DEGREE ANALYSIS for n={n}")
    print(f"{'='*60}")
    print(f"Maximum degree: {max_order}")
    print(f"Coefficients per degree: {', '.join(f'd{k}:{len(by_order[k])}' for k in sorted(by_order.keys()))}")

    # Check: is the degree related to floor((n-1)/2)*2?
    expected_walsh_deg = 2 * ((n-1)//2)
    print(f"Expected Walsh degree (2*floor((n-1)/2)): {expected_walsh_deg}")
    print(f"Actual max degree: {max_order}")
    if max_order == expected_walsh_deg:
        print("  MATCH with Walsh degree formula!")
    elif max_order < expected_walsh_deg:
        print(f"  BELOW Walsh degree (gap = {expected_walsh_deg - max_order})")
    else:
        print(f"  ABOVE Walsh degree (excess = {max_order - expected_walsh_deg})")

    # Apex tile analysis
    apex_idx = None
    for i, (x,y) in enumerate(tiles):
        if x == n and y == 1:
            apex_idx = i
            break
    if apex_idx is not None:
        apex_coeffs = [(S, c) for S, c in coeffs.items() if S & (1 << apex_idx)]
        print(f"\nApex tile ({n},1) participates in {len(apex_coeffs)} nonzero coefficients")
        apex_by_order = defaultdict(list)
        for S, c in apex_coeffs:
            apex_by_order[bin(S).count('1')].append(c)
        for k in sorted(apex_by_order.keys()):
            vals = apex_by_order[k]
            print(f"  Order {k}: {len(vals)} coeffs, sum={sum(vals)}, "
                  f"values={sorted(set(vals))}")

    return by_order

def main():
    for n in [5, 6]:
        print(f"\n\nComputing H for all tilings at n={n}...")
        H, tiles, m = compute_all_H(n)
        print(f"H range: [{H.min()}, {H.max()}], mean={H.mean():.1f}")
        print(f"H(transitive) = H(0) = {H[0]}")
        print(f"H(anti-transitive) = H({(1<<m)-1}) = {H[(1<<m)-1]}")

        print(f"\nExtracting multilinear coefficients...")
        coeffs = mobius_inversion(H, m)
        by_order = analyze_coefficients(coeffs, tiles, m, n)

    # n=7 is bigger but doable (2^15 = 32768)
    n = 7
    print(f"\n\nComputing H for all tilings at n={n} (2^{n*(n-1)//2 - (n-1)} = {1 << (n*(n-1)//2 - (n-1))})...")
    H, tiles, m = compute_all_H(n)
    print(f"H range: [{H.min()}, {H.max()}], mean={H.mean():.1f}")

    print(f"\nExtracting multilinear coefficients (this may take a moment)...")
    coeffs = mobius_inversion(H, m)
    by_order = analyze_coefficients(coeffs, tiles, m, n)

if __name__ == "__main__":
    main()
