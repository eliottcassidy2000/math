"""
boundary_polynomial_pattern.py — opus-2026-04-04-S4

Study how the boundary correction Δ_n evolves with n.

The boundary tiles at n+2 are:
  Bottom: (3,1), (4,1), ..., (n+2,1)  [n tiles total, skips 2,...,n+1]
  Top:    (n+2,2), (n+2,3), ..., (n+2,n)  [n-1 tiles, skips n,...,2]
  Apex:   (n+2,1)  [1 tile, skip n+1]

Wait, in the n+2 tournament with vertices 1,...,n+2:
  Bottom: (x,1) for x=3,...,n+2 → n tiles, skips 2,...,n+1
  Top: (n+2,y) for y=1,...,n → n tiles, skips n+1,...,2
  But (n+2,1) is both bottom and top! It's the APEX.

Actually: bottom = (x,1) for x=3,...,n+1 → n-1 tiles (not including apex)
         top = (n+2,y) for y=2,...,n → n-1 tiles (not including apex)
         apex = (n+2,1) → 1 tile
         Total: (n-1)+(n-1)+1 = 2n-1 ✓

KEY QUESTION: Does the PURE BOUNDARY polynomial (no inner tiles) have a
pattern that's independent of the inner structure?

SECOND QUESTION: Do the bottom-top interactions follow the disjoint pair
formula from THM-287?

THIRD QUESTION: The boundary polynomial at n+2 includes tiles with the SAME
skip structure as an (n-1)-tournament's tiles. Is there a deeper connection?
"""
import numpy as np
from collections import defaultdict
from recursive_n_n2 import make_tiles, tiling_to_tournament, count_hp

def compute_all_H(n):
    tiles = make_tiles(n)
    m = len(tiles)
    H = np.zeros(1 << m, dtype=np.int64)
    for bits in range(1 << m):
        T = tiling_to_tournament(n, tiles, bits)
        H[bits] = count_hp(T, n)
    return H, tiles, m

def mobius_inversion(H, m):
    coeffs = {}
    for S in range(1 << m):
        val = 0
        T = S
        while True:
            sign = (-1) ** bin(S ^ T).count('1')
            val += sign * int(H[T])
            if T == 0: break
            T = (T - 1) & S
        if val != 0:
            coeffs[S] = val
    return coeffs

def extract_pure_boundary(n_outer):
    """Extract the pure boundary multilinear polynomial at n_outer."""
    n_inner = n_outer - 2
    tiles_outer = make_tiles(n_outer)
    tiles_inner = make_tiles(n_inner)
    m_outer = len(tiles_outer)
    m_inner = len(tiles_inner)

    # Identify boundary tiles
    outer_tile_idx = {t: i for i, t in enumerate(tiles_outer)}
    inner_outer_set = set()
    for i, (x, y) in enumerate(tiles_inner):
        if (x+1, y+1) in outer_tile_idx:
            inner_outer_set.add(outer_tile_idx[(x+1, y+1)])

    boundary_idx = sorted(i for i in range(m_outer) if i not in inner_outer_set)
    n_bnd = len(boundary_idx)

    # The pure boundary polynomial: H_{n_outer} restricted to inner=0
    # This means: set all inner tiles to 0, vary boundary tiles
    # H_{n_outer}(0_inner, t_bnd) = Σ_S c_S Π t_i (summing over subsets of boundary tiles)

    # Compute by brute force: for each boundary tiling, compute H with inner=0
    H_bnd = np.zeros(1 << n_bnd, dtype=np.int64)
    for bnd_bits in range(1 << n_bnd):
        # Construct full outer tiling: inner=0, boundary from bnd_bits
        outer_bits = 0
        for j, b_idx in enumerate(boundary_idx):
            if bnd_bits & (1 << j):
                outer_bits |= (1 << b_idx)
        T = tiling_to_tournament(n_outer, tiles_outer, outer_bits)
        H_bnd[bnd_bits] = count_hp(T, n_outer)

    # Möbius inversion on boundary tiles
    bnd_coeffs = {}
    for S in range(1 << n_bnd):
        val = 0
        T = S
        while True:
            sign = (-1) ** bin(S ^ T).count('1')
            val += sign * int(H_bnd[T])
            if T == 0: break
            T = (T - 1) & S
        if val != 0:
            bnd_coeffs[S] = val

    bnd_tiles = [tiles_outer[i] for i in boundary_idx]

    return H_bnd, bnd_coeffs, bnd_tiles, boundary_idx

def compare_boundary_across_n():
    """Compare pure boundary polynomials across different n."""
    print("PURE BOUNDARY POLYNOMIAL EVOLUTION")
    print("="*60)

    for n_outer in [5, 6, 7]:
        n_inner = n_outer - 2
        H_bnd, bnd_coeffs, bnd_tiles, bnd_idx = extract_pure_boundary(n_outer)
        n_bnd = len(bnd_tiles)

        print(f"\nn_outer={n_outer}, boundary tiles={n_bnd}")
        print(f"  Tiles: {bnd_tiles}")
        print(f"  H range with inner=0: [{H_bnd.min()}, {H_bnd.max()}], mean={H_bnd.mean():.1f}")
        print(f"  H(0,0,...,0) = {H_bnd[0]} (= H_{n_inner}(0) = 1)")

        # Group boundary tiles
        bottom = [(i, t) for i, t in enumerate(bnd_tiles) if t[1] == 1 and t[0] < n_outer]
        top = [(i, t) for i, t in enumerate(bnd_tiles) if t[0] == n_outer and t[1] > 1]
        apex = [(i, t) for i, t in enumerate(bnd_tiles) if t[0] == n_outer and t[1] == 1]

        print(f"  Bottom ({len(bottom)}): {[t for _,t in bottom]}")
        print(f"  Top ({len(top)}): {[t for _,t in top]}")
        print(f"  Apex ({len(apex)}): {[t for _,t in apex]}")

        # Linear coefficients
        print(f"\n  Linear coefficients (pure boundary):")
        for j in range(n_bnd):
            c = bnd_coeffs.get(1 << j, 0)
            tile = bnd_tiles[j]
            skip = tile[0] - tile[1]
            expected = 2**(skip-1)
            role = 'bottom' if tile[1]==1 and tile[0]<n_outer else ('top' if tile[0]==n_outer else 'apex')
            match = '✓' if c == expected else f'✗ ({expected})'
            print(f"    {tile} skip={skip} [{role}]: c={c} {match}")

        # KEY: The pure boundary polynomial H(0, t_bnd) — what IS this?
        # It's H of a tournament where inner vertices form a transitive tournament
        # and boundary vertices connect to them via the boundary tile settings.
        # This is the "transitive inner + boundary wiring" tournament.

        # For this tournament:
        # - Inner vertices {2,...,n-1} form a transitive tournament (all descending)
        # - Vertex 1 connects to inner via bottom wiring tiles
        # - Vertex n+2 connects to inner via top wiring tiles
        # - Vertex 1 and n+2 connect via apex

        # The HP count for this tournament depends ONLY on the boundary tile settings!

        # OBSERVATION: The boundary polynomial H(0, t_bnd) looks like a SMALLER
        # tournament's H polynomial, possibly with a different vertex structure.

        # Check: does H(0, t_bnd) equal anything known?
        # H(0, 0...0) = 1 (all transitive)
        # H(0, all boundary = 1) = ?

        all_bnd = (1 << n_bnd) - 1
        print(f"\n  H(0_inner, 1_boundary) = {H_bnd[all_bnd]}")
        print(f"  H(0_inner, 0_boundary) = {H_bnd[0]}")

        # The boundary polynomial degree
        max_order_bnd = max(bin(S).count('1') for S in bnd_coeffs.keys()) if bnd_coeffs else 0
        by_order = defaultdict(int)
        for S, c in bnd_coeffs.items():
            by_order[bin(S).count('1')] += 1
        print(f"\n  Boundary polynomial degree: {max_order_bnd}")
        print(f"  Coefficients per order: {dict(sorted(by_order.items()))}")

def analyze_bottom_top_structure(n_outer=7):
    """
    The bottom wiring tiles have skips 2,3,...,n (n-1 tiles).
    The top wiring tiles have skips n,...,3,2 (n-1 tiles, reversed).
    The apex has skip n+1 (max).

    The bottom tiles form a "staircase column" with the same skip distribution
    as tiles of an (n-1)-tournament!

    Specifically: bottom tiles (x,1) for x=3,...,n+1 have skips 2,...,n.
    An (n-1)-tournament on vertices 1,...,n-1 has base path (n-1)→...→1 and tiles
    (x,y) with x-y≥2. The tiles (x,1) with x=3,...,n-1 have skips 2,...,n-2.

    But the bottom wiring has skips 2,...,n, which is LARGER. The extra tiles
    have skips n-1 and n, connecting vertex 1 to the top of the inner tournament.

    Actually, the bottom wiring tiles (x,1) for x=3,...,n+2 in the (n+2)-tournament
    have skips 2,...,n+1. The number of bottom tiles is n (excluding apex when
    n+2 is counted separately).

    Hmm, this doesn't directly correspond to a smaller tournament's tile set.
    But the SKIP SEQUENCE is 2,3,...,n+1 (for n+2-tournament boundary).

    Let me check if the boundary polynomial has any FACTORIZATION.
    """
    n_inner = n_outer - 2
    H_bnd, bnd_coeffs, bnd_tiles, bnd_idx = extract_pure_boundary(n_outer)
    n_bnd = len(bnd_tiles)

    print(f"\n{'='*60}")
    print(f"BOTTOM-TOP STRUCTURE at n_outer={n_outer}")
    print(f"{'='*60}")

    # Identify bottom and top indices within boundary
    bottom_idx = []
    top_idx = []
    apex_idx = None
    for j, t in enumerate(bnd_tiles):
        if t[1] == 1 and t[0] == n_outer:
            apex_idx = j
        elif t[1] == 1:
            bottom_idx.append(j)
        elif t[0] == n_outer:
            top_idx.append(j)

    print(f"Bottom tiles: {[bnd_tiles[j] for j in bottom_idx]}")
    print(f"Top tiles: {[bnd_tiles[j] for j in top_idx]}")
    print(f"Apex tile: {bnd_tiles[apex_idx] if apex_idx is not None else None}")

    # Check: with apex=0, does the boundary FACTORIZE as bottom × top?
    # i.e., H(0_inner, t_bottom, t_top, 0_apex) =? f(t_bottom) × g(t_top)

    n_bottom = len(bottom_idx)
    n_top = len(top_idx)

    H_bottom_only = np.zeros(1 << n_bottom, dtype=np.int64)
    H_top_only = np.zeros(1 << n_top, dtype=np.int64)

    # H with only bottom tiles varied (top=0, apex=0, inner=0)
    for b_bits in range(1 << n_bottom):
        bnd_bits = 0
        for j, b_j in enumerate(bottom_idx):
            if b_bits & (1 << j):
                bnd_bits |= (1 << b_j)
        H_bottom_only[b_bits] = H_bnd[bnd_bits]

    # H with only top tiles varied (bottom=0, apex=0, inner=0)
    for t_bits in range(1 << n_top):
        bnd_bits = 0
        for j, t_j in enumerate(top_idx):
            if t_bits & (1 << j):
                bnd_bits |= (1 << t_j)
        H_top_only[t_bits] = H_bnd[bnd_bits]

    # Check factorization: H(bottom, top, 0_apex) =? H(bottom, 0, 0) × H(0, top, 0) / H(0,0,0)
    factorization_holds = True
    for b_bits in range(1 << n_bottom):
        for t_bits in range(1 << n_top):
            # Construct combined boundary bits
            bnd_bits = 0
            for j, b_j in enumerate(bottom_idx):
                if b_bits & (1 << j):
                    bnd_bits |= (1 << b_j)
            for j, t_j in enumerate(top_idx):
                if t_bits & (1 << j):
                    bnd_bits |= (1 << t_j)
            # apex = 0

            h_combined = H_bnd[bnd_bits]
            h_bottom = H_bottom_only[b_bits]
            h_top = H_top_only[t_bits]
            h_00 = H_bnd[0]

            if h_00 > 0:
                product = h_bottom * h_top
                actual = h_combined * h_00
                if product != actual:
                    factorization_holds = False
                    if b_bits < 4 and t_bits < 4:
                        print(f"  Factorization fail: bottom={b_bits}, top={t_bits}, "
                              f"H_bt={h_combined}, H_b={h_bottom}, H_t={h_top}, "
                              f"H_bt·H_00={actual}, H_b·H_t={product}")

    if factorization_holds:
        print(f"\n  *** BOTTOM-TOP FACTORIZATION HOLDS (apex=0)! ***")
        print(f"  H(bottom, top, 0_apex) × H(0,0,0) = H(bottom, 0, 0) × H(0, top, 0)")
        print(f"  The boundary polynomial FACTORS into bottom and top when apex=0!")
    else:
        print(f"\n  Bottom-top factorization: FAILS")

    # The bottom-only polynomial
    print(f"\n  Bottom-only polynomial:")
    bottom_tiles = [bnd_tiles[j] for j in bottom_idx]
    for b_bits in range(1 << n_bottom):
        t_list = [bottom_tiles[j] for j in range(n_bottom) if b_bits & (1 << j)]
        if b_bits < 8 or b_bits == (1 << n_bottom) - 1:
            print(f"    bits={b_bits:0{n_bottom}b}: H={H_bottom_only[b_bits]} "
                  f"tiles flipped: {t_list}")

    # The top-only polynomial
    print(f"\n  Top-only polynomial:")
    top_tiles = [bnd_tiles[j] for j in top_idx]
    for t_bits in range(1 << n_top):
        t_list = [top_tiles[j] for j in range(n_top) if t_bits & (1 << j)]
        if t_bits < 8 or t_bits == (1 << n_top) - 1:
            print(f"    bits={t_bits:0{n_top}b}: H={H_top_only[t_bits]} "
                  f"tiles flipped: {t_list}")

if __name__ == "__main__":
    compare_boundary_across_n()
    analyze_bottom_top_structure(7)
