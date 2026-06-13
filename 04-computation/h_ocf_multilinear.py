"""
h_ocf_multilinear.py — opus-2026-04-04-S2

CREATIVE DIRECTION: Decompose H(t) via OCF into cycle-set contributions.

H(T) = I(Ω(T), 2) = Σ_{I independent} 2^|I|

where I ranges over independent sets of vertex-disjoint odd directed cycles.

As a function of tile variables:
  H(t) = Σ_{I independent} 2^|I| · Π_{C ∈ I} χ_C(t)

where χ_C(t) = 1 iff all arcs of cycle C are present in T(t).

Each χ_C(t) is a MULTILINEAR FUNCTION of the tiles involved in C:
  χ_C(t) = Π_{arcs of C that are tile arcs} (t_k or (1-t_k))

This gives a COMBINATORIAL INTERPRETATION of each multilinear coefficient:
  c_S = Σ_{I independent} 2^|I| · [coefficient of Π_{i∈S} t_i in Π_{C∈I} χ_C(t)]

The goal: compute this decomposition and verify it reproduces the
Möbius-computed coefficients exactly.
"""
import numpy as np
from itertools import combinations
from collections import defaultdict
from fractions import Fraction

def make_tiles(n):
    """Tiles in explorer order."""
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

def find_all_odd_cycles(n):
    """Find ALL potential directed odd cycles on vertices {0,...,n-1}.
    Return as list of tuples (vertex sequence in cycle order).
    For each cycle, we pick the canonical representative (lexicographically smallest rotation)."""
    cycles = []

    # 3-cycles: (a,b,c) with a→b→c→a
    for a in range(n):
        for b in range(n):
            if b == a: continue
            for c in range(n):
                if c == a or c == b: continue
                # Canonical: start with smallest vertex
                cycle = (a, b, c)
                rotations = [cycle, (b,c,a), (c,a,b)]
                canon = min(rotations)
                if cycle == canon:
                    cycles.append(cycle)

    # 5-cycles: (a,b,c,d,e) with a→b→c→d→e→a
    if n >= 5:
        for combo in combinations(range(n), 5):
            from itertools import permutations
            for perm in permutations(combo):
                rotations = [tuple(perm[i:] + perm[:i]) for i in range(5)]
                canon = min(rotations)
                if tuple(perm) == canon:
                    cycles.append(tuple(perm))

    # 7-cycles: only at n>=7
    if n >= 7:
        for combo in combinations(range(n), 7):
            from itertools import permutations
            for perm in permutations(combo):
                rotations = [tuple(perm[i:] + perm[:i]) for i in range(7)]
                canon = min(rotations)
                if tuple(perm) == canon:
                    cycles.append(tuple(perm))

    return cycles

def cycle_indicator(cycle, n, tiles):
    """
    Compute the indicator χ_C(t) as a multilinear polynomial.

    χ_C(t) = Π_{arc (u,v) in C} [arc u→v present in T(t)]

    Each arc (u,v) is either:
    - Base path: u = v+1 (0-indexed), always present → factor 1
    - Tile forward: u > v+1, u→v present when t_k=0 → factor (1-t_k)
    - Tile backward: v > u+1, v→u direction, present when t_k=1 → factor t_k
      Actually: tile (x,y) with x-1=max(u,v), y-1=min(u,v).
      If u > v (forward): arc u→v present when t_k=0, i.e., (1-t_k)
      If u < v: arc u→v, but tile (v+1,u+1) has skip v-u.
        Tile (v+1,u+1) forward means (v+1-1)→(u+1-1) = v→u, NOT u→v.
        So u→v is the BACKWARD direction of tile (v+1,u+1), present when t_k=1.

    Returns: list of (coefficient, tile_subset) pairs representing the polynomial.
    """
    tile_to_idx = {t: i for i, t in enumerate(tiles)}
    m = len(tiles)

    # For each arc in the cycle, determine the factor
    L = len(cycle)
    factors = []  # list of (tile_index, is_backward)
    for i in range(L):
        u = cycle[i]
        v = cycle[(i+1) % L]
        # Is this a base-path arc?
        if u == v + 1:
            # Base path: u→v always present
            continue
        elif v == u + 1:
            # Reverse base path: u→v = u→u+1, never a base path arc
            # Actually base path goes k+1→k. So u→u+1 is NOT base path.
            pass

        # It's a tile arc
        high = max(u, v)
        low = min(u, v)
        tile_label = (high + 1, low + 1)  # 1-indexed

        if tile_label not in tile_to_idx:
            # This shouldn't happen for valid tiles
            return None  # Cycle uses a base-path arc incorrectly

        k = tile_to_idx[tile_label]

        if u > v:
            # Forward direction of tile: u→v present when t_k=0 → factor (1-t_k)
            factors.append((k, False))  # factor = (1-t_k)
        else:
            # Backward direction: u→v present when t_k=1 → factor t_k
            factors.append((k, True))  # factor = t_k

    # Now expand the product of factors
    # Each factor is (1-t_k) [is_backward=False] or t_k [is_backward=True]
    # The product is a multilinear polynomial in the involved tile variables

    # Start with constant 1
    # For each factor, multiply:
    #   (1-t_k): each existing term gets two copies: one unchanged, one with -t_k
    #   t_k: each existing term gets multiplied by t_k

    # Represent as dict: frozenset of tile indices → coefficient
    poly = {frozenset(): 1}

    for k, is_backward in factors:
        new_poly = {}
        if is_backward:
            # Multiply by t_k
            for S, coeff in poly.items():
                new_S = S | {k}
                if new_S in new_poly:
                    new_poly[new_S] += coeff
                else:
                    new_poly[new_S] = coeff
        else:
            # Multiply by (1 - t_k)
            for S, coeff in poly.items():
                # 1 * term
                if S in new_poly:
                    new_poly[S] += coeff
                else:
                    new_poly[S] = coeff
                # -t_k * term
                new_S = S | {k}
                if new_S in new_poly:
                    new_poly[new_S] -= coeff
                else:
                    new_poly[new_S] = -coeff
        poly = {S: c for S, c in new_poly.items() if c != 0}

    return poly

def check_cycle_valid(cycle, t, n, tiles):
    """Check if cycle exists in tournament T(t)."""
    tile_to_idx = {t_label: i for i, t_label in enumerate(tiles)}
    L = len(cycle)
    for i in range(L):
        u = cycle[i]
        v = cycle[(i+1) % L]
        # Check arc u→v exists
        if u == v + 1:
            continue  # base path arc, always exists
        high = max(u, v)
        low = min(u, v)
        tile_label = (high + 1, low + 1)
        if tile_label not in tile_to_idx:
            return False
        k = tile_to_idx[tile_label]
        bit = (t >> k) & 1
        if u > v:
            # Forward: present when bit=0
            if bit != 0:
                return False
        else:
            # Backward: present when bit=1
            if bit != 1:
                return False
    return True

def independent_sets(conflict_adj, cycle_indices):
    """Find all independent sets in the conflict graph among given cycle indices."""
    # conflict_adj[i] = set of j that conflict with i
    n = len(cycle_indices)
    results = [frozenset()]  # empty set is independent

    for i in range(n):
        ci = cycle_indices[i]
        new_results = []
        for IS in results:
            new_results.append(IS)
            # Can we add ci?
            if all(ci not in conflict_adj.get(cj, set()) for cj in IS):
                new_results.append(IS | {ci})
        results = new_results

    return results

def ocf_multilinear_decomposition(n):
    """
    Decompose H(t) = Σ_{I independent} 2^|I| · Π_{C∈I} χ_C(t)
    and extract the multilinear coefficients.
    """
    tiles = make_tiles(n)
    m = len(tiles)
    tile_to_idx = {t: i for i, t in enumerate(tiles)}

    print(f"\n{'='*60}")
    print(f"OCF MULTILINEAR DECOMPOSITION n={n}")
    print(f"{'='*60}")

    # Find all odd cycles and their indicators
    # For n≤7, only 3-cycles matter for the polynomial (5-cycles too complex)
    # Actually, let me start with just 3-cycles for simplicity
    all_3cycles = []
    for a in range(n):
        for b in range(n):
            if b == a: continue
            for c in range(n):
                if c == a or c == b: continue
                verts = tuple(sorted([a,b,c]))
                if verts not in [tuple(sorted(x[:3])) for x in all_3cycles]:
                    # Check both orientations
                    if True:  # a→b→c→a
                        all_3cycles.append((a,b,c))

    # Deduplicate by vertex set
    seen_vsets = set()
    unique_3cycles = []
    for cyc in all_3cycles:
        vset = frozenset(cyc)
        if vset not in seen_vsets:
            seen_vsets.add(vset)
            # There are exactly 2 directed 3-cycles on any 3 vertices (the two orientations)
            unique_3cycles.append(cyc)
            unique_3cycles.append((cyc[0], cyc[2], cyc[1]))  # reverse orientation

    # For each directed 3-cycle, compute its indicator polynomial
    cycle_polys = {}
    for cyc in unique_3cycles:
        poly = cycle_indicator(cyc, n, tiles)
        if poly is not None:
            cycle_polys[cyc] = poly

    print(f"Total directed 3-cycles considered: {len(cycle_polys)}")

    # Build conflict graph (cycles sharing a vertex)
    cycle_list = list(cycle_polys.keys())
    conflict = defaultdict(set)
    for i, c1 in enumerate(cycle_list):
        for j, c2 in enumerate(cycle_list):
            if i >= j: continue
            if set(c1) & set(c2):  # share a vertex
                conflict[i].add(j)
                conflict[j].add(i)

    # Find all independent sets
    print(f"Finding independent sets...")
    ind_sets = independent_sets(conflict, list(range(len(cycle_list))))
    print(f"Total independent sets: {len(ind_sets)}")

    # For each independent set I, compute the product polynomial Π χ_C
    # and weight by 2^|I|
    total_poly = defaultdict(int)

    for IS in ind_sets:
        weight = 2 ** len(IS)

        # Product of cycle indicators for cycles in IS
        prod = {frozenset(): 1}
        for idx in IS:
            cyc = cycle_list[idx]
            chi = cycle_polys[cyc]
            new_prod = defaultdict(int)
            for S1, c1 in prod.items():
                for S2, c2 in chi.items():
                    S = S1 | S2
                    new_prod[S] += c1 * c2
            prod = {S: c for S, c in new_prod.items() if c != 0}

        # Add weighted product to total
        for S, c in prod.items():
            total_poly[S] += weight * c

    # Convert tile-set keys to bitmasks
    coeffs_from_ocf = {}
    for S, c in total_poly.items():
        if c != 0:
            bitmask = 0
            for idx in S:
                bitmask |= (1 << idx)
            coeffs_from_ocf[bitmask] = c

    # Compare with direct computation
    from h_multilinear_complete import compute_all_H, mobius_inversion
    H, _, _ = compute_all_H(n)
    coeffs_direct = mobius_inversion(H, m)

    print(f"\n--- COMPARISON ---")
    print(f"OCF coefficients: {len(coeffs_from_ocf)} nonzero")
    print(f"Direct coefficients: {len(coeffs_direct)} nonzero")

    # Check if they match
    all_keys = set(coeffs_from_ocf.keys()) | set(coeffs_direct.keys())
    mismatches = 0
    for key in sorted(all_keys):
        c_ocf = coeffs_from_ocf.get(key, 0)
        c_dir = coeffs_direct.get(key, 0)
        if c_ocf != c_dir:
            mismatches += 1
            tile_set = [tiles[i] for i in range(m) if key & (1 << i)]
            if mismatches <= 10:
                print(f"  MISMATCH at S={key:0{m}b} tiles={tile_set}: "
                      f"OCF={c_ocf}, direct={c_dir}, diff={c_dir-c_ocf}")

    if mismatches == 0:
        print(f"  PERFECT MATCH! OCF decomposition reproduces all coefficients.")
    else:
        print(f"  {mismatches} mismatches (3-cycles only; 5-cycles needed)")

    # Analyze which coefficients come from which independent set sizes
    by_IS_size = defaultdict(lambda: defaultdict(int))
    for IS in ind_sets:
        weight = 2 ** len(IS)
        prod = {frozenset(): 1}
        for idx in IS:
            cyc = cycle_list[idx]
            chi = cycle_polys[cyc]
            new_prod = defaultdict(int)
            for S1, c1 in prod.items():
                for S2, c2 in chi.items():
                    S = S1 | S2
                    new_prod[S] += c1 * c2
            prod = {S: c for S, c in new_prod.items() if c != 0}

        for S, c in prod.items():
            bitmask = 0
            for idx in S:
                bitmask |= (1 << idx)
            by_IS_size[len(IS)][bitmask] += weight * c

    print(f"\n--- CONTRIBUTION BY INDEPENDENT SET SIZE ---")
    for k in sorted(by_IS_size.keys()):
        entries = by_IS_size[k]
        nz = sum(1 for c in entries.values() if c != 0)
        total = sum(entries.values())
        print(f"  |I|={k}: contributes to {nz} coefficients, total weight sum={total}")

    return coeffs_from_ocf, coeffs_direct

if __name__ == "__main__":
    # Start with n=4 (simple)
    ocf_multilinear_decomposition(4)

    # n=5
    ocf_multilinear_decomposition(5)
