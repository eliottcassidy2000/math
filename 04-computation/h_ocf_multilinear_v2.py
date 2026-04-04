"""
h_ocf_multilinear_v2.py — opus-2026-04-04-S2

Complete OCF multilinear decomposition including ALL odd cycles (3,5,7,...).

H(t) = Σ_{I independent in Ω(t)} 2^|I| = Σ_{I vertex-disjoint odd cycles} 2^|I| · Π_{C∈I} χ_C(t)

This should reproduce H(t) EXACTLY.
"""
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict

def make_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

def enumerate_directed_cycles(n, k):
    """Enumerate all directed k-cycles on {0,...,n-1}.
    A directed k-cycle (v₀,...,v_{k-1}) has arcs v_i → v_{(i+1) mod k}.
    Two sequences are the same cycle iff one is a rotation of the other.
    Return canonical representatives (lexicographically smallest rotation)."""
    cycles = set()
    for verts in combinations(range(n), k):
        for perm in permutations(verts):
            # Check canonical (lex-smallest rotation)
            rotations = [perm[i:] + perm[:i] for i in range(k)]
            canon = min(rotations)
            if perm == canon:
                cycles.add(perm)
    return list(cycles)

def cycle_indicator(cycle, n, tiles, tile_to_idx):
    """Compute χ_C(t) = Π_{arcs of C} [arc present in T(t)].

    Returns polynomial as dict: frozenset of tile indices → coefficient.
    Returns None if any arc is a REVERSE base-path arc (never exists → cycle impossible).
    """
    m = len(tiles)
    L = len(cycle)

    # Start with constant 1
    poly = {frozenset(): 1}

    for i in range(L):
        u = cycle[i]
        v = cycle[(i+1) % L]

        # Determine what kind of arc u→v is
        if u == v + 1:
            # Base path arc (u→v = (k+1)→k): always present
            continue
        elif v == u + 1:
            # Reverse base path: u→v = k→(k+1). This is NEVER present
            # (base path goes the other way, and there's no tile for consecutive vertices).
            return None  # Cycle can never exist

        # Non-consecutive: it's a tile arc
        high = max(u, v)
        low = min(u, v)
        tile_label = (high + 1, low + 1)  # Convert to 1-indexed

        if tile_label not in tile_to_idx:
            # Shouldn't happen for non-consecutive vertices
            return None

        k = tile_to_idx[tile_label]

        # Multiply polynomial by the appropriate factor
        new_poly = defaultdict(int)
        if u > v:
            # Forward direction: arc u→v present when t_k = 0 → factor (1 - t_k)
            for S, c in poly.items():
                new_poly[S] += c           # × 1
                new_poly[S | {k}] -= c     # × (-t_k)
        else:
            # Backward direction: arc u→v (with u < v) present when t_k = 1 → factor t_k
            for S, c in poly.items():
                new_poly[S | {k}] += c     # × t_k

        poly = {S: c for S, c in new_poly.items() if c != 0}

    return poly

def ocf_decomposition(n, max_cycle_len=None):
    """Full OCF multilinear decomposition."""
    tiles = make_tiles(n)
    m = len(tiles)
    tile_to_idx = {t: i for i, t in enumerate(tiles)}

    if max_cycle_len is None:
        max_cycle_len = n if n % 2 == 1 else n - 1

    print(f"\n{'='*60}")
    print(f"OCF MULTILINEAR DECOMPOSITION n={n}, max cycle length={max_cycle_len}")
    print(f"m={m} tiles")
    print(f"{'='*60}")

    # Enumerate all odd-length directed cycles with valid indicators
    all_cycles = []
    for k in range(3, max_cycle_len + 1, 2):
        raw = enumerate_directed_cycles(n, k)
        valid = 0
        for cyc in raw:
            chi = cycle_indicator(cyc, n, tiles, tile_to_idx)
            if chi is not None:
                all_cycles.append((cyc, chi))
                valid += 1
        print(f"  {k}-cycles: {len(raw)} enumerated, {valid} have valid indicators")

    print(f"  Total valid directed odd cycles: {len(all_cycles)}")

    # Build conflict graph (cycles sharing a vertex)
    nc = len(all_cycles)
    conflict = [set() for _ in range(nc)]
    for i in range(nc):
        vi = set(all_cycles[i][0])
        for j in range(i+1, nc):
            vj = set(all_cycles[j][0])
            if vi & vj:
                conflict[i].add(j)
                conflict[j].add(i)

    # Find all independent sets via backtracking
    print(f"  Finding independent sets...", flush=True)
    ind_sets = []

    def backtrack(start, current_IS, forbidden):
        ind_sets.append(frozenset(current_IS))
        for i in range(start, nc):
            if i in forbidden:
                continue
            new_forbidden = forbidden | conflict[i]
            backtrack(i + 1, current_IS + [i], new_forbidden)

    backtrack(0, [], set())
    print(f"  Total independent sets: {len(ind_sets)}")

    # Compute H(t) = Σ_I 2^|I| · Π_{C∈I} χ_C(t)
    total_poly = defaultdict(int)
    contribution_by_size = defaultdict(lambda: defaultdict(int))

    for IS in ind_sets:
        weight = 2 ** len(IS)

        # Product of indicators
        prod = {frozenset(): 1}
        for idx in IS:
            _, chi = all_cycles[idx]
            new_prod = defaultdict(int)
            for S1, c1 in prod.items():
                for S2, c2 in chi.items():
                    new_prod[S1 | S2] += c1 * c2
            prod = {S: c for S, c in new_prod.items() if c != 0}

        for S, c in prod.items():
            total_poly[S] += weight * c
            contribution_by_size[len(IS)][S] += weight * c

    # Convert to bitmask representation
    coeffs_ocf = {}
    for S, c in total_poly.items():
        if c != 0:
            bitmask = 0
            for idx in S:
                bitmask |= (1 << idx)
            coeffs_ocf[bitmask] = coeffs_ocf.get(bitmask, 0) + c

    # Compare with direct computation
    from h_multilinear_complete import compute_all_H, mobius_inversion
    H, _, _ = compute_all_H(n)
    coeffs_direct = mobius_inversion(H, m)

    all_keys = set(coeffs_ocf.keys()) | set(coeffs_direct.keys())
    mismatches = 0
    for key in sorted(all_keys):
        c_ocf = coeffs_ocf.get(key, 0)
        c_dir = coeffs_direct.get(key, 0)
        if c_ocf != c_dir:
            mismatches += 1

    print(f"\n--- COMPARISON ---")
    print(f"OCF coefficients: {len(coeffs_ocf)} nonzero")
    print(f"Direct coefficients: {len(coeffs_direct)} nonzero")
    if mismatches == 0:
        print(f"  *** PERFECT MATCH! *** OCF decomposition reproduces ALL coefficients.")
    else:
        print(f"  {mismatches} mismatches")
        for key in sorted(all_keys):
            c_ocf = coeffs_ocf.get(key, 0)
            c_dir = coeffs_direct.get(key, 0)
            if c_ocf != c_dir:
                ts = [tiles[i] for i in range(m) if key & (1 << i)]
                print(f"    S={key:0{m}b} tiles={ts}: OCF={c_ocf}, direct={c_dir}")
                if mismatches > 10:
                    break

    # Analyze contributions by independent set size
    print(f"\n--- CONTRIBUTIONS BY |I| ---")
    for k in sorted(contribution_by_size.keys()):
        entries = contribution_by_size[k]
        nz = sum(1 for c in entries.values() if c != 0)
        # Convert to bitmask and check which multilinear orders
        by_order = defaultdict(int)
        for S, c in entries.items():
            order = len(S)
            by_order[order] += abs(c)
        print(f"  |I|={k}: contributes to {nz} tile subsets. "
              f"By multilinear order: {dict(sorted(by_order.items()))}")

    # THE KEY INSIGHT: What cycles contribute to each multilinear coefficient?
    print(f"\n--- CYCLE ANATOMY OF LINEAR COEFFICIENTS ---")
    for i in range(m):
        bitmask = 1 << i
        tile = tiles[i]
        skip = tile[0] - tile[1]
        c_direct = coeffs_direct.get(bitmask, 0)

        # Which single cycles contribute to this tile's linear coefficient?
        contributions = []
        for IS in ind_sets:
            if len(IS) != 1:
                continue
            idx = list(IS)[0]
            cyc, chi = all_cycles[idx]
            # Does chi have a term with S = {i}?
            for S, c in chi.items():
                if S == frozenset({i}):
                    contributions.append((cyc, c, 2))  # weight = 2^1 = 2

        total_from_cycles = sum(w * c for _, c, w in contributions)
        print(f"  tile {tile} skip={skip}: c₁={c_direct}")
        print(f"    From single cycles: {total_from_cycles} "
              f"({'match' if total_from_cycles == c_direct else 'MISMATCH, higher IS contribute'})")
        if len(contributions) <= 8:
            for cyc, c, w in contributions:
                print(f"      cycle {cyc}: chi coefficient for {{tile {i}}} = {c}, "
                      f"contribution = {w*c}")

    return coeffs_ocf, coeffs_direct

if __name__ == "__main__":
    ocf_decomposition(4)
    ocf_decomposition(5)
    # n=6 might be slow due to 5-cycles on 6 vertices (C(6,5)*4!=720 potential cycles)
    # but let's try
    ocf_decomposition(6, max_cycle_len=5)
