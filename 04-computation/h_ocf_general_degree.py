"""
h_ocf_general_degree.py — opus-2026-04-04-S2

The complete OCF picture for ALL multilinear coefficients.

THM-287: degree-k coefficient c_S = Σ_{|I|=1}^{k} 2^|I| · [product contribution]

Key question: at each degree k, what is the relative contribution of
different independent set sizes? And can we derive CLOSED FORMULAS
for the coefficients?

SPECIFIC GOAL: Derive the formula for same-end shared-vertex pairs
c_{ij} = -2^max(1,|s₁-s₂|-1) from the OCF.
"""
import numpy as np
from collections import defaultdict
from h_ocf_multilinear_v2 import make_tiles, enumerate_directed_cycles, cycle_indicator
from fractions import Fraction

def same_end_formula(n, max_cycle_len=None):
    """Derive the same-end shared-vertex formula from OCF."""
    tiles = make_tiles(n)
    m = len(tiles)
    tile_to_idx = {t: i for i, t in enumerate(tiles)}

    if max_cycle_len is None:
        max_cycle_len = n if n % 2 == 1 else n - 1

    print(f"\n{'='*60}")
    print(f"SAME-END FORMULA DERIVATION n={n}")
    print(f"{'='*60}")

    # Enumerate cycles
    all_cycles = []
    for k in range(3, max_cycle_len + 1, 2):
        for cyc in enumerate_directed_cycles(n, k):
            chi = cycle_indicator(cyc, n, tiles, tile_to_idx)
            if chi is not None:
                all_cycles.append((cyc, chi))

    # For each same-end pair, find which cycles contribute
    for i in range(m):
        for j in range(i+1, m):
            t1, t2 = tiles[i], tiles[j]
            x1, y1 = t1
            x2, y2 = t2
            s1, s2 = x1-y1, x2-y2

            # Check: same upper or same lower vertex
            shared_type = None
            if x1 == x2:
                shared_type = f"upper={x1}"
            elif y1 == y2:
                shared_type = f"lower={y1}"
            else:
                continue  # not same-end

            # Compute c_{ij} from OCF (single cycles only, since shared vertex)
            c_ij = 0
            contributing_cycles = []
            for cyc, chi in all_cycles:
                c = chi.get(frozenset({i, j}), 0)
                if c != 0:
                    c_ij += 2 * c
                    # Classify the cycle
                    L = len(cyc)
                    # How many tile arcs?
                    non_base = sum(1 for k in range(L) if abs(cyc[k]-cyc[(k+1)%L]) != 1)
                    contributing_cycles.append((cyc, c, L, non_base))

            expected = -2**max(1, abs(s1-s2)-1)
            match = "✓" if c_ij == expected else f"✗ ({expected})"

            print(f"\n  ({x1},{y1})-({x2},{y2}) s=({s1},{s2}) {shared_type}:")
            print(f"    c₂ = {c_ij} {match}")
            print(f"    Contributing cycles ({len(contributing_cycles)}):")
            for cyc, c, L, nb in contributing_cycles:
                print(f"      {L}-cycle {cyc}: chi coeff = {c:+d}, non-base arcs = {nb}")

            # WHY is c_ij = -2^max(1,|Δs|-1)?
            # Each contributing cycle gives ±1. The total is c_ij/2.
            # So c_ij/2 = net count = (#positive) - (#negative)
            net = c_ij // 2
            pos = sum(1 for _,c,_,_ in contributing_cycles if c > 0)
            neg = sum(1 for _,c,_,_ in contributing_cycles if c < 0)
            total = pos + neg
            print(f"    Net count = {net} = {pos}(+) - {neg}(-) = {pos-neg}")

def cross_end_formula(n, max_cycle_len=None):
    """Derive the cross-end formula c₂ = +2^(s₁+s₂-2) from OCF."""
    tiles = make_tiles(n)
    m = len(tiles)
    tile_to_idx = {t: i for i, t in enumerate(tiles)}

    if max_cycle_len is None:
        max_cycle_len = n if n % 2 == 1 else n - 1

    print(f"\n{'='*60}")
    print(f"CROSS-END FORMULA DERIVATION n={n}")
    print(f"{'='*60}")

    all_cycles = []
    for k in range(3, max_cycle_len + 1, 2):
        for cyc in enumerate_directed_cycles(n, k):
            chi = cycle_indicator(cyc, n, tiles, tile_to_idx)
            if chi is not None:
                all_cycles.append((cyc, chi))

    for i in range(m):
        for j in range(i+1, m):
            t1, t2 = tiles[i], tiles[j]
            x1, y1 = t1
            x2, y2 = t2
            s1, s2 = x1-y1, x2-y2

            # Cross-end: x1 = y2 or y1 = x2
            if x1 == y2:
                shared_type = f"x1=y2={x1}"
            elif y1 == x2:
                shared_type = f"y1=x2={y1}"
            else:
                continue

            c_ij = 0
            contributing_cycles = []
            for cyc, chi in all_cycles:
                c = chi.get(frozenset({i, j}), 0)
                if c != 0:
                    c_ij += 2 * c
                    contributing_cycles.append((cyc, c, len(cyc)))

            expected = 2**(s1+s2-2)
            match = "✓" if c_ij == expected else f"✗ ({expected})"

            print(f"\n  ({x1},{y1})-({x2},{y2}) s=({s1},{s2}) {shared_type}:")
            print(f"    c₂ = {c_ij} {match}")
            net = c_ij // 2
            print(f"    Net cycle count = {net} = 2^({s1}+{s2}-2)/2 = 2^{s1+s2-3}")
            print(f"    Contributing cycles ({len(contributing_cycles)}):")
            for cyc, c, L in contributing_cycles[:6]:
                print(f"      {L}-cycle {cyc}: chi coeff = {c:+d}")

def degree_k_source_analysis(n, max_cycle_len=None):
    """For each degree k, breakdown of contributions by |I| and cycle length."""
    tiles = make_tiles(n)
    m = len(tiles)
    tile_to_idx = {t: i for i, t in enumerate(tiles)}

    if max_cycle_len is None:
        max_cycle_len = n if n % 2 == 1 else n - 1

    all_cycles = []
    for k in range(3, max_cycle_len + 1, 2):
        for cyc in enumerate_directed_cycles(n, k):
            chi = cycle_indicator(cyc, n, tiles, tile_to_idx)
            if chi is not None:
                all_cycles.append((cyc, chi))

    nc = len(all_cycles)
    conflict = [set() for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if set(all_cycles[i][0]) & set(all_cycles[j][0]):
                conflict[i].add(j)
                conflict[j].add(i)

    print(f"\n{'='*60}")
    print(f"DEGREE-k SOURCE ANALYSIS n={n}")
    print(f"{'='*60}")

    from h_multilinear_complete import compute_all_H, mobius_inversion
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    # Group coefficients by degree
    by_degree = defaultdict(list)
    for S, c in coeffs.items():
        by_degree[bin(S).count('1')].append((S, c))

    for deg in sorted(by_degree.keys()):
        entries = by_degree[deg]
        if deg == 0:
            print(f"\n  Degree {deg}: constant = {entries[0][1]}")
            continue

        # Classify each coefficient by which |I| sizes contribute
        source_types = defaultdict(int)
        for S_mask, c_target in entries:
            S = frozenset(k for k in range(m) if S_mask & (1 << k))

            # |I|=1
            c1 = 0
            for idx in range(nc):
                c1 += 2 * all_cycles[idx][1].get(S, 0)

            # |I|=2 (if deg >= 2)
            c2 = 0
            if deg >= 2:
                for i1 in range(nc):
                    for i2 in range(i1+1, nc):
                        if i2 in conflict[i1]: continue
                        chi1 = all_cycles[i1][1]
                        chi2 = all_cycles[i2][1]
                        for S1, v1 in chi1.items():
                            S2 = S - S1
                            if S1 & S2 == frozenset() and S1 | S2 == S:
                                v2 = chi2.get(S2, 0)
                                c2 += 4 * v1 * v2

            if c1 == c_target:
                source_types['|I|=1 only'] += 1
            elif c2 == c_target and c1 == 0:
                source_types['|I|=2 only'] += 1
            elif c1 + c2 == c_target:
                source_types['|I|=1+2'] += 1
            else:
                source_types['UNKNOWN (need |I|>=3?)'] += 1

        print(f"\n  Degree {deg}: {len(entries)} coefficients")
        for src, count in sorted(source_types.items()):
            print(f"    {src}: {count}")

if __name__ == "__main__":
    same_end_formula(6, max_cycle_len=5)
    cross_end_formula(6, max_cycle_len=5)
    degree_k_source_analysis(5)
    degree_k_source_analysis(6, max_cycle_len=5)
