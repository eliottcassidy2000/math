"""
h_ocf_max_degree.py — opus-2026-04-04-S2

Analyze the maximum-degree multilinear coefficients via the OCF decomposition.

THM-287 shows: degree-k coefficient c_S gets contributions from independent
sets of size 1,...,k only. At maximum degree d_max = 2*floor((n-1)/2):

Each contributing independent set I has |I| ≤ d_max, and each cycle in I
contributes ≥1 to the degree. So we need sum of per-cycle degrees = d_max.

Since d_max is the max degree and the minimum cycle degree is 1 (one tile arc),
the most efficient packing is |I| = d_max single-tile cycles. But this requires
d_max vertex-disjoint odd cycles, each using only 1 non-base-path arc.

A 3-cycle with exactly 1 non-base-path arc uses 2 base-path arcs (consecutive
vertices) and 1 tile arc. The tile arc "jumps" skip ≥ 2 vertices.

The question: which independent set sizes actually contribute to degree d_max?
"""
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
from h_ocf_multilinear_v2 import make_tiles, enumerate_directed_cycles, cycle_indicator

def max_degree_anatomy(n, max_cycle_len=None):
    tiles = make_tiles(n)
    m = len(tiles)
    tile_to_idx = {t: i for i, t in enumerate(tiles)}

    if max_cycle_len is None:
        max_cycle_len = n if n % 2 == 1 else n - 1

    d_max = 2 * ((n-1) // 2)
    print(f"\n{'='*60}")
    print(f"MAX DEGREE ANATOMY n={n}, d_max={d_max}")
    print(f"{'='*60}")

    # Enumerate cycles and their indicators
    all_cycles = []
    for k in range(3, max_cycle_len + 1, 2):
        for cyc in enumerate_directed_cycles(n, k):
            chi = cycle_indicator(cyc, n, tiles, tile_to_idx)
            if chi is not None:
                # Compute the degree (max order of chi's terms)
                if chi:
                    max_deg = max(len(S) for S in chi.keys())
                else:
                    max_deg = 0
                all_cycles.append((cyc, chi, max_deg))

    print(f"Total valid cycles: {len(all_cycles)}")

    # Distribution of cycle polynomial degrees
    deg_dist = defaultdict(int)
    for _, _, d in all_cycles:
        deg_dist[d] += 1
    print(f"Cycle indicator degree distribution: {dict(sorted(deg_dist.items()))}")

    # For each degree-d_max tile subset S that has nonzero c_S,
    # determine which independent sets contribute.
    from h_multilinear_complete import compute_all_H, mobius_inversion
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    max_order_coeffs = {S: c for S, c in coeffs.items() if bin(S).count('1') == d_max}
    print(f"\nDegree-{d_max} coefficients: {len(max_order_coeffs)} nonzero")
    if not max_order_coeffs:
        return

    # Build conflict graph
    nc = len(all_cycles)
    conflict = [set() for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if set(all_cycles[i][0]) & set(all_cycles[j][0]):
                conflict[i].add(j)
                conflict[j].add(i)

    # For a sample of max-degree tile subsets, trace which IS contribute
    sample = list(max_order_coeffs.items())[:5]

    for S_mask, c_target in sample:
        S = frozenset(k for k in range(m) if S_mask & (1 << k))
        tile_names = [tiles[k] for k in sorted(S)]
        print(f"\n  S = {tile_names}, c_S = {c_target}")

        # Find all independent sets that contribute to this coefficient
        # We need: Π_{C∈I} χ_C has a term matching exactly the tile set S
        # And since each χ_C has no constant term, |I| ≤ |S| = d_max

        # For efficiency, only check IS up to size d_max
        # For each IS, compute the product and check for the S term

        contributions = defaultdict(int)  # |I| → total contribution

        # |I|=1: single cycle
        for idx, (cyc, chi, deg) in enumerate(all_cycles):
            c = chi.get(S, 0)
            if c != 0:
                contributions[1] += 2 * c

        # |I|=2: pairs
        for idx1 in range(nc):
            for idx2 in range(idx1+1, nc):
                if idx2 in conflict[idx1]:
                    continue
                chi1 = all_cycles[idx1][1]
                chi2 = all_cycles[idx2][1]

                # Coefficient of Π_{k∈S} t_k in χ₁·χ₂
                # = Σ_{S₁∪S₂=S, S₁∩S₂=∅} chi1[S₁] * chi2[S₂]
                c = 0
                for S1, c1 in chi1.items():
                    S2 = S - S1
                    if S1 & S2 == frozenset() and S1 | S2 == S:
                        c2 = chi2.get(S2, 0)
                        c += c1 * c2
                if c != 0:
                    contributions[2] += 4 * c

        # |I|=3: triples (only if n >= 9)
        if n >= 9:
            pass  # Skip for now

        total = sum(contributions.values())
        print(f"    Contributions by |I|: {dict(sorted(contributions.items()))}")
        print(f"    Total = {total}, {'match' if total == c_target else f'MISMATCH (expected {c_target})'}")

    # IMPORTANT QUESTION: at d_max, do |I|=1 or |I|=2 dominate?
    print(f"\n--- DOMINANCE ANALYSIS ---")
    is1_count = 0
    is2_count = 0
    both_count = 0

    for S_mask, c_target in max_order_coeffs.items():
        S = frozenset(k for k in range(m) if S_mask & (1 << k))

        c_from_1 = 0
        for idx, (cyc, chi, deg) in enumerate(all_cycles):
            c_from_1 += 2 * chi.get(S, 0)

        c_from_2 = 0
        for idx1 in range(nc):
            for idx2 in range(idx1+1, nc):
                if idx2 in conflict[idx1]:
                    continue
                chi1 = all_cycles[idx1][1]
                chi2 = all_cycles[idx2][1]
                c = 0
                for S1, c1 in chi1.items():
                    S2 = S - S1
                    if S1 & S2 == frozenset() and S1 | S2 == S:
                        c += c1 * chi2.get(S2, 0)
                c_from_2 += 4 * c

        if c_from_1 != 0 and c_from_2 == 0:
            is1_count += 1
        elif c_from_1 == 0 and c_from_2 != 0:
            is2_count += 1
        elif c_from_1 != 0 and c_from_2 != 0:
            both_count += 1

    total_max = len(max_order_coeffs)
    print(f"  |I|=1 only: {is1_count}/{total_max}")
    print(f"  |I|=2 only: {is2_count}/{total_max}")
    print(f"  Both |I|=1 and |I|=2: {both_count}/{total_max}")

if __name__ == "__main__":
    max_degree_anatomy(5)
    max_degree_anatomy(6, max_cycle_len=5)
    max_degree_anatomy(7, max_cycle_len=5)  # Skip 7-cycles for speed
