"""
h_ie_even_n_cancel.py — opus-2026-04-04-S3

THE KEY OPEN QUESTION: At even n=6, why are degree-5 coefficients zero?

Degree-5 contributions come from:
  Source A: Single 5-cycles with 5 tile arcs (cancel by reversal, 5 odd)
  Source B: Single 6-cycles with 5 tile arcs (1 base-path arc) — does NOT cancel!
  Source C: Pairs of 3-cycles with combined 5 tile arcs — weight 4

Hypothesis: Source B and Source C cancel each other perfectly.
This would be a DEEP cancellation between single long cycles and cycle PAIRS.
"""
import numpy as np
from collections import defaultdict
from h_ocf_multilinear_v2 import make_tiles, enumerate_directed_cycles, cycle_indicator

def analyze_degree5_cancellation():
    n = 6
    tiles = make_tiles(n)
    m = len(tiles)
    tile_to_idx = {t: i for i, t in enumerate(tiles)}
    d_target = 5  # Degree we're investigating

    print(f"EVEN-n DEGREE CANCELLATION ANALYSIS")
    print(f"n={n}, investigating degree {d_target}")
    print(f"{'='*60}")

    # Enumerate all cycles
    all_cycles = []
    for k in [3, 5]:
        for cyc in enumerate_directed_cycles(n, k):
            chi = cycle_indicator(cyc, n, tiles, tile_to_idx)
            if chi is not None:
                max_deg = max(len(S) for S in chi.keys()) if chi else 0
                all_cycles.append((cyc, chi, k, max_deg))

    # Now enumerate 6-cycles (Hamiltonian cycles at n=6)
    ham_cycles = []
    for cyc in enumerate_directed_cycles(n, 6):
        chi = cycle_indicator(cyc, n, tiles, tile_to_idx)
        if chi is not None:
            max_deg = max(len(S) for S in chi.keys()) if chi else 0
            ham_cycles.append((cyc, chi, 6, max_deg))
            all_cycles.append((cyc, chi, 6, max_deg))

    print(f"Valid cycles: {len(all_cycles)}")
    print(f"  3-cycles: {sum(1 for _,_,k,_ in all_cycles if k==3)}")
    print(f"  5-cycles: {sum(1 for _,_,k,_ in all_cycles if k==5)}")
    print(f"  6-cycles: {sum(1 for _,_,k,_ in all_cycles if k==6)}")

    # Degree distribution of 6-cycle indicators
    ham_deg_dist = defaultdict(int)
    for _, _, _, d in ham_cycles:
        ham_deg_dist[d] += 1
    print(f"\n6-cycle indicator degrees: {dict(sorted(ham_deg_dist.items()))}")

    # Source B: 6-cycles with degree >= 5
    source_b_cycles = [(cyc, chi) for cyc, chi, k, d in ham_cycles if d >= 5]
    print(f"\nSource B: {len(source_b_cycles)} 6-cycles with degree >= 5")

    # Collect degree-5 contributions from ALL sources
    # Build conflict graph
    nc = len(all_cycles)
    conflict = [set() for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if set(all_cycles[i][0]) & set(all_cycles[j][0]):
                conflict[i].add(j)
                conflict[j].add(i)

    # For each 5-tile subset, compute contributions from:
    # |I|=1 (single cycles): weight 2
    # |I|=2 (cycle pairs): weight 4

    from h_multilinear_complete import compute_all_H, mobius_inversion
    H, _, _ = compute_all_H(n)
    coeffs_direct = mobius_inversion(H, m)

    # Check all degree-5 subsets
    from itertools import combinations
    all_d5_subsets = []
    for combo in combinations(range(m), d_target):
        S = frozenset(combo)
        S_mask = sum(1 << k for k in combo)

        # Direct coefficient (should be 0)
        c_direct = coeffs_direct.get(S_mask, 0)

        # Source A: 5-cycles with 5 tile arcs (odd, cancel by reversal)
        source_a = 0
        for idx, (cyc, chi, k, d) in enumerate(all_cycles):
            if k != 5: continue
            c = chi.get(S, 0)
            source_a += 2 * c

        # Source B: 6-cycles contributing degree-5 term
        source_b = 0
        for idx, (cyc, chi, k, d) in enumerate(all_cycles):
            if k != 6: continue
            c = chi.get(S, 0)
            source_b += 2 * c

        # Source C: 3-cycle pairs with combined degree 5
        source_c = 0
        three_cycles = [(idx, cyc, chi) for idx, (cyc, chi, k, d) in enumerate(all_cycles) if k == 3]
        for i1, (idx1, cyc1, chi1) in enumerate(three_cycles):
            for i2, (idx2, cyc2, chi2) in enumerate(three_cycles):
                if i2 <= i1: continue
                if idx2 in conflict[idx1]: continue  # vertex overlap
                # Product coefficient for S
                c = 0
                for S1, v1 in chi1.items():
                    S2 = S - S1
                    if S1.isdisjoint(S2) and S1 | S2 == S:
                        v2 = chi2.get(S2, 0)
                        c += v1 * v2
                source_c += 4 * c

        # Source D: pairs involving 5-cycles or 6-cycles (unlikely to contribute at degree 5)
        # Actually: a 5-cycle uses 5 vertices, paired with nothing useful at n=6 (only 1 vertex left)
        # A 6-cycle uses all vertices, can't pair with anything
        # So Source C is the only pair source

        total = source_a + source_b + source_c

        all_d5_subsets.append({
            'S': combo, 'direct': c_direct,
            'source_a': source_a, 'source_b': source_b,
            'source_c': source_c, 'total': total
        })

    # Display results
    print(f"\n{'Tiles':<35} {'Direct':>7} {'SrcA(5cy)':>10} {'SrcB(6cy)':>10} {'SrcC(pair)':>10} {'Total':>7} {'Match':>6}")
    print("-" * 95)

    match_count = 0
    nonzero_count = 0
    for entry in all_d5_subsets:
        match = entry['direct'] == entry['total']
        if match: match_count += 1
        if entry['source_b'] != 0 or entry['source_c'] != 0:
            nonzero_count += 1
            ts = [tiles[k] for k in entry['S']]
            print(f"{str(ts):<35} {entry['direct']:>7} {entry['source_a']:>10} "
                  f"{entry['source_b']:>10} {entry['source_c']:>10} "
                  f"{entry['total']:>7} {'✓' if match else '✗':>6}")

    print(f"\n  Total degree-5 subsets: {len(all_d5_subsets)}")
    print(f"  Match direct: {match_count}/{len(all_d5_subsets)}")
    print(f"  Subsets with nonzero Source B or C: {nonzero_count}")

    # KEY QUESTION: Does Source B + Source C = 0 for ALL degree-5 subsets?
    all_cancel = all(e['source_b'] + e['source_c'] == 0 for e in all_d5_subsets)
    print(f"\n  *** Source B + Source C = 0 always: {all_cancel} ***")

    if all_cancel:
        print(f"  THEOREM: At n=6, degree-5 cancellation comes from 6-cycle")
        print(f"  contributions perfectly canceling against 3-cycle pair contributions!")

    # Check degree 6 too
    print(f"\n{'='*60}")
    print(f"DEGREE 6 ANALYSIS")
    print(f"{'='*60}")

    for combo in combinations(range(m), 6):
        S = frozenset(combo)
        S_mask = sum(1 << k for k in combo)
        c_direct = coeffs_direct.get(S_mask, 0)

        source_6cy = 0
        for idx, (cyc, chi, k, d) in enumerate(all_cycles):
            if k != 6: continue
            c = chi.get(S, 0)
            source_6cy += 2 * c

        source_pair = 0
        for i1, (idx1, cyc1, chi1) in enumerate(three_cycles):
            for i2, (idx2, cyc2, chi2) in enumerate(three_cycles):
                if i2 <= i1: continue
                if idx2 in conflict[idx1]: continue
                c = 0
                for S1, v1 in chi1.items():
                    S2 = S - S1
                    if S1.isdisjoint(S2) and S1 | S2 == S:
                        v2 = chi2.get(S2, 0)
                        c += v1 * v2
                source_pair += 4 * c

        total = source_6cy + source_pair
        if source_6cy != 0 or source_pair != 0:
            ts = [tiles[k] for k in combo]
            print(f"  {ts}: 6cy={source_6cy}, pair={source_pair}, total={total}")

    # Summary
    d6_subsets = list(combinations(range(m), 6))
    d6_cancel = []
    for combo in d6_subsets:
        S = frozenset(combo)
        s6 = sum(2 * all_cycles[idx][1].get(S, 0) for idx, (cyc, chi, k, d) in enumerate(all_cycles) if k == 6)
        sp = 0
        for i1, (idx1, cyc1, chi1) in enumerate(three_cycles):
            for i2, (idx2, cyc2, chi2) in enumerate(three_cycles):
                if i2 <= i1: continue
                if idx2 in conflict[idx1]: continue
                c = 0
                for S1, v1 in chi1.items():
                    S2 = S - S1
                    if S1.isdisjoint(S2) and S1 | S2 == S:
                        v2 = chi2.get(S2, 0)
                        c += v1 * v2
                sp += 4 * c
        d6_cancel.append(s6 + sp)

    all_d6_zero = all(c == 0 for c in d6_cancel)
    print(f"\n  All degree-6 (6cy + 3cy-pair) = 0: {all_d6_zero}")

if __name__ == "__main__":
    analyze_degree5_cancellation()
