"""
h_ocf_quadratic_anatomy.py — opus-2026-04-04-S2

Decompose quadratic multilinear coefficients via OCF.

c_{ij} gets contributions from:
  TYPE A: Single cycles C where χ_C has a quadratic term in (t_i, t_j)
  TYPE B: Pairs {C₁, C₂} of disjoint cycles where C₁ has t_i and C₂ has t_j

At n≤5: only Type A (no disjoint cycle pairs).
At n≥6: both Type A and Type B.

This should EXPLAIN:
1. Why shared-vertex pairs have c₂ = -2^max(1,|Δs|-1) (same-end)
2. Why cross-end pairs have c₂ = +2^(s₁+s₂-2)
3. Why disjoint pairs have complex non-power-of-2 values at n≥6
"""
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
from h_ocf_multilinear_v2 import make_tiles, enumerate_directed_cycles, cycle_indicator

def quadratic_anatomy(n, max_cycle_len=None):
    tiles = make_tiles(n)
    m = len(tiles)
    tile_to_idx = {t: i for i, t in enumerate(tiles)}

    if max_cycle_len is None:
        max_cycle_len = n if n % 2 == 1 else n - 1

    print(f"\n{'='*60}")
    print(f"QUADRATIC COEFFICIENT ANATOMY n={n}")
    print(f"{'='*60}")

    # Enumerate all valid directed odd cycles with their indicators
    all_cycles = []
    for k in range(3, max_cycle_len + 1, 2):
        for cyc in enumerate_directed_cycles(n, k):
            chi = cycle_indicator(cyc, n, tiles, tile_to_idx)
            if chi is not None:
                all_cycles.append((cyc, chi))

    # Build conflict graph
    nc = len(all_cycles)
    conflict = [set() for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if set(all_cycles[i][0]) & set(all_cycles[j][0]):
                conflict[i].add(j)
                conflict[j].add(i)

    # For each tile pair, decompose c₂ into Type A (single cycle) and Type B (pair)
    print(f"\n{'Tile pair':<25} {'Type':<8} {'s1,s2':<8} {'c₂':>6} {'TypeA':>8} {'TypeB':>8}")
    print("-" * 75)

    for i in range(m):
        for j in range(i+1, m):
            t1, t2 = tiles[i], tiles[j]
            x1, y1 = t1
            x2, y2 = t2
            s1, s2 = x1-y1, x2-y2

            # Vertex relationship
            shared = set()
            if x1 == x2: shared.add('upper')
            if y1 == y2: shared.add('lower')
            if x1 == y2: shared.add('cross1')
            if y1 == x2: shared.add('cross2')
            verts1, verts2 = {x1,y1}, {x2,y2}
            if shared:
                rel = 'shared'
            elif not (verts1 & verts2):
                rel = 'disjoint'
            else:
                rel = 'other'

            # TYPE A: single cycles with quadratic term in (t_i, t_j)
            type_a = 0
            for idx, (cyc, chi) in enumerate(all_cycles):
                # Look for term with exactly {i, j}
                for S, c in chi.items():
                    if S == frozenset({i, j}):
                        type_a += 2 * c  # weight 2^1 = 2

            # TYPE B: pairs of disjoint cycles where one has t_i and other has t_j
            type_b = 0
            for idx1 in range(nc):
                for idx2 in range(idx1+1, nc):
                    if idx2 in conflict[idx1]:
                        continue  # Not independent
                    cyc1, chi1 = all_cycles[idx1]
                    cyc2, chi2 = all_cycles[idx2]
                    # Check: chi1 has linear term in t_i, chi2 has linear term in t_j
                    c1_i = chi1.get(frozenset({i}), 0)
                    c2_j = chi2.get(frozenset({j}), 0)
                    c1_j = chi1.get(frozenset({j}), 0)
                    c2_i = chi2.get(frozenset({i}), 0)

                    # Contribution from pair (C1 gives t_i, C2 gives t_j)
                    type_b += 4 * c1_i * c2_j  # weight 2^2 = 4
                    type_b += 4 * c1_j * c2_i  # weight 2^2 = 4

                    # Also: chi1 has t_i*t_j term? No — that's Type A with |I|=1.
                    # For Type B, we need the PRODUCT of linear terms.

            # Hmm, actually the formula is more subtle. Let me redo this.
            # The OCF gives: c_S = Σ_{I independent} 2^|I| · [coeff of prod_{k∈S} t_k in Π_{C∈I} χ_C]
            #
            # For S = {i,j}:
            # |I|=0: contributes 0 (empty product is 1, no t_i t_j term)
            # |I|=1: single cycle C, coefficient of t_i t_j in χ_C, × weight 2
            # |I|=2: pair {C₁,C₂}, coefficient of t_i t_j in χ_{C₁} · χ_{C₂}, × weight 4

            type_a = 0
            for idx in range(nc):
                chi = all_cycles[idx][1]
                c_ij = chi.get(frozenset({i, j}), 0)
                type_a += 2 * c_ij

            type_b = 0
            for idx1 in range(nc):
                for idx2 in range(idx1+1, nc):
                    if idx2 in conflict[idx1]:
                        continue
                    chi1 = all_cycles[idx1][1]
                    chi2 = all_cycles[idx2][1]

                    # Product χ₁·χ₂ coefficient for {i,j}
                    # = Σ_{S₁∪S₂={i,j}} chi1[S₁] * chi2[S₂]
                    # Possible splits: ({i,j},∅), ({i},{j}), ({j},{i}), (∅,{i,j})
                    c = 0
                    c += chi1.get(frozenset({i,j}), 0) * chi2.get(frozenset(), 0)
                    c += chi1.get(frozenset({i}), 0) * chi2.get(frozenset({j}), 0)
                    c += chi1.get(frozenset({j}), 0) * chi2.get(frozenset({i}), 0)
                    c += chi1.get(frozenset(), 0) * chi2.get(frozenset({i,j}), 0)

                    type_b += 4 * c

            total = type_a + type_b
            if total != 0:
                print(f"({x1},{y1})-({x2},{y2})  {rel:<8} ({s1},{s2})  {total:>6} {type_a:>8} {type_b:>8}")

    # Detailed Type B analysis for disjoint tile pairs at n=6
    if n >= 6:
        print(f"\n--- TYPE B DETAILS (disjoint pairs, n={n}) ---")
        for i in range(m):
            for j in range(i+1, m):
                t1, t2 = tiles[i], tiles[j]
                x1, y1 = t1
                x2, y2 = t2
                if {x1,y1} & {x2,y2}:
                    continue  # shared vertex

                s1, s2 = x1-y1, x2-y2

                # Type B contributions
                pair_contribs = []
                for idx1 in range(nc):
                    for idx2 in range(idx1+1, nc):
                        if idx2 in conflict[idx1]:
                            continue
                        chi1 = all_cycles[idx1][1]
                        chi2 = all_cycles[idx2][1]

                        c = 0
                        c += chi1.get(frozenset({i,j}), 0) * chi2.get(frozenset(), 0)
                        c += chi1.get(frozenset({i}), 0) * chi2.get(frozenset({j}), 0)
                        c += chi1.get(frozenset({j}), 0) * chi2.get(frozenset({i}), 0)
                        c += chi1.get(frozenset(), 0) * chi2.get(frozenset({i,j}), 0)

                        if c != 0:
                            cyc1 = all_cycles[idx1][0]
                            cyc2 = all_cycles[idx2][0]
                            pair_contribs.append((cyc1, cyc2, 4*c))

                type_b = sum(c for _,_,c in pair_contribs)
                if type_b != 0 and len(pair_contribs) <= 5:
                    # Compute type_a too
                    type_a = 0
                    for idx in range(nc):
                        chi = all_cycles[idx][1]
                        type_a += 2 * chi.get(frozenset({i, j}), 0)

                    print(f"\n  ({x1},{y1})-({x2},{y2}) s=({s1},{s2}): "
                          f"TypeA={type_a}, TypeB={type_b}, total={type_a+type_b}")
                    for cyc1, cyc2, c in pair_contribs:
                        print(f"    pair {cyc1} + {cyc2}: contrib={c}")

if __name__ == "__main__":
    quadratic_anatomy(5)
    quadratic_anatomy(6, max_cycle_len=5)
