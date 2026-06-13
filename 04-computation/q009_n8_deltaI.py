"""
Q-009: Compute DeltaI from conflict graph at n=8 and decompose.

ΔI = I(Ω(T),2) - I(Ω(T'),2)
   = 2·Δ(#cycles) + 4·Δ(#VD pairs)

The old formula gives: -2·Σ s_x·H(B_x) + 2·Σ_{L≥5}(DL-CL)

At n=6, this was derived from:
  2·Δ(#cycles) = 2·(-Σ s_x + D5-C5)
  4·Δ(#VD 3-3 pairs) = -4·Σ s_x·c(B_x)
  Combining: -2·Σ s_x·(1+2c(B_x)) + 2·(D5-C5) = -2·Σ s_x·H(B_x) + 2·(D5-C5)

At n=8, the VD pairs include 3-5 pairs, which are NOT captured.
So the residual = 4·Δ(#VD pairs with at least one ≥5 cycle).

But does the formula -2·Σ s_x·H(B_x) at n=8 already implicitly capture
more than just 3-3 VD pairs? Let me check.

Author: opus-2026-03-05-S2
"""

import random
from itertools import combinations, permutations


def hamiltonian_path_count_dp(T, verts):
    n = len(verts)
    if n <= 1:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for k in range(n):
        dp[1 << k][k] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[verts[v]][verts[u]]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


def find_all_directed_odd_cycles(T, n):
    cycles = []
    for L in range(3, n + 1, 2):
        for verts in combinations(range(n), L):
            v0 = verts[0]
            rest = list(verts[1:])
            for perm in permutations(rest):
                cycle = (v0,) + perm
                valid = True
                for k in range(L):
                    if T[cycle[k]][cycle[(k+1) % L]] != 1:
                        valid = False
                        break
                if valid:
                    cycles.append(cycle)
    return cycles


def count_vd_pairs(cycles):
    nc = len(cycles)
    vertex_sets = [frozenset(c) for c in cycles]
    count = 0
    for a in range(nc):
        for b in range(a + 1, nc):
            if not (vertex_sets[a] & vertex_sets[b]):
                count += 1
    return count


def count_vd_pairs_by_type(cycles):
    from collections import Counter
    nc = len(cycles)
    vertex_sets = [frozenset(c) for c in cycles]
    type_counts = Counter()
    for a in range(nc):
        for b in range(a + 1, nc):
            if not (vertex_sets[a] & vertex_sets[b]):
                la, lb = len(cycles[a]), len(cycles[b])
                key = tuple(sorted([la, lb]))
                type_counts[key] += 1
    return type_counts


def count_cycles_using_arc_by_length(T, n, i, j, min_length=3):
    others = [u for u in range(n) if u != i and u != j]
    counts = {}
    for L in range(min_length, n + 1, 2):
        extra_needed = L - 2
        if extra_needed > len(others):
            break
        count = 0
        for extra in combinations(others, extra_needed):
            verts = list(extra)
            nv = len(verts)
            if nv == 0:
                # L=2, not odd
                continue
            dp = [[0] * nv for _ in range(1 << nv)]
            for k in range(nv):
                if T[j][verts[k]]:
                    dp[1 << k][k] = 1
            full = (1 << nv) - 1
            for mask in range(1, 1 << nv):
                for v in range(nv):
                    if not (mask & (1 << v)) or dp[mask][v] == 0:
                        continue
                    for u in range(nv):
                        if mask & (1 << u):
                            continue
                        if T[verts[v]][verts[u]]:
                            dp[mask | (1 << u)][u] += dp[mask][v]
            for v in range(nv):
                if dp[full][v] > 0 and T[verts[v]][i]:
                    count += dp[full][v]
        counts[L] = count
    return counts


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def random_tournament(n):
    T = [[0]*n for _ in range(n)]
    for a in range(n):
        for b in range(a+1, n):
            if random.random() < 0.5:
                T[a][b] = 1
            else:
                T[b][a] = 1
    return T


if __name__ == "__main__":
    random.seed(42)
    n = 8
    all_verts = list(range(n))

    for trial in range(8):
        T = random_tournament(n)
        i, j = random.sample(range(n), 2)
        if T[i][j] == 0:
            i, j = j, i

        T2 = flip_arc(T, i, j)
        others = [x for x in range(n) if x != i and x != j]

        # Ground truth
        ht = hamiltonian_path_count_dp(T, all_verts)
        ht2 = hamiltonian_path_count_dp(T2, all_verts)
        delta_H = ht - ht2

        # Old formula
        sx_hbx = 0
        for x in others:
            sx = 1 - T[x][i] - T[j][x]
            bx = [u for u in others if u != x]
            hbx = hamiltonian_path_count_dp(T, bx)
            sx_hbx += sx * hbx

        dc_T = count_cycles_using_arc_by_length(T, n, i, j, min_length=5)
        dc_T2 = count_cycles_using_arc_by_length(T2, n, j, i, min_length=5)
        long_cycle_term = sum(dc_T.get(L, 0) - dc_T2.get(L, 0)
                              for L in range(5, n+1, 2))
        old_formula = -2 * sx_hbx + 2 * long_cycle_term

        # Conflict graph approach: DeltaI = 2*Delta(#cycles) + 4*Delta(#VD pairs)
        cycles_T = find_all_directed_odd_cycles(T, n)
        cycles_T2 = find_all_directed_odd_cycles(T2, n)
        delta_cycles = len(cycles_T) - len(cycles_T2)
        vd_T = count_vd_pairs(cycles_T)
        vd_T2 = count_vd_pairs(cycles_T2)
        delta_vd = vd_T - vd_T2
        delta_I = 2 * delta_cycles + 4 * delta_vd

        # Decompose delta_cycles
        dc3_T = count_cycles_using_arc_by_length(T, n, i, j, min_length=3)
        dc3_T2 = count_cycles_using_arc_by_length(T2, n, j, i, min_length=3)

        # Check: all destroyed/created cycles use arc i->j or j->i
        # So delta_cycles = sum over L of (dc3_T[L] - dc3_T2[L])
        delta_cycles_from_arc = sum(dc3_T.get(L, 0) - dc3_T2.get(L, 0)
                                    for L in range(3, n+1, 2))
        # And D3-C3 should = -sum(s_x)
        sum_sx = sum(1 - T[x][i] - T[j][x] for x in others)
        d3_c3 = dc3_T.get(3, 0) - dc3_T2.get(3, 0)

        # VD pairs by type
        vd_type_T = count_vd_pairs_by_type(cycles_T)
        vd_type_T2 = count_vd_pairs_by_type(cycles_T2)

        residual = delta_H - old_formula

        print(f"\nTrial {trial}: flip {i}->{j}")
        print(f"  DH={delta_H}, old_formula={old_formula}, residual={residual}")
        print(f"  DeltaI_from_conflict={delta_I}, match_DH={delta_I == delta_H}")
        print(f"  delta_cycles={delta_cycles} (from_arc={delta_cycles_from_arc})")
        print(f"  D3-C3={d3_c3}, -sum_sx={-sum_sx}")
        print(f"  delta_vd={delta_vd}")
        print(f"  VD_T: {dict(vd_type_T)}")
        print(f"  VD_T': {dict(vd_type_T2)}")

        # Check what the formula components give vs actual
        # The formula's "cycle" contribution: 2*sum(DL-CL) for L>=5
        # The actual cycle contribution to DeltaI: 2*delta_cycles = 2*(D3-C3 + D5-C5 + D7-C7)
        # The formula also encodes 2*(D3-C3) inside -2*sum(s_x) (the "1" in H(B_x))
        # So the formula's total cycle contribution = -2*sum(s_x) + 2*(D5-C5+D7-C7) = 2*delta_cycles ✓

        # The formula's "VD" contribution: the remaining part of -2*sx_hbx
        # -2*sx_hbx = -2*sum(s_x) - 2*sum(s_x*(H(B_x)-1))
        # H(B_x) - 1 = H(B_x) - 1. For 3-vertex tournament, H-1 = 2*c(B_x)
        # For 5-vertex tournament, H-1 = ???
        formula_vd_part = -2 * (sx_hbx - sum_sx)  # = -2*sum(s_x*(H(B_x)-1))
        actual_vd_part = 4 * delta_vd

        print(f"  formula_vd_part={formula_vd_part}, actual_4*dvd={actual_vd_part}, diff={actual_vd_part - formula_vd_part}")

        # Decompose further: what does -2*sum(s_x*(H(B_x)-1)) capture?
        # At n=6: B_x has 3 vertices, H(B_x)-1 = 2*c(B_x), so this = -4*sum(s_x*c(B_x))
        # which = 4*delta(#VD 3-3 pairs). QED for n=6.
        # At n=8: B_x has 5 vertices. H(B_x)-1 is more complex.
        # What IS H(B_x)-1 for a 5-vertex tournament?
        for x in others[:2]:
            sx = 1 - T[x][i] - T[j][x]
            bx = [u for u in others if u != x]
            hbx = hamiltonian_path_count_dp(T, bx)
            # Count odd cycles in B_x
            bx_cycles = find_all_directed_odd_cycles(
                T, n  # Use full T but only look at B_x vertices
            )
            # Actually need sub-tournament cycles
            bx_set = set(bx)
            bx_cycles = [c for c in cycles_T if set(c) <= bx_set]
            bx_vd = count_vd_pairs(bx_cycles)
            print(f"    x={x}: sx={sx}, H(B_x)={hbx}, H-1={hbx-1}, "
                  f"#cycles_in_Bx={len(bx_cycles)}, #VD_in_Bx={bx_vd}, "
                  f"I(Omega_Bx,2)-1={2*len(bx_cycles)+4*bx_vd}")
