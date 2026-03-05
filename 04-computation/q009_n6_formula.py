"""
Q-009: Explicit ΔI formula at n=6.

At n=6, max independent set in Ω(T) has size 2 (pairs of VD 3-cycles).
So: I(Ω,2) = 1 + 2·|cycles| + 4·|{VD 3-cycle pairs}|

When flipping i→j to j→i:
  ΔI = 2·(|D| - |C|) + 4·Δ(#VD pairs)

Key derivation:
  #D3 - #C3 = sum_x (T[j][x]·T[x][i] - T[i][x]·T[x][j])
            = -sum_x s_x   where s_x = 1 - T[x][i] - T[j][x]

  Δ(#VD pairs) = sum_x s_x · c(B_x)
  where c(B_x) = 1 if V\{i,j,x} is a cyclic tournament

Therefore:
  ΔI = -2·sum_x s_x + 2·(#D5 - #C5) + 4·sum_x s_x · c(B_x)
     = 2·sum_x s_x · (2c(B_x) - 1) + 2·(#D5 - #C5)

Note: 2c(B_x)-1 = H(B_x)-2 where H(B_x) ∈ {1,3}.
So: ΔI = 2·sum_x s_x · (H(B_x) - 2) + 2·(#D5 - #C5)

We need to verify: ΔH = ΔI for all flips at n=6.

Author: opus-2026-03-05-S2
"""

import sys
sys.path.insert(0, '03-artifacts/code')
from tournament_lib import (
    all_tournaments, hamiltonian_path_count, find_odd_cycles
)
from itertools import permutations


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def is_cyclic_3(T, a, b, c):
    """Is the 3-vertex sub-tournament on {a,b,c} cyclic?"""
    return (T[a][b] and T[b][c] and T[c][a]) or (T[a][c] and T[c][b] and T[b][a])


def count_5cycles_using_arc(T, n, i, j):
    """Count directed 5-cycles in T that use arc i→j."""
    verts = list(range(n))
    count = 0
    # 5-cycles on 5 of 6 vertices, using i→j.
    # For each excluded vertex x (x != i,j):
    for x in verts:
        if x == i or x == j:
            continue
        subset = [u for u in verts if u != x]
        # Count directed 5-cycles on subset starting at i with i→j
        # Cycle: i→j→a→b→c→i (the 3 other verts after j must form path back to i)
        others = [u for u in subset if u != i and u != j]
        for perm in permutations(others):
            # Cycle: i→j→perm[0]→perm[1]→perm[2]→i
            if (T[j][perm[0]] and T[perm[0]][perm[1]] and
                T[perm[1]][perm[2]] and T[perm[2]][i]):
                count += 1
    # Also: i→j might not be at the start. The cycle can be
    # a→i→j→b→c→a for various positions.
    # Actually, a directed 5-cycle has 5 arcs. i→j is one of them.
    # The cycle visits all 5 vertices. Starting from i→j:
    # i→j→?→?→?→i (3 vertices in some order from j back to i)
    # This is what I computed above. But i→j could also appear at other
    # positions: e.g., a→i→j→b→c→a = same cycle, just different starting point.
    # Since I fix the representation as i→j→...→i, each cycle is counted once.
    # Actually no: the same cycle visiting 5 vertices can be written starting
    # from ANY of the 5 vertices. By starting from i, I get a unique representation.
    return count


def adj_count(T, i, j):
    """Count Ham paths of T where i immediately precedes j."""
    n = len(T)
    count = 0
    for perm in permutations(range(n)):
        valid = True
        for k in range(n - 1):
            if T[perm[k]][perm[k + 1]] != 1:
                valid = False
                break
        if valid:
            for k in range(n - 1):
                if perm[k] == i and perm[k + 1] == j:
                    count += 1
    return count


if __name__ == "__main__":
    n = 6
    print(f"Q-009: Verify explicit ΔI formula at n={n}")
    print("=" * 70)

    total = 0
    formula_match = 0
    delta_h_match = 0
    five_cycle_issue = 0

    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= 500:
            break
        for i in range(n):
            for j in range(i + 1, n):
                if T[i][j] == 0:
                    continue
                total += 1
                T2 = flip_arc(T, i, j)
                others = [x for x in range(n) if x != i and x != j]

                # Compute s_x, c(B_x), H(B_x) for each x
                sx_vals = []
                hbx_vals = []
                for x in others:
                    sx = 1 - T[x][i] - T[j][x]
                    bx = [u for u in range(n) if u not in {i, j, x}]
                    cx = is_cyclic_3(T, *bx)
                    hbx = 3 if cx else 1
                    sx_vals.append(sx)
                    hbx_vals.append(hbx)

                # 5-cycle counts
                d5 = count_5cycles_using_arc(T, n, i, j)
                c5 = count_5cycles_using_arc(T2, n, j, i)

                # Formula: ΔI = 2·sum_x s_x·(H(B_x)-2) + 2·(d5-c5)
                formula_delta_I = (2 * sum(s * (h - 2) for s, h in zip(sx_vals, hbx_vals))
                                   + 2 * (d5 - c5))

                # Ground truth
                ht = hamiltonian_path_count(T)
                ht2 = hamiltonian_path_count(T2)
                delta_H = ht - ht2

                # Also compute ΔI from independence polynomial
                cycles_T = find_odd_cycles(T)
                cycles_T2 = find_odd_cycles(T2)

                def count_vd_pairs(cycles):
                    count = 0
                    for a in range(len(cycles)):
                        for b in range(a + 1, len(cycles)):
                            if len(cycles[a]) == 3 and len(cycles[b]) == 3:
                                if not (set(cycles[a]) & set(cycles[b])):
                                    count += 1
                    return count

                I_T = 1 + 2 * len(cycles_T) + 4 * count_vd_pairs(cycles_T)
                I_T2 = 1 + 2 * len(cycles_T2) + 4 * count_vd_pairs(cycles_T2)
                delta_I_gt = I_T - I_T2

                if delta_H == delta_I_gt:
                    delta_h_match += 1

                if formula_delta_I == delta_I_gt:
                    formula_match += 1
                elif total <= 20:
                    print(f"  Formula mismatch T#{t_idx} flip {i}->{j}:")
                    print(f"    formula_ΔI={formula_delta_I}, gt_ΔI={delta_I_gt}")
                    print(f"    s_x={sx_vals}, H(B_x)={hbx_vals}")
                    print(f"    d5={d5}, c5={c5}")

    print(f"\nTotal: {total}")
    print(f"ΔH = ΔI (ground truth): {delta_h_match}/{total}")
    print(f"Formula ΔI = gt ΔI: {formula_match}/{total}")

    if formula_match == total:
        print("\n*** FORMULA VERIFIED! ***")
        print("ΔI = 2·Σ_x s_x·(H(B_x)-2) + 2·(#D5-#C5)")
        print("where s_x = 1 - T[x][i] - T[j][x], B_x = V\\{i,j,x}")
        print("\nThis reduces Q-009 to proving:")
        print("adj(i,j) - adj'(j,i) = 2·Σ_x s_x·(H(B_x)-2) + 2·(#D5-#C5)")

    # If formula works, let's also understand the 5-cycle term
    if formula_match == total or formula_match > total * 0.95:
        print(f"\n--- Analyzing 5-cycle contribution ---")
        d5_dist = {}
        c5_dist = {}
        for t_idx, T in enumerate(all_tournaments(n)):
            if t_idx >= 100:
                break
            for i in range(n):
                for j in range(i + 1, n):
                    if T[i][j] == 0:
                        continue
                    T2 = flip_arc(T, i, j)
                    d5 = count_5cycles_using_arc(T, n, i, j)
                    c5 = count_5cycles_using_arc(T2, n, j, i)
                    d5_dist[d5] = d5_dist.get(d5, 0) + 1
                    c5_dist[c5] = c5_dist.get(c5, 0) + 1

        print(f"  #destroyed 5-cycles distribution: {dict(sorted(d5_dist.items()))}")
        print(f"  #created 5-cycles distribution: {dict(sorted(c5_dist.items()))}")

        # Can we express #D5 in terms of local tournament structure?
        # D5 = # of 5-cycles on V\{x} using i→j, for each excluded x
        # = sum_x #{5-cycles on V\{x} that use i→j}
        # Wait, that's what I computed. But at n=6, V\{x} has 5 vertices.
        # A 5-cycle on 5 vertices using i→j = a directed Hamilton cycle of
        # T[V\{x}] using arc i→j.
        # = # Ham cycles of the 5-vertex tournament T[V\{x}] using i→j.
        #
        # A 5-vertex tournament has either 0 or 1 directed Hamilton cycle
        # (since it's either cyclic or not... actually no, a 5-tournament
        # can have multiple Hamilton cycles).
        #
        # Actually, #Ham_cycles(T[V\{x}]) using arc i→j =
        # #{paths j→...→i in T[V\{x}\{i,j}] (3 vertices)} that are valid.
        # = #{Ham paths from j to i in the 3-vertex tournament T[others\{x}]}

        print(f"\n  5-cycle decomposition:")
        print(f"  #D5 = sum_x #{'{'}Ham paths j→...→i in T[V\\{'{'}i,j,x{'}'}]{'}'})")
        for t_idx, T in enumerate(all_tournaments(n)):
            if t_idx >= 5:
                break
            for i in range(n):
                for j in range(i + 1, n):
                    if T[i][j] == 0:
                        continue
                    others = [x for x in range(n) if x != i and x != j]
                    d5_total = 0
                    d5_per_x = []
                    for x in others:
                        triple = [u for u in others if u != x]
                        # Count Ham paths from j to i in T[triple]
                        hp_ji = 0
                        for perm in permutations(triple):
                            if (T[j][perm[0]] and T[perm[0]][perm[1]] and
                                T[perm[1]][perm[2]] and T[perm[2]][i]):
                                hp_ji += 1
                        d5_per_x.append(hp_ji)
                        d5_total += hp_ji
                    if d5_total > 0:
                        print(f"  T#{t_idx} {i}->{j}: D5 per x: {dict(zip(others, d5_per_x))}, total={d5_total}")
                        break
                else:
                    continue
                break
