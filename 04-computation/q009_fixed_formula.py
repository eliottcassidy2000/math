"""
Q-009: CORRECTED ΔI formula at n=6.

ΔI = I(Ω(T),2) - I(Ω(T'),2)
   = 2·(#cycles_T - #cycles_T') + 4·(#VD_T - #VD_T')

With correct signs:
  #cycles_T - #cycles_T' = -Σ s_x + (D5 - C5)
  #VD_T - #VD_T' = -Σ s_x · c(B_x)

So:
  ΔI = -2·Σ s_x · (1 + 2·c(B_x)) + 2·(D5 - C5)

Since H(B_x) = 1 + 2·c(B_x) for a 3-vertex tournament:

  *** ΔI = -2·Σ_x s_x · H(T[V\{i,j,x}]) + 2·(D5 - C5) ***

And OCF (ΔH = ΔI) reduces to proving:
  adj(i,j) - adj'(j,i) = -2·Σ_x s_x · H(B_x) + 2·(D5 - C5)

Author: opus-2026-03-05-S2
"""

import sys
sys.path.insert(0, '03-artifacts/code')
from tournament_lib import (
    all_tournaments, hamiltonian_path_count, find_odd_cycles,
    conflict_graph, independence_poly_at
)
from itertools import permutations


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def is_cyclic_3(T, a, b, c):
    return (T[a][b] and T[b][c] and T[c][a]) or (T[a][c] and T[c][b] and T[b][a])


def sub_H(T, verts):
    """H of sub-tournament T[verts]."""
    verts = sorted(verts)
    m = len(verts)
    if m <= 1:
        return 1
    count = 0
    for perm in permutations(verts):
        valid = True
        for k in range(m - 1):
            if T[perm[k]][perm[k + 1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count


def count_5cycles_ij(T, n, i, j):
    """Count 5-cycles in T using arc i→j."""
    count = 0
    for x in range(n):
        if x == i or x == j:
            continue
        others = [u for u in range(n) if u not in {i, j, x}]
        for perm in permutations(others):
            if (T[j][perm[0]] and T[perm[0]][perm[1]] and
                T[perm[1]][perm[2]] and T[perm[2]][i]):
                count += 1
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
    print(f"Q-009: Corrected ΔI formula at n={n}")
    print("=" * 70)

    total = 0
    formula_ok = 0
    ocf_ok = 0

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

                # Compute formula ΔI
                sx_hbx_sum = 0
                for x in others:
                    sx = 1 - T[x][i] - T[j][x]
                    bx = [u for u in others if u != x]
                    hbx = sub_H(T, bx)
                    sx_hbx_sum += sx * hbx

                d5 = count_5cycles_ij(T, n, i, j)
                c5 = count_5cycles_ij(T2, n, j, i)

                formula_dI = -2 * sx_hbx_sum + 2 * (d5 - c5)

                # Ground truth
                ht = hamiltonian_path_count(T)
                ht2 = hamiltonian_path_count(T2)
                delta_H = ht - ht2

                # Check ΔI formula against ground truth I computation
                cycles_T = find_odd_cycles(T)
                cycles_T2 = find_odd_cycles(T2)
                if cycles_T:
                    cg = conflict_graph(cycles_T)
                    I_T = independence_poly_at(cg, 2)
                else:
                    I_T = 1
                if cycles_T2:
                    cg2 = conflict_graph(cycles_T2)
                    I_T2 = independence_poly_at(cg2, 2)
                else:
                    I_T2 = 1
                gt_dI = I_T - I_T2

                if formula_dI == gt_dI:
                    formula_ok += 1
                elif total <= 5:
                    print(f"  FORMULA MISS: T#{t_idx} {i}->{j}: "
                          f"formula={formula_dI}, gt={gt_dI}")

                if delta_H == gt_dI:
                    ocf_ok += 1

    print(f"\nTotal: {total}")
    print(f"Formula ΔI matches ground truth: {formula_ok}/{total}")
    print(f"ΔH = ΔI (OCF consistency): {ocf_ok}/{total}")

    if formula_ok == total:
        print("\n*** FORMULA VERIFIED! ***")
        print("ΔI = -2·Σ_x s_x·H(B_x) + 2·(D5-C5)")
        print(f"\nSince ΔH = ΔI = {ocf_ok}/{total}, OCF reduces to proving:")
        print("adj(i,j) - adj'(j,i) = -2·Σ_x s_x·H(B_x) + 2·(D5-C5)")
        print("\nwhere s_x = 1 - T[x][i] - T[j][x] ∈ {-1, 0, 1}")
        print("and B_x = V\\{i,j,x} (3-vertex sub-tournament)")
        print("and D5, C5 = destroyed/created 5-cycles using arc i→j / j→i")

        # Explore: can D5, C5 also be expressed using s_x and sub-tournament H?
        print("\n--- Can we express D5-C5 in terms of s_x and sub-tournament H? ---")
        # D5 = Σ_x #{5-cycles on V\{x} using i→j}
        # = Σ_x #{Ham paths from j to i in T[V\{i,j,x}]}
        # = Σ_x h_{j→i}(B_x)
        # where h_{j→i}(B_x) = # Ham paths from j to i in T[B_x]... no,
        # the 5-cycle is on V\{x}, including i and j. The cycle goes
        # i→j→(path through B_x)→i. So the path through B_x goes j→...→i
        # which is a path from j to i through the 3 vertices of B_x = V\{i,j,x}.
        # Wait, B_x has 3 vertices (none being i or j).
        # The path j→p0→p1→p2→i goes from j to i through B_x = {p0,p1,p2}.
        # But j and i are NOT in B_x. The arcs used: j→p0 (from T, not restricted
        # to B_x), p0→p1 (within B_x), p1→p2 (within B_x), p2→i (from T).
        #
        # So D5 = Σ_x #{permutations (p0,p1,p2) of B_x such that
        #          T[j][p0]=1 AND T[p0][p1]=1 AND T[p1][p2]=1 AND T[p2][i]=1}

        # Similarly: C5 = Σ_x #{perms of B_x: T'[i][p0] AND T[p0][p1] AND T[p1][p2] AND T'[p2][j]}
        # Since T'=T except T'[j][i]=1, T'[i][j]=0:
        # T'[i][p0] = T[i][p0] (p0 not i or j), T'[p2][j] = T[p2][j]
        # So C5 = Σ_x #{perms: T[i][p0] AND T[p0][p1] AND T[p1][p2] AND T[p2][j]}

        # Define: for a permutation (p0,p1,p2) of B_x:
        # f_D(x,perm) = T[j][p0] · T[p2][i] · Π_{internal arcs}
        # f_C(x,perm) = T[i][p0] · T[p2][j] · Π_{internal arcs}

        # Note: Π_{internal arcs} = T[p0][p1] · T[p1][p2] is the same for D and C.
        # So:
        # D5 - C5 = Σ_x Σ_{perm} Π_{internal} · (T[j][p0]·T[p2][i] - T[i][p0]·T[p2][j])

        # For a VALID Ham path p0→p1→p2 in B_x:
        # D5 - C5 = Σ_x Σ_{valid path P in B_x} (T[j][P[0]]·T[P[2]][i] - T[i][P[0]]·T[P[2]][j])

        print("D5-C5 = Σ_x Σ_{path P in B_x} (T[j][P[0]]·T[P[2]][i] - T[i][P[0]]·T[P[2]][j])")
        print("This is a sum over Ham paths of B_x weighted by boundary arcs.")

        # Full OCF equation:
        # adj(i,j) - adj'(j,i) = -2·Σ_x s_x·H(B_x)
        #   + 2·Σ_x Σ_{path P in B_x} (T[j][P[0]]·T[P[2]][i] - T[i][P[0]]·T[P[2]][j])

        # SIMPLIFY: adj(i,j) counts 6-vertex paths ...→i→j→...
        # Can split as: before-i part (ends at i) and after-j part (starts at j).
        # The 4 other vertices split between before and after.
        # For each x (excluded from the non-{i,j} split):
        #   Wait, this is different...

        # Actually, let me just verify the full equation numerically
        print("\n--- Full equation verification ---")
        eq_ok = 0
        for t_idx, T in enumerate(all_tournaments(n)):
            if t_idx >= 500:
                break
            for i in range(n):
                for j in range(i + 1, n):
                    if T[i][j] == 0:
                        continue
                    T2 = flip_arc(T, i, j)
                    others = [x for x in range(n) if x != i and x != j]

                    lhs = adj_count(T, i, j) - adj_count(T2, j, i)

                    rhs = 0
                    for x in others:
                        sx = 1 - T[x][i] - T[j][x]
                        bx = [u for u in others if u != x]
                        hbx = sub_H(T, bx)
                        rhs -= 2 * sx * hbx

                        # 5-cycle contribution for this x
                        for perm in permutations(bx):
                            internal = T[perm[0]][perm[1]] * T[perm[1]][perm[2]]
                            if internal:
                                rhs += 2 * (T[j][perm[0]] * T[perm[2]][i] -
                                            T[i][perm[0]] * T[perm[2]][j])

                    if lhs == rhs:
                        eq_ok += 1

        print(f"Full equation adj(i,j)-adj'(j,i) = RHS: {eq_ok}/{total}")
