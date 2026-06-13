"""
Q-009: Derive and verify the CORRECT general formula for DeltaH.

Analysis shows:
  DeltaH = 2*Delta(alpha_1) + 4*Delta(alpha_2)

  Delta(alpha_1) = sum_L (DL - CL) = -sum(s_x) + sum_{L>=5}(DL-CL)

  Delta(alpha_2) = sum over L of [sum_{C of length L using i->j} alpha_1(complement(C))
                                  - sum_{C' of length L using j->i} alpha_1(complement(C'))]

At n<=7: the alpha_2 contribution from L=3 gives -sum(s_x * alpha_1(B_x)),
         and alpha_1(B_x) = #odd cycles in B_x = (H(B_x)-1)/2  (by OCF on smaller n).
         The L>=5 contribution to alpha_2 is 0 (complement too small).
         So: DeltaH = -2*sum(s_x*H(B_x)) + 2*sum_{L>=5}(DL-CL)  <-- verified formula

At n=8: the L=5 contribution to alpha_2 is NONZERO.
  Missing term = 4 * sum over 5-cycles C using i->j in T of alpha_1(V\V(C))
               - 4 * sum over 5-cycles C' using j->i in T' of alpha_1(V\V(C'))
  where alpha_1(S) = #directed odd cycles in T[S].

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


def count_odd_cycles_in_subtournament(T, verts):
    """Count directed odd cycles in T restricted to verts."""
    n = len(verts)
    count = 0
    for L in range(3, n + 1, 2):
        for subset in combinations(verts, L):
            v0 = subset[0]
            rest = list(subset[1:])
            for perm in permutations(rest):
                cycle = (v0,) + perm
                valid = True
                for k in range(L):
                    if T[cycle[k]][cycle[(k+1) % L]] != 1:
                        valid = False
                        break
                if valid:
                    count += 1
    return count


def count_cycles_using_arc_with_complement_info(T, n, i, j, min_length=3):
    """For each L-cycle using arc i->j, compute alpha_1(complement).
    Returns dict: L -> sum of alpha_1(complement) over all L-cycles using i->j."""
    others = [u for u in range(n) if u != i and u != j]
    all_verts = set(range(n))
    results = {}
    for L in range(min_length, n + 1, 2):
        extra_needed = L - 2
        if extra_needed > len(others):
            break
        total_alpha1 = 0
        total_count = 0
        for extra in combinations(others, extra_needed):
            verts = list(extra)
            nv = len(verts)
            # Count paths j -> verts -> i using DP
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
            num_cycles_on_this_set = 0
            for v in range(nv):
                if dp[full][v] > 0 and T[verts[v]][i]:
                    num_cycles_on_this_set += dp[full][v]

            if num_cycles_on_this_set > 0:
                # Complement vertices
                cycle_verts = set(extra) | {i, j}
                complement = sorted(all_verts - cycle_verts)
                alpha1_comp = count_odd_cycles_in_subtournament(T, complement)
                total_alpha1 += num_cycles_on_this_set * alpha1_comp
                total_count += num_cycles_on_this_set

        results[L] = (total_count, total_alpha1)
    return results


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

    match_old = 0
    match_new = 0
    total = 0

    for trial in range(15):
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

        # Old formula: -2*sum(s_x*H(B_x)) + 2*sum_{L>=5}(DL-CL)
        sx_hbx = 0
        for x in others:
            sx = 1 - T[x][i] - T[j][x]
            bx = [u for u in others if u != x]
            hbx = hamiltonian_path_count_dp(T, bx)
            sx_hbx += sx * hbx

        # Cycle info with complement alpha_1
        info_T = count_cycles_using_arc_with_complement_info(T, n, i, j)
        info_T2 = count_cycles_using_arc_with_complement_info(T2, n, j, i)

        # DL, CL for each length
        long_cycle_term = 0
        for L in range(5, n+1, 2):
            DL = info_T.get(L, (0,0))[0]
            CL = info_T2.get(L, (0,0))[0]
            long_cycle_term += DL - CL

        old_formula = -2 * sx_hbx + 2 * long_cycle_term

        # New correction: 4 * sum_{L>=5} [sum_{C_L} alpha_1(comp(C_L)) in T - same in T']
        correction = 0
        for L in range(5, n+1, 2):
            alpha1_sum_T = info_T.get(L, (0,0))[1]
            alpha1_sum_T2 = info_T2.get(L, (0,0))[1]
            correction += alpha1_sum_T - alpha1_sum_T2

        new_formula = old_formula + 4 * correction

        total += 1
        if old_formula == delta_H:
            match_old += 1
        if new_formula == delta_H:
            match_new += 1
        else:
            print(f"  Trial {trial}: DH={delta_H}, old={old_formula}, new={new_formula}, "
                  f"correction={4*correction}, residual_new={delta_H - new_formula}")

        if (trial + 1) % 5 == 0:
            print(f"Progress: {trial+1}, old={match_old}/{total}, new={match_new}/{total}")

    print(f"\nOld formula: {match_old}/{total}")
    print(f"New formula: {match_new}/{total}")
