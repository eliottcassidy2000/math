"""
Q-009: CORRECT general formula for DeltaH, valid for ALL n.

DeltaH = DeltaI(Omega, 2) = sum_{k>=1} 2^k * Delta(alpha_k)

where alpha_k = # independent sets of size k in the conflict graph Omega.

Delta(alpha_1) = sum_L (DL - CL)  [cycles destroyed minus created]

Delta(alpha_2) = sum_L [sum_{C using i->j in T} alpha_1(complement(C))
                       - sum_{C' using j->i in T'} alpha_1(complement(C'))]

For k=1,2 this gives:
  DeltaH = 2*sum_L(DL-CL) + 4*[sum_L sum_C alpha_1(comp(C)) - ...]

Practically, at n<=8: alpha_k = 0 for k>=3, so:
  DeltaH = 2*Delta(alpha_1) + 4*Delta(alpha_2)

The COLLAPSED form (using OCF on sub-tournaments):
  DeltaH = -2*sum_x s_x*H(B_x) + 2*sum_{L>=5}(DL-CL)
           + 4*sum_{L>=5}[sum_C alpha_1(comp(C)) - sum_C' alpha_1(comp(C'))]

where the last term is nonzero only when n >= 8.

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


def count_odd_cycles_in_sub(T, verts):
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
                    if T[cycle[k]][cycle[(k + 1) % L]] != 1:
                        valid = False
                        break
                if valid:
                    count += 1
    return count


def count_vd_pairs_in_sub(T, verts):
    """Count VD odd-cycle pairs in T restricted to verts."""
    # First find all cycles
    cycles = []
    n = len(verts)
    for L in range(3, n + 1, 2):
        for subset in combinations(verts, L):
            v0 = subset[0]
            rest = list(subset[1:])
            for perm in permutations(rest):
                cycle = (v0,) + perm
                valid = True
                for k in range(L):
                    if T[cycle[k]][cycle[(k + 1) % L]] != 1:
                        valid = False
                        break
                if valid:
                    cycles.append(frozenset(subset))
    # Count VD pairs
    count = 0
    for a in range(len(cycles)):
        for b in range(a + 1, len(cycles)):
            if not (cycles[a] & cycles[b]):
                count += 1
    return count


def cycles_using_arc_info(T, n, i, j):
    """For each L, find cycles using arc i->j. Return:
    - count: number of such cycles
    - alpha1_sum: sum of alpha_1(complement(C)) over such cycles
    - alpha2_sum: sum of alpha_2(complement(C)) over such cycles
    """
    others = [u for u in range(n) if u != i and u != j]
    all_verts = set(range(n))
    results = {}
    for L in range(3, n + 1, 2):
        extra_needed = L - 2
        if extra_needed > len(others):
            break
        count = 0
        a1_sum = 0
        a2_sum = 0
        for extra in combinations(others, extra_needed):
            verts = list(extra)
            nv = len(verts)
            # Count paths j -> verts -> i
            if nv == 0:
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
            n_cycles = sum(dp[full][v] for v in range(nv) if T[verts[v]][i])
            if n_cycles > 0:
                comp = sorted(all_verts - set(extra) - {i, j})
                a1 = count_odd_cycles_in_sub(T, comp)
                a2 = count_vd_pairs_in_sub(T, comp)
                count += n_cycles
                a1_sum += n_cycles * a1
                a2_sum += n_cycles * a2
        results[L] = (count, a1_sum, a2_sum)
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

    for n in [4, 5, 6, 7, 8]:
        n_trials = {4: 200, 5: 100, 6: 50, 7: 30, 8: 15}[n]
        n_match = 0
        all_verts = list(range(n))

        for trial in range(n_trials):
            T = random_tournament(n)
            i, j = random.sample(range(n), 2)
            if T[i][j] == 0:
                i, j = j, i

            T2 = flip_arc(T, i, j)

            # Ground truth
            ht = hamiltonian_path_count_dp(T, all_verts)
            ht2 = hamiltonian_path_count_dp(T2, all_verts)
            delta_H = ht - ht2

            # Full formula: DeltaH = sum_{k>=1} 2^k * Delta(alpha_k)
            info_T = cycles_using_arc_info(T, n, i, j)
            info_T2 = cycles_using_arc_info(T2, n, j, i)

            # Delta(alpha_1) = sum_L (DL - CL)
            delta_a1 = sum(info_T.get(L, (0,0,0))[0] - info_T2.get(L, (0,0,0))[0]
                          for L in range(3, n+1, 2))

            # Delta(alpha_2) = sum_L [sum_C alpha_1(comp(C)) - sum_C' alpha_1(comp(C'))]
            delta_a2 = sum(info_T.get(L, (0,0,0))[1] - info_T2.get(L, (0,0,0))[1]
                          for L in range(3, n+1, 2))

            # Delta(alpha_3) = sum_L [sum_C alpha_2(comp(C)) - sum_C' alpha_2(comp(C'))]
            delta_a3 = sum(info_T.get(L, (0,0,0))[2] - info_T2.get(L, (0,0,0))[2]
                          for L in range(3, n+1, 2))

            formula = 2 * delta_a1 + 4 * delta_a2 + 8 * delta_a3

            if formula == delta_H:
                n_match += 1
            else:
                print(f"  n={n} trial={trial}: DH={delta_H}, formula={formula}, "
                      f"da1={delta_a1}, da2={delta_a2}, da3={delta_a3}")

        print(f"n={n}: {n_match}/{n_trials}")

    print("\nFormula: DH = 2*Da1 + 4*Da2 + 8*Da3 + ...")
    print("= DeltaI(Omega, 2)")
