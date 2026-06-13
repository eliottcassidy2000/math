"""
Full decomposition of arc reversal effects on Claim A.

When we flip arc i->j (i,j != v), the three quantities change:
  delta_H = H(T) - H(T')
  delta_Hv = H(T-v) - H(T'-v)
  delta_mu = sum mu_T(C) - sum mu_{T'}(C)

Claim A says H(T) - H(T-v) = 2*sum mu(C), so:
  delta_H - delta_Hv = 2*delta_mu

Let's verify this and find formulas for each delta.

Known: H(T) - H(T') = adj_T(i,j) - adj_{T'}(j,i) (adjacency formula)
Similarly: H(T-v) - H(T'-v) = adj_{T-v}(i,j) - adj_{T'-v}(j,i)

So delta_H - delta_Hv = [adj_T(i,j) - adj_{T'}(j,i)] - [adj_{T-v}(i,j) - adj_{T'-v}(j,i)]

And this should equal 2*delta_mu.

Author: opus-2026-03-05-S2
"""

import sys
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (
    all_tournaments, hamiltonian_path_count, delete_vertex,
    find_odd_cycles, mu
)
from itertools import permutations


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def adj_count_dp(T, i, j):
    """Count Ham paths of T where i immediately precedes j. Uses DP."""
    n = len(T)
    if n <= 1 or T[i][j] != 1:
        return 0
    # dp[mask][last] = # of paths through vertices in mask ending at last
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    count = 0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
                    # If this completes the path and v=i, u=j, count it
                    # Actually we need to count all paths where i->j appears
    # Recount: for each complete path, check if i immediately before j
    count = 0
    for perm in permutations(range(n)):
        valid = True
        for k in range(n - 1):
            if T[perm[k]][perm[k+1]] != 1:
                valid = False
                break
        if valid:
            for k in range(n - 1):
                if perm[k] == i and perm[k+1] == j:
                    count += 1
    return count


def compute_sum_mu(T, v):
    """Compute sum_{C through v} mu(C)."""
    all_cycles = find_odd_cycles(T)
    cycles_v = [c for c in all_cycles if v in set(c)]
    if not cycles_v:
        return 0
    Tv, old_labels = delete_vertex(T, v)
    tv_cycles = find_odd_cycles(Tv)
    cache = (Tv, old_labels, tv_cycles)
    return sum(mu(T, v, c, _tv_cache=cache) for c in cycles_v)


if __name__ == "__main__":
    n = 5
    print(f"Arc-reversal decomposition at n={n}")
    print("="*70)

    # Verify: H(T) - H(T') = adj_T(i,j) - adj_{T'}(j,i)
    adj_formula_ok = 0
    decomp_ok = 0
    total = 0

    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= 200:
            break
        for v in range(n):
            others = [u for u in range(n) if u != v]
            for i in others:
                for j in others:
                    if i == j or T[i][j] == 0:
                        continue
                    total += 1
                    T2 = flip_arc(T, i, j)

                    ht = hamiltonian_path_count(T)
                    ht2 = hamiltonian_path_count(T2)
                    delta_H = ht - ht2

                    # Adjacency formula for full tournament
                    adj_ij = adj_count_dp(T, i, j)
                    adj_ji2 = adj_count_dp(T2, j, i)
                    if delta_H == adj_ij - adj_ji2:
                        adj_formula_ok += 1

                    # For T-v
                    Tv, _ = delete_vertex(T, v)
                    Tv2, _ = delete_vertex(T2, v)
                    htv = hamiltonian_path_count(Tv)
                    htv2 = hamiltonian_path_count(Tv2)
                    delta_Hv = htv - htv2

                    # mu sums
                    mu_sum = compute_sum_mu(T, v)
                    mu_sum2 = compute_sum_mu(T2, v)
                    delta_mu = mu_sum - mu_sum2

                    # Check: delta_H - delta_Hv = 2 * delta_mu
                    if delta_H - delta_Hv == 2 * delta_mu:
                        decomp_ok += 1

    print(f"Total arc flips: {total}")
    print(f"H(T)-H(T') = adj(i,j) - adj'(j,i): {adj_formula_ok}/{total}")
    print(f"delta_H - delta_Hv = 2*delta_mu: {decomp_ok}/{total}")

    # Now explore: what IS the relationship between adj and mu?
    print(f"\n--- Decomposing delta_mu ---")
    # delta_mu = sum mu_T(C) - sum mu_{T'}(C)
    # Split into:
    #   cycles destroyed by flip (exist in T, not T')
    #   cycles created by flip (exist in T', not T)
    #   cycles that persist but mu changes

    n = 5
    destroy_total = 0
    create_total = 0
    persist_change = 0
    persist_same = 0

    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= 100:
            break
        for v in range(n):
            others = [u for u in range(n) if u != v]
            for i in others:
                for j in others:
                    if i == j or T[i][j] == 0:
                        continue
                    T2 = flip_arc(T, i, j)
                    cycles_T = set(tuple(c) for c in find_odd_cycles(T)
                                   if v in set(c))
                    cycles_T2 = set(tuple(c) for c in find_odd_cycles(T2)
                                    if v in set(c))

                    destroyed = cycles_T - cycles_T2
                    created = cycles_T2 - cycles_T
                    persisted = cycles_T & cycles_T2

                    destroy_total += len(destroyed)
                    create_total += len(created)

                    for c in persisted:
                        Tv, ol = delete_vertex(T, v)
                        tc = find_odd_cycles(Tv)
                        mu_T = mu(T, v, c, (Tv, ol, tc))
                        Tv2, ol2 = delete_vertex(T2, v)
                        tc2 = find_odd_cycles(Tv2)
                        mu_T2 = mu(T2, v, c, (Tv2, ol2, tc2))
                        if mu_T != mu_T2:
                            persist_change += 1
                        else:
                            persist_same += 1

    print(f"Cycles destroyed by flip: {destroy_total}")
    print(f"Cycles created by flip: {create_total}")
    print(f"Persisting cycles, mu changed: {persist_change}")
    print(f"Persisting cycles, mu same: {persist_same}")

    # At n=5 all mu=1, so "mu changed" should be 0
    if persist_change == 0:
        print("(All mu=1 at n=5, so no mu changes for persisting cycles)")

    print("\nConclusion: delta_mu = sum_destroyed mu_T - sum_created mu_T'"
          " + sum_persist (mu_T - mu_T')")
    print("At n=5: delta_mu = #destroyed - #created (since all mu=1)")
    print("The arc-reversal proof needs: delta_H - delta_Hv = "
          "2*(#destroyed - #created)")
