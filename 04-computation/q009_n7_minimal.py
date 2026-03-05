"""
Q-009: Minimal test of ΔI formula at n=7.

Skip the I(Ω,2) verification. Just test:
  ΔH = -2·Σ s_x·H(B_x) + 2·Σ_L(DL-CL)

All using DP, no find_odd_cycles.

Author: opus-2026-03-05-S2
"""

import random
from itertools import combinations


def hamiltonian_path_count_dp(T, verts):
    """Count Ham paths of T[verts] using bitmask DP."""
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


def count_cycles_using_arc(T, n, i, j):
    """Count directed odd cycles using arc i→j, by length.
    Cycle: i→j→(path through subset back to i).
    Uses DP on the internal vertices."""
    others = [u for u in range(n) if u != i and u != j]
    counts = {}
    for L in range(3, n + 1, 2):
        extra_needed = L - 2
        if extra_needed > len(others):
            break
        count = 0
        for extra in combinations(others, extra_needed):
            verts = list(extra)
            nv = len(verts)
            # Count paths j→verts[...]→i
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
    total = 0
    match = 0

    for n in [5, 6, 7]:
        n_trials = {5: 200, 6: 100, 7: 50}[n]
        n_match = 0
        n_total = 0
        all_verts = list(range(n))

        for _ in range(n_trials):
            T = random_tournament(n)
            i, j = random.sample(range(n), 2)
            if T[i][j] == 0:
                i, j = j, i

            T2 = flip_arc(T, i, j)
            others = [x for x in range(n) if x != i and x != j]
            n_total += 1
            total += 1

            # -2·Σ s_x·H(B_x)
            sx_hbx = 0
            for x in others:
                sx = 1 - T[x][i] - T[j][x]
                bx = [u for u in others if u != x]
                hbx = hamiltonian_path_count_dp(T, bx)
                sx_hbx += sx * hbx

            # Cycle destruction/creation terms
            dc_T = count_cycles_using_arc(T, n, i, j)
            dc_T2 = count_cycles_using_arc(T2, n, j, i)
            cycle_term = 0
            for L in range(3, n + 1, 2):
                cycle_term += dc_T.get(L, 0) - dc_T2.get(L, 0)

            formula = -2 * sx_hbx + 2 * cycle_term

            # Ground truth
            ht = hamiltonian_path_count_dp(T, all_verts)
            ht2 = hamiltonian_path_count_dp(T2, all_verts)
            delta_H = ht - ht2

            if formula == delta_H:
                n_match += 1
                match += 1
            else:
                print(f"  n={n} MISS: formula={formula}, ΔH={delta_H}, "
                      f"sx_hbx={sx_hbx}, cycle={cycle_term}")

        print(f"n={n}: {n_match}/{n_total}")

    print(f"\nTotal: {match}/{total}")
    if match == total:
        print("\n*** FORMULA VERIFIED FOR n=5,6,7! ***")
        print("ΔH = -2·Σ_x s_x·H(B_x) + 2·Σ_L(DL-CL)")
        print("This is the GENERAL arc-flip formula for OCF!")
