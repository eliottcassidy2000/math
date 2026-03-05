"""
Q-009: Test general arc-flip formula at n=8.

ΔH = -2·Σ_x s_x·H(B_x) + 2·Σ_{L≥5,odd}(DL-CL)

n=8 is tractable with bitmask DP (2^8=256 states).
At n=8, max independent set size in Ω is still 2 (need n≥9 for 3 VD 3-cycles).

Author: opus-2026-03-05-S2
"""

import random
from itertools import combinations


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


def count_cycles_using_arc_by_length(T, n, i, j, min_length=5):
    """Count directed odd cycles of length >= min_length using arc i->j."""
    others = [u for u in range(n) if u != i and u != j]
    total = 0
    for L in range(min_length, n + 1, 2):
        extra_needed = L - 2
        if extra_needed > len(others):
            break
        for extra in combinations(others, extra_needed):
            verts = list(extra)
            nv = len(verts)
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
                    total += dp[full][v]
    return total


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

    for n in [8]:
        n_trials = 30
        n_match = 0
        n_total = 0
        all_verts = list(range(n))

        for trial in range(n_trials):
            T = random_tournament(n)
            i, j = random.sample(range(n), 2)
            if T[i][j] == 0:
                i, j = j, i

            T2 = flip_arc(T, i, j)
            others = [x for x in range(n) if x != i and x != j]
            n_total += 1

            # -2·Σ s_x · H(B_x)
            sx_hbx = 0
            for x in others:
                sx = 1 - T[x][i] - T[j][x]
                bx = [u for u in others if u != x]
                hbx = hamiltonian_path_count_dp(T, bx)
                sx_hbx += sx * hbx

            # Cycle terms: L>=5 ONLY
            long_D = count_cycles_using_arc_by_length(T, n, i, j, min_length=5)
            long_C = count_cycles_using_arc_by_length(T2, n, j, i, min_length=5)

            formula = -2 * sx_hbx + 2 * (long_D - long_C)

            # Ground truth
            ht = hamiltonian_path_count_dp(T, all_verts)
            ht2 = hamiltonian_path_count_dp(T2, all_verts)
            delta_H = ht - ht2

            if formula == delta_H:
                n_match += 1
            else:
                print(f"  MISS trial={trial}: formula={formula}, DH={delta_H}, "
                      f"sx_hbx={sx_hbx}, longD={long_D}, longC={long_C}")

            if (trial + 1) % 10 == 0:
                print(f"  Progress: {trial+1}/{n_trials}, match={n_match}/{n_total}")

        print(f"n={n}: {n_match}/{n_total}")
        if n_match == n_total:
            print("*** FORMULA VERIFIED AT n=8! ***")
