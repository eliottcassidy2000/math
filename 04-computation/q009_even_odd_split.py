"""
Q-009: Test the "even-odd split" conjecture.

Conjecture: In the adj decomposition delta = sum_S Delta(S, others\S),
the sum over even-sized S equals the sum over odd-sized S.
If true, this would mean delta = 2 * (odd-S sum), connecting directly
to the cycle formula (which only involves odd cycle lengths).

This is a KEY structural observation that could lead to a general proof.

Author: opus-2026-03-05-S4
"""

import random
from itertools import permutations, combinations


def random_tournament(n):
    T = [[0]*n for _ in range(n)]
    for a in range(n):
        for b in range(a+1, n):
            if random.random() < 0.5:
                T[a][b] = 1
            else:
                T[b][a] = 1
    return T


def h_dp(T, verts):
    """Count Hamiltonian paths."""
    n = len(verts)
    if n <= 1:
        return 1
    dp = [[0]*n for _ in range(1<<n)]
    for k in range(n):
        dp[1<<k][k] = 1
    full = (1<<n)-1
    for mask in range(1, 1<<n):
        for v in range(n):
            if not (mask & (1<<v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1<<u):
                    continue
                if T[verts[v]][verts[u]]:
                    dp[mask|(1<<u)][u] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


def h_end_dp(T, verts, target):
    """Count Ham paths on verts ending at target."""
    n = len(verts)
    if n == 1:
        return 1
    idx = {v: k for k, v in enumerate(verts)}
    ti = idx[target]
    dp = [[0]*n for _ in range(1<<n)]
    for k in range(n):
        dp[1<<k][k] = 1
    full = (1<<n)-1
    for mask in range(1, 1<<n):
        for v in range(n):
            if not (mask & (1<<v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1<<u):
                    continue
                if T[verts[v]][verts[u]]:
                    dp[mask|(1<<u)][u] += dp[mask][v]
    return dp[full][ti]


def h_start_dp(T, verts, source):
    """Count Ham paths on verts starting at source."""
    n = len(verts)
    if n == 1:
        return 1
    idx = {v: k for k, v in enumerate(verts)}
    si = idx[source]
    dp = [[0]*n for _ in range(1<<n)]
    dp[1<<si][si] = 1
    full = (1<<n)-1
    for mask in range(1, 1<<n):
        for v in range(n):
            if not (mask & (1<<v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1<<u):
                    continue
                if T[verts[v]][verts[u]]:
                    dp[mask|(1<<u)][u] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


def compute_Delta(T, i, j, S, R):
    """Compute Delta(S, R) = L_j(S)*R_i(R) - L_i(S)*R_j(R).

    L_j(S) = sum_{a in S} h_end(S, a) * T[a][j]
    L_i(S) = sum_{a in S} h_end(S, a) * T[a][i]
    R_i(R) = sum_{b in R} h_start(R, b) * T[i][b]
    R_j(R) = sum_{b in R} h_start(R, b) * T[j][b]
    """
    S_list = sorted(S)
    R_list = sorted(R)

    if not S_list:
        L_j = 1  # empty prefix: direct to j
        L_i = 1  # empty prefix: direct to i
    else:
        L_j = sum(h_end_dp(T, S_list, a) * T[a][j] for a in S_list)
        L_i = sum(h_end_dp(T, S_list, a) * T[a][i] for a in S_list)

    if not R_list:
        R_i = 1
        R_j = 1
    else:
        R_i = sum(h_start_dp(T, R_list, b) * T[i][b] for b in R_list)
        R_j = sum(h_start_dp(T, R_list, b) * T[j][b] for b in R_list)

    return L_j * R_i - L_i * R_j


def compute_cycle_contribution(T, i, j, S, R):
    """Compute [g(S) - l(S)] * H(R) for a subset S.

    g(S) = #{Ham paths on S: i -> sigma -> j}
    l(S) = #{Ham paths on S: j -> sigma -> i}
    """
    S_list = sorted(S)
    R_list = sorted(R)

    H_R = h_dp(T, R_list)

    if not S_list:
        # 2-cycle (i,j): gained iff T[i][j] is the flipped arc... skip for now
        return 0

    g = 0
    l = 0
    for perm in permutations(S_list):
        # Gained: i -> perm -> j
        valid = T[i][perm[0]]
        for k in range(len(perm)-1):
            valid *= T[perm[k]][perm[k+1]]
        valid *= T[perm[-1]][j]
        g += valid

        # Lost: j -> perm -> i
        valid = T[j][perm[0]]
        for k in range(len(perm)-1):
            valid *= T[perm[k]][perm[k+1]]
        valid *= T[perm[-1]][i]
        l += valid

    return (g - l) * H_R


if __name__ == "__main__":
    random.seed(42)

    for n in [5, 6, 7, 8]:
        print(f"\n{'='*60}")
        print(f"n = {n}")

        for trial in range(20 if n <= 7 else 10):
            T = random_tournament(n)
            i, j = random.sample(range(n), 2)
            if T[i][j] == 0:
                i, j = j, i

            others = [x for x in range(n) if x != i and x != j]
            m = len(others)

            even_sum = 0
            odd_sum = 0
            cycle_sum = 0

            for size in range(m + 1):
                for S in combinations(others, size):
                    R = tuple(x for x in others if x not in S)
                    d = compute_Delta(T, i, j, S, R)

                    if size % 2 == 0:
                        even_sum += d
                    else:
                        odd_sum += d

                    # Cycle contribution (only for odd |S| >= 1,
                    # corresponding to odd cycle length |S|+2 >= 3)
                    if size >= 1 and size % 2 == 1:
                        cycle_sum += compute_cycle_contribution(T, i, j, S, R)

            delta_total = even_sum + odd_sum

            # Also compute delta_H directly
            T2 = [row[:] for row in T]
            T2[i][j] = 0
            T2[j][i] = 1
            H_T = h_dp(T, list(range(n)))
            H_T2 = h_dp(T2, list(range(n)))
            delta_H = H_T2 - H_T

            if trial < 5 or even_sum != odd_sum:
                print(f"  Trial {trial}: delta_H={delta_H}, delta_decomp={delta_total}, "
                      f"even={even_sum}, odd={odd_sum}, "
                      f"even==odd: {even_sum==odd_sum}, "
                      f"2*cycle_sum={2*cycle_sum}, match: {delta_H==2*cycle_sum}")

        # Summary
        print(f"  Checking even==odd for n={n}...")
        failures = 0
        for trial in range(100 if n <= 7 else 30):
            T = random_tournament(n)
            i, j = random.sample(range(n), 2)
            if T[i][j] == 0:
                i, j = j, i
            others = [x for x in range(n) if x != i and x != j]
            m = len(others)

            even_s = 0
            odd_s = 0
            for size in range(m + 1):
                for S in combinations(others, size):
                    R = tuple(x for x in others if x not in S)
                    d = compute_Delta(T, i, j, S, R)
                    if size % 2 == 0:
                        even_s += d
                    else:
                        odd_s += d
            if even_s != odd_s:
                failures += 1
        total = 100 if n <= 7 else 30
        print(f"  Even==Odd: {total-failures}/{total}")
