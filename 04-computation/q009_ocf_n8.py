"""
Q-009: Verify OCF (H(T) = I(Omega(T), 2)) at n=8.

Need to verify the foundation before fixing the delta formula.

Author: opus-2026-03-05-S2
"""

import random
from itertools import combinations, permutations


def hamiltonian_path_count_dp(T, n):
    verts = list(range(n))
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
                if T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


def find_all_directed_odd_cycles(T, n):
    """Find all directed odd cycles. Each cycle represented as tuple starting from min vertex."""
    cycles = []
    for L in range(3, n + 1, 2):
        for verts in combinations(range(n), L):
            v0 = verts[0]
            rest = list(verts[1:])
            for perm in permutations(rest):
                cycle = (v0,) + perm
                # Check directed cycle
                valid = True
                for k in range(L):
                    if T[cycle[k]][cycle[(k+1) % L]] != 1:
                        valid = False
                        break
                if valid:
                    cycles.append(cycle)
    return cycles


def independence_poly_at_2(cycles):
    """Compute I(Omega, 2) where Omega is the conflict graph on cycles."""
    nc = len(cycles)
    vertex_sets = [frozenset(c) for c in cycles]

    # alpha_0 = 1
    # alpha_1 = nc
    # alpha_k = # independent sets of size k
    # I(Omega, 2) = sum_{k>=0} alpha_k * 2^k

    # For efficiency, compute alpha_k up to whatever size is needed.
    # At n=8, max independent set size is 2 (3+5=8 or 3+3=6<=8).
    # Actually, could have 3+3+... but 3+3+3=9>8, so max is 2.

    alpha_0 = 1
    alpha_1 = nc
    alpha_2 = 0
    for a in range(nc):
        for b in range(a + 1, nc):
            if not (vertex_sets[a] & vertex_sets[b]):
                alpha_2 += 1

    # At n=8, can we have independent sets of size 3?
    # Need 3 mutually VD odd cycles. Min total vertices = 3+3+3=9 > 8. No.
    # So alpha_k = 0 for k >= 3.

    return alpha_0 + 2 * alpha_1 + 4 * alpha_2


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
    n_match = 0
    n_total = 0

    for trial in range(20):
        T = random_tournament(n)
        n_total += 1

        H = hamiltonian_path_count_dp(T, n)
        cycles = find_all_directed_odd_cycles(T, n)
        I2 = independence_poly_at_2(cycles)

        if H == I2:
            n_match += 1
        else:
            nc = len(cycles)
            # Count by length
            by_len = {}
            for c in cycles:
                L = len(c)
                by_len[L] = by_len.get(L, 0) + 1
            print(f"  Trial {trial}: H={H}, I(2)={I2}, diff={H-I2}, "
                  f"#cycles={nc}, by_len={by_len}")

        if (trial + 1) % 10 == 0:
            print(f"Progress: {n_match}/{n_total}")

    print(f"\nn={n}: OCF holds for {n_match}/{n_total} random tournaments")
    if n_match == n_total:
        print("*** OCF VERIFIED AT n=8 ***")
