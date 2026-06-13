"""
Q-009: Test general formula at n=9 where alpha_3 first becomes nonzero.

Three VD 3-cycles need 9 vertices, so n=9 is the first size where
alpha_3(Omega) can be positive, adding an 8*Delta(alpha_3) term.

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
    """Count directed odd cycles in T[verts]."""
    n = len(verts)
    count = 0
    for L in range(3, n + 1, 2):
        for subset in combinations(verts, L):
            v0 = subset[0]
            for perm in permutations(subset[1:]):
                cycle = (v0,) + perm
                valid = True
                for k in range(L):
                    if T[cycle[k]][cycle[(k+1) % L]] != 1:
                        valid = False
                        break
                if valid:
                    count += 1
    return count


def count_vd_pairs_in_sub(T, verts):
    """Count VD odd-cycle pairs in T[verts]."""
    cycles = []
    n = len(verts)
    for L in range(3, n + 1, 2):
        for subset in combinations(verts, L):
            v0 = subset[0]
            for perm in permutations(subset[1:]):
                cycle = (v0,) + perm
                valid = True
                for k in range(L):
                    if T[cycle[k]][cycle[(k+1) % L]] != 1:
                        valid = False
                        break
                if valid:
                    cycles.append(frozenset(subset))
    count = 0
    for a in range(len(cycles)):
        for b in range(a + 1, len(cycles)):
            if not (cycles[a] & cycles[b]):
                count += 1
    return count


def count_vd_triples_in_sub(T, verts):
    """Count VD odd-cycle triples in T[verts]."""
    cycles = []
    n = len(verts)
    for L in range(3, n + 1, 2):
        for subset in combinations(verts, L):
            v0 = subset[0]
            for perm in permutations(subset[1:]):
                cycle = (v0,) + perm
                valid = True
                for k in range(L):
                    if T[cycle[k]][cycle[(k+1) % L]] != 1:
                        valid = False
                        break
                if valid:
                    cycles.append(frozenset(subset))
    count = 0
    nc = len(cycles)
    for a in range(nc):
        for b in range(a + 1, nc):
            if cycles[a] & cycles[b]:
                continue
            for c in range(b + 1, nc):
                if not (cycles[a] & cycles[c]) and not (cycles[b] & cycles[c]):
                    count += 1
    return count


def cycles_using_arc_info(T, n, i, j):
    """For L-cycles using arc i->j, compute count, alpha_1(comp), alpha_2(comp)."""
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
            n_cyc = sum(dp[full][v] for v in range(nv) if T[verts[v]][i])
            if n_cyc > 0:
                comp = sorted(all_verts - set(extra) - {i, j})
                a1 = count_odd_cycles_in_sub(T, comp)
                a2 = count_vd_pairs_in_sub(T, comp)
                count += n_cyc
                a1_sum += n_cyc * a1
                a2_sum += n_cyc * a2
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
    n = 9
    n_trials = 5
    n_match = 0
    all_verts = list(range(n))

    for trial in range(n_trials):
        print(f"Trial {trial}...", flush=True)
        T = random_tournament(n)
        i, j = random.sample(range(n), 2)
        if T[i][j] == 0:
            i, j = j, i

        T2 = flip_arc(T, i, j)

        ht = hamiltonian_path_count_dp(T, all_verts)
        ht2 = hamiltonian_path_count_dp(T2, all_verts)
        delta_H = ht - ht2

        info_T = cycles_using_arc_info(T, n, i, j)
        info_T2 = cycles_using_arc_info(T2, n, j, i)

        delta_a1 = sum(info_T.get(L, (0,0,0))[0] - info_T2.get(L, (0,0,0))[0]
                      for L in range(3, n+1, 2))
        delta_a2 = sum(info_T.get(L, (0,0,0))[1] - info_T2.get(L, (0,0,0))[1]
                      for L in range(3, n+1, 2))
        delta_a3 = sum(info_T.get(L, (0,0,0))[2] - info_T2.get(L, (0,0,0))[2]
                      for L in range(3, n+1, 2))

        formula = 2 * delta_a1 + 4 * delta_a2 + 8 * delta_a3

        if formula == delta_H:
            n_match += 1
            print(f"  MATCH: DH={delta_H}, da1={delta_a1}, da2={delta_a2}, da3={delta_a3}")
        else:
            print(f"  MISS: DH={delta_H}, formula={formula}, "
                  f"da1={delta_a1}, da2={delta_a2}, da3={delta_a3}")

    print(f"\nn={n}: {n_match}/{n_trials}")
    if n_match == n_trials:
        print("*** FORMULA VERIFIED AT n=9 ***")
