"""
Q-009: Verify the algebraic structure of the proof at n=5.

At n=4: D(a) + D(n-2-a) = -(s_a+s_b) (gamma cancels).
At n=5: D(a) contributions should sum to -2*sum(s_x*H(B_x)) + 2*(D5-C5)?

Test: decompose delta into D(a) for a=0,...,3 and check each.
Also check: does D(0)+D(3) and D(1)+D(2) have nice forms?

Author: opus-2026-03-05-S2
"""

import random
from itertools import combinations


def h_end_dp(T, verts, target):
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
            if not (mask & (1<<v)) or dp[mask][v]==0:
                continue
            for u in range(n):
                if mask & (1<<u): continue
                if T[verts[v]][verts[u]]:
                    dp[mask|(1<<u)][u] += dp[mask][v]
    return dp[full][ti]


def h_start_dp(T, verts, source):
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
            if not (mask & (1<<v)) or dp[mask][v]==0:
                continue
            for u in range(n):
                if mask & (1<<u): continue
                if T[verts[v]][verts[u]]:
                    dp[mask|(1<<u)][u] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


def hamiltonian_path_count_dp(T, verts):
    n = len(verts)
    if n <= 1:
        return 1
    dp = [[0]*n for _ in range(1<<n)]
    for k in range(n):
        dp[1<<k][k] = 1
    full = (1<<n)-1
    for mask in range(1, 1<<n):
        for v in range(n):
            if not (mask & (1<<v)) or dp[mask][v]==0:
                continue
            for u in range(n):
                if mask & (1<<u): continue
                if T[verts[v]][verts[u]]:
                    dp[mask|(1<<u)][u] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


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
    n = 5

    for trial in range(10):
        T = random_tournament(n)
        i, j = random.sample(range(n), 2)
        if T[i][j] == 0:
            i, j = j, i

        T2 = flip_arc(T, i, j)
        others = [x for x in range(n) if x != i and x != j]
        all_verts = list(range(n))

        # D(a) for each position
        D = {}
        for a in range(n - 1):  # |S| = a, from 0 to n-2
            val = 0
            for S in combinations(others, a):
                S_set = set(S)
                R = [x for x in others if x not in S_set]
                he_i = h_end_dp(T, sorted(list(S) + [i]), i) if a > 0 else 1
                hs_j = h_start_dp(T, sorted([j] + R), j) if R else 1
                he_j = h_end_dp(T, sorted(list(S) + [j]), j) if a > 0 else 1
                hs_i = h_start_dp(T, sorted([i] + R), i) if R else 1
                val += he_i * hs_j - he_j * hs_i
            D[a] = val

        delta = sum(D.values())

        # RHS
        sx_hbx = sum((1-T[x][i]-T[j][x]) * hamiltonian_path_count_dp(T, [u for u in others if u!=x])
                     for x in others)
        # Cycle term
        dc = 0
        for L in range(5, n+1, 2):
            for extra in combinations(others, L-2):
                verts = list(extra)
                nv = len(verts)
                dp = [[0]*nv for _ in range(1<<nv)]
                for k in range(nv):
                    if T[j][verts[k]]:
                        dp[1<<k][k] = 1
                full = (1<<nv)-1
                for mask in range(1, 1<<nv):
                    for v in range(nv):
                        if not (mask & (1<<v)) or dp[mask][v]==0: continue
                        for u in range(nv):
                            if mask & (1<<u): continue
                            if T[verts[v]][verts[u]]:
                                dp[mask|(1<<u)][u] += dp[mask][v]
                dc += sum(dp[full][v] for v in range(nv) if T[verts[v]][i])
            for extra in combinations(others, L-2):
                verts = list(extra)
                nv = len(verts)
                dp = [[0]*nv for _ in range(1<<nv)]
                for k in range(nv):
                    if T2[i][verts[k]]:
                        dp[1<<k][k] = 1
                full = (1<<nv)-1
                for mask in range(1, 1<<nv):
                    for v in range(nv):
                        if not (mask & (1<<v)) or dp[mask][v]==0: continue
                        for u in range(nv):
                            if mask & (1<<u): continue
                            if T2[verts[v]][verts[u]]:
                                dp[mask|(1<<u)][u] += dp[mask][v]
                dc -= sum(dp[full][v] for v in range(nv) if T2[verts[v]][j])

        rhs = -2 * sx_hbx + 2 * dc

        # Symmetric pairing
        s_vals = {x: 1-T[x][i]-T[j][x] for x in others}
        h_vals = {x: hamiltonian_path_count_dp(T, [u for u in others if u!=x]) for x in others}

        print(f"\nTrial {trial}: flip {i}->{j}, delta={delta}, rhs={rhs}")
        print(f"  D(a): {D}")
        print(f"  D(0)+D({n-2})={D[0]+D[n-2]}, D(1)+D({n-3})={D[1]+D[n-3]}")
        if n >= 6:
            print(f"  D(2)+D({n-4})={D[2]+D[n-4]}")
        print(f"  s_x: {s_vals}")
        print(f"  H(B_x): {h_vals}")
        print(f"  -2*sum(s*H): {-2*sx_hbx}, 2*dc: {2*dc}")
