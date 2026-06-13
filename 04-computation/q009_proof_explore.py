"""
Q-009: Explore proof of adj(i,j) - adj'(j,i) = DeltaI(Omega, 2).

Strategy: Find a SIGNED INVOLUTION on Hamiltonian paths that cancels
all paths except those "explained" by the cycle structure.

Key observation: a Ham path pi of T using arc i->j can be "paired" with
a modified path using other arcs. The unpaired paths correspond to
cycle-related terms.

Alternative: express adj(i,j) directly in terms of sub-tournament data
and show it matches the RHS.

adj(i,j) = sum_{S subset V\{i,j}} h_end(T[S+{i}], i) * h_start(T[{j}+R], j)
          where R = V\{i,j}\S

This decomposes by "split position" |S| = 0,1,...,n-2.

Author: opus-2026-03-05-S2
"""

import random
from itertools import combinations, permutations


def h_end_dp(T, verts, target):
    """# Ham paths of T[verts] ending at target."""
    idx = {v: k for k, v in enumerate(verts)}
    n = len(verts)
    ti = idx[target]
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
    return dp[full][ti]


def h_start_dp(T, verts, source):
    """# Ham paths of T[verts] starting at source."""
    idx = {v: k for k, v in enumerate(verts)}
    n = len(verts)
    si = idx[source]
    dp = [[0] * n for _ in range(1 << n)]
    dp[1 << si][si] = 1
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
    n = 6

    # For n=6, adj(i,j) - adj'(j,i) = -2*sum(s_x*H(B_x)) + 2*(D5-C5)
    #
    # Key: adj(i,j) = sum_{S} h_end(S+{i}, i) * h_start({j}+R, j)
    #      adj'(j,i) = sum_{S} h_end(S+{j}, j) * h_start({i}+R, i)
    #
    # Difference by split size a = |S|:
    # For each a and each S of size a:
    #   delta_S = h_end(S+{i}, i)*h_start({j}+R, j) - h_end(S+{j}, j)*h_start({i}+R, i)
    #
    # For a=0: S=empty, R=others.
    #   delta = 1*h_start({j}+others, j) - 1*h_start({i}+others, i)
    #   These are "starting Hamiltonian path counts" from j vs i in (n-1)-vertex tournaments
    #
    # For a=1: S={x}, R = others\{x} = B_x
    #   delta_x = T[x][i]*h_start({j}+B_x, j) - T[x][j]*h_start({i}+B_x, i)
    #
    # For a=n-3: S has n-3 vertices = B_x for some x
    #   delta_x = h_end(B_x+{i}, i)*T[j][x] - h_end(B_x+{j}, j)*T[i][x]

    # Let me compute per-x contributions combining all split sizes that involve x
    # on a specific side.

    print(f"Exploring adj decomposition at n={n}")
    print("=" * 70)

    for trial in range(3):
        T = random_tournament(n)
        # Pick an arc
        for i in range(n):
            for j in range(n):
                if i != j and T[i][j]:
                    break
            if T[i][j]:
                break

        others = [x for x in range(n) if x != i and x != j]

        # Compute adj(i,j) and adj'(j,i) by split
        adj_by_split = {}
        adj2_by_split = {}
        for a in range(n - 1):
            adj_a = 0
            adj2_a = 0
            for S in combinations(others, a):
                S_set = set(S)
                R = [x for x in others if x not in S_set]
                he_i = h_end_dp(T, sorted(list(S) + [i]), i)
                hs_j = h_start_dp(T, sorted([j] + R), j)
                he_j = h_end_dp(T, sorted(list(S) + [j]), j)
                hs_i = h_start_dp(T, sorted([i] + R), i)
                adj_a += he_i * hs_j
                adj2_a += he_j * hs_i
            adj_by_split[a] = adj_a
            adj2_by_split[a] = adj2_a

        delta_by_split = {a: adj_by_split[a] - adj2_by_split[a] for a in adj_by_split}
        total_delta = sum(delta_by_split.values())

        # The RHS
        sx_hbx = sum((1 - T[x][i] - T[j][x]) * hamiltonian_path_count_dp(T, [u for u in others if u != x])
                     for x in others)

        # D5-C5
        T2 = flip_arc(T, i, j)
        d5c5 = 0
        for x in others:
            bx = [u for u in others if u != x]
            for perm in permutations(bx):
                if T[perm[0]][perm[1]] * T[perm[1]][perm[2]]:
                    d5c5 += T[j][perm[0]] * T[perm[2]][i] - T[i][perm[0]] * T[perm[2]][j]

        rhs = -2 * sx_hbx + 2 * d5c5

        print(f"\nTrial {trial}: flip {i}->{j}")
        print(f"  Delta by split: {delta_by_split}")
        print(f"  Total: {total_delta}, RHS: {rhs}, match: {total_delta == rhs}")

        # Now check: is there a PAIRING between split a and split (n-2-a)?
        # a=0 pairs with a=n-2 (=4 at n=6)
        # a=1 pairs with a=3
        # a=2 is self-paired
        for a in range(n - 1):
            b = n - 2 - a
            if a <= b:
                pair_sum = delta_by_split[a] + (delta_by_split[b] if b != a else 0)
                print(f"  Pair ({a},{b}): delta[{a}]={delta_by_split[a]}, "
                      f"delta[{b}]={delta_by_split.get(b,0)}, sum={pair_sum}")

        # Check the a=1 + a=3 pair per x
        print(f"  --- Per x analysis ---")
        for x in others:
            bx = [u for u in others if u != x]
            sx = 1 - T[x][i] - T[j][x]
            hbx = hamiltonian_path_count_dp(T, bx)

            # a=1 with S={x}: T[x][i]*h_start({j}+B_x,j) - T[x][j]*h_start({i}+B_x,i)
            a1_x = T[x][i] * h_start_dp(T, sorted([j] + bx), j) - \
                   T[x][j] * h_start_dp(T, sorted([i] + bx), i)

            # a=3 with S=B_x (excluded x): h_end(B_x+{i},i)*T[j][x] - h_end(B_x+{j},j)*T[i][x]
            a3_x = h_end_dp(T, sorted(bx + [i]), i) * T[j][x] - \
                   h_end_dp(T, sorted(bx + [j]), j) * T[i][x]

            # What about a=2? There are C(3,2)=3 ways to pick S of size 2 from others\{x}
            # and each involves x in R. Not a clean per-x decomposition.

            print(f"    x={x}: sx={sx}, H(Bx)={hbx}, a1={a1_x}, a3={a3_x}, "
                  f"a1+a3={a1_x+a3_x}, -2*sx*H={-2*sx*hbx}")

        # Key question: does a1_x + a3_x = -2*s_x*H(B_x)?
        # Or is there a term from a=0 and a=4 that contributes?
        # And what is a=2?
