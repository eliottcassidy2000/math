"""
Q-009: Connect nadj terms to cycle formula.

Key identity to prove:
  adj(i,j) - adj'(j,i) = #U_T - #U_T'
  where #U_T = sum_{x:s=-1} [nadj(x,i,j) + nadj(i,j,x)]
             - sum_{x,y:s=-1} nadj(x,i,j,y)
  and #U_T' similarly for s=+1.

Questions:
1. What is nadj(x,i,j) + nadj(i,j,x) in terms of H(B_x)?
2. What is nadj(x,i,j,y) in terms of cycles?
3. Can we match the RHS = -2*sum(s_x*H(B_x)) + 2*sum(DL-CL) term by term?

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


def nadj_count(T, n, seq):
    """Count Ham paths of T containing 'seq' as consecutive subsequence."""
    # seq = (a, b, c, ...) must appear consecutively in the path
    seq_set = set(seq)
    others = [v for v in range(n) if v not in seq_set]
    k = len(seq)
    # Check that seq arcs are valid
    for p in range(k - 1):
        if T[seq[p]][seq[p+1]] != 1:
            return 0
    # The path has form: (L, seq[0], seq[1], ..., seq[-1], R)
    # L uses some subset S of others, R uses the rest
    bx = others  # B = V \ seq
    total = 0
    for S_size in range(len(bx) + 1):
        for S in combinations(bx, S_size):
            S_set = set(S)
            R = [v for v in bx if v not in S_set]
            # h_end(S + {seq[0]}, seq[0]) * h_start({seq[-1]} + R, seq[-1])
            # But need to handle first/last vertex connections
            he = h_end_dp(T, sorted(list(S) + [seq[0]]), seq[0]) if S_size > 0 or True else 1
            hs = h_start_dp(T, sorted([seq[-1]] + R), seq[-1]) if R or True else 1
            total += he * hs
    return total


def h_end_dp(T, verts, target):
    n = len(verts)
    if n == 1:
        return 1
    idx = {v: k for k, v in enumerate(verts)}
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
    n = len(verts)
    if n == 1:
        return 1
    idx = {v: k for k, v in enumerate(verts)}
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

    for n in [5, 6, 7]:
        print(f"\n{'='*60}")
        print(f"n = {n}")

        for trial in range(5):
            T = random_tournament(n)
            i, j = random.sample(range(n), 2)
            if T[i][j] == 0:
                i, j = j, i

            T2 = flip_arc(T, i, j)
            others = [x for x in range(n) if x != i and x != j]

            # Compute per-vertex nadj terms
            for x in others:
                sx = 1 - T[x][i] - T[j][x]
                if sx == 0:
                    continue

                bx = [u for u in others if u != x]
                hbx = hamiltonian_path_count_dp(T, bx)

                if sx == -1:
                    # T-side: nadj(x,i,j) and nadj(i,j,x)
                    n1 = nadj_count(T, n, (x, i, j))
                    n2 = nadj_count(T, n, (i, j, x))
                    total_nadj = n1 + n2

                    # Also compute the symmetric T' version (for "swapped")
                    # Since sx=-1: T[x][j]=0, so nadj_{T'}(x,j,i) = 0 (x can't precede j)
                    # and T[j][x]=1, T[i][x]? : for succ, need i→x in T', i.e., T[i][x].
                    # T[x][i]=1, so T[i][x]=0. nadj_{T'}(j,i,x) needs T[i][x]=1. No.
                    # So T'-side contribution from x is 0. ✓

                elif sx == 1:
                    # T'-side: nadj_{T'}(x,j,i) and nadj_{T'}(j,i,x)
                    n1 = nadj_count(T2, n, (x, j, i))
                    n2 = nadj_count(T2, n, (j, i, x))
                    total_nadj = n1 + n2

                if trial < 2:
                    print(f"  Trial {trial}, x={x}: s={sx}, H(B_x)={hbx}, "
                          f"nadj_pred={n1}, nadj_succ={n2}, total={total_nadj}")

            # Now compute the full formula comparison
            # LHS = adj(i,j) - adj'(j,i)
            adj_ij = nadj_count(T, n, (i, j))
            adj_ji = nadj_count(T2, n, (j, i))
            delta = adj_ij - adj_ji

            # #U_T and #U_T' from nadj terms
            U_T_pred = sum(nadj_count(T, n, (x, i, j)) for x in others if 1-T[x][i]-T[j][x] == -1)
            U_T_succ = sum(nadj_count(T, n, (i, j, x)) for x in others if 1-T[x][i]-T[j][x] == -1)
            U_T_both = 0
            s_neg = [x for x in others if 1-T[x][i]-T[j][x] == -1]
            for a in range(len(s_neg)):
                for b in range(len(s_neg)):
                    if a != b:
                        U_T_both += nadj_count(T, n, (s_neg[a], i, j, s_neg[b]))
            U_T = U_T_pred + U_T_succ - U_T_both

            U_T2_pred = sum(nadj_count(T2, n, (x, j, i)) for x in others if 1-T[x][i]-T[j][x] == 1)
            U_T2_succ = sum(nadj_count(T2, n, (j, i, x)) for x in others if 1-T[x][i]-T[j][x] == 1)
            U_T2_both = 0
            s_pos = [x for x in others if 1-T[x][i]-T[j][x] == 1]
            for a in range(len(s_pos)):
                for b in range(len(s_pos)):
                    if a != b:
                        U_T2_both += nadj_count(T2, n, (s_pos[a], j, i, s_pos[b]))
            U_T2 = U_T2_pred + U_T2_succ - U_T2_both

            # RHS
            sx_hbx = sum((1-T[x][i]-T[j][x]) * hamiltonian_path_count_dp(T, [u for u in others if u!=x])
                        for x in others)
            # D5-C5 and longer cycles
            d5c5 = 0
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
                            if not (mask & (1<<v)) or dp[mask][v]==0:
                                continue
                            for u in range(nv):
                                if mask & (1<<u):
                                    continue
                                if T[verts[v]][verts[u]]:
                                    dp[mask|(1<<u)][u] += dp[mask][v]
                    d5c5 += sum(dp[full][v] for v in range(nv) if T[verts[v]][i])
                    # Same for T2, j->i
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
                            if not (mask & (1<<v)) or dp[mask][v]==0:
                                continue
                            for u in range(nv):
                                if mask & (1<<u):
                                    continue
                                if T2[verts[v]][verts[u]]:
                                    dp[mask|(1<<u)][u] += dp[mask][v]
                    d5c5 -= sum(dp[full][v] for v in range(nv) if T2[verts[v]][j])

            rhs = -2 * sx_hbx + 2 * d5c5

            if trial < 2 or delta != U_T - U_T2:
                print(f"  Trial {trial}: delta={delta}, U_T-U_T2={U_T-U_T2}, "
                      f"U_T={U_T} (p={U_T_pred},s={U_T_succ},b={U_T_both}), "
                      f"U_T2={U_T2} (p={U_T2_pred},s={U_T2_succ},b={U_T2_both})")
                print(f"    RHS={rhs}, match: {delta==rhs}")

                # Per-vertex analysis: is nadj_pred+nadj_succ related to H(B_x)?
                for x in s_neg:
                    n_pred = nadj_count(T, n, (x, i, j))
                    n_succ = nadj_count(T, n, (i, j, x))
                    hbx = hamiltonian_path_count_dp(T, [u for u in others if u!=x])
                    print(f"    x={x} (s=-1): pred={n_pred}, succ={n_succ}, "
                          f"sum={n_pred+n_succ}, H(Bx)={hbx}")

                for x in s_pos:
                    n_pred = nadj_count(T2, n, (x, j, i))
                    n_succ = nadj_count(T2, n, (j, i, x))
                    hbx = hamiltonian_path_count_dp(T, [u for u in others if u!=x])
                    print(f"    x={x} (s=+1): pred'={n_pred}, succ'={n_succ}, "
                          f"sum'={n_pred+n_succ}, H(Bx)={hbx}")
