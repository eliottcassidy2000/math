#!/usr/bin/env python3
"""
Express adj(i,j) - adj'(j,i) as a subset convolution and test if
it matches Delta_I term by term in the subset size stratification.

Key decomposition:
  adj(i,j) = sum_{S subset V0} f_i(S) * g_j(R)
  adj'(j,i) = sum_{S subset V0} f_j(S) * g_i(R)
  where R = V0 \ S, V0 = V \ {i,j}

  f_i(S) = #{paths through S ending adjacent to i}  (+ 1 if S=empty)
  g_j(R) = #{paths through R starting from j-neighbor} (+ 1 if R=empty)

Stratify by |S|=k:
  Delta_H(k) = sum_{|S|=k} [f_i(S)*g_j(R) - f_j(S)*g_i(R)]
  Delta_H = sum_k Delta_H(k)

Test: does Delta_H(k) match a specific part of Delta_I?

Instance: opus-2026-03-05-S3
"""

import sys
sys.path.insert(0, '.')
from tournament_lib import (
    hamiltonian_path_count, find_odd_cycles, independence_poly_at_fast,
    all_tournaments
)
from itertools import permutations, combinations
import random


def compute_fi_gj(T, i, j):
    """Compute f_i(S) and g_j(R) for all subsets S of V0."""
    n = len(T)
    V0 = [v for v in range(n) if v != i and v != j]
    m = len(V0)
    v0_idx = {v: k for k, v in enumerate(V0)}

    fi = {}  # fi[mask] = f_i(S)
    gj = {}  # gj[mask] = g_j(R)
    fj = {}  # fj[mask] = f_j(S)
    gi = {}  # gi[mask] = g_i(R)

    for mask in range(1 << m):
        S = [V0[k] for k in range(m) if (mask >> k) & 1]
        R = [V0[k] for k in range(m) if not ((mask >> k) & 1)]

        if len(S) == 0:
            fi[mask] = 1
            fj[mask] = 1
        else:
            # f_i(S) = #{Ham paths through S ending at u with T[u][i]=1}
            fi_val = 0
            for u in S:
                if T[u][i]:
                    # Count paths through S ending at u
                    fi_val += h_end_sub(T, S, u)
            fi[mask] = fi_val

            # f_j(S) = #{Ham paths through S ending at u with T[u][j]=1}
            fj_val = 0
            for u in S:
                if T[u][j]:
                    fj_val += h_end_sub(T, S, u)
            fj[mask] = fj_val

        rmask = ((1 << m) - 1) ^ mask
        if len(R) == 0:
            gj[rmask] = 1
            gi[rmask] = 1
        else:
            # g_j(R) = #{Ham paths through R starting at w with T[j][w]=1}
            gj_val = 0
            for w in R:
                if T[j][w]:
                    gj_val += h_start_sub(T, R, w)
            gj[rmask] = gj_val

            # g_i(R) = #{Ham paths through R starting at w with T[i][w]=1}
            gi_val = 0
            for w in R:
                if T[i][w]:
                    gi_val += h_start_sub(T, R, w)
            gi[rmask] = gi_val

    return fi, gj, fj, gi, V0


def h_end_sub(T, verts, target):
    """#Ham paths through verts ending at target."""
    if len(verts) == 1:
        return 1 if target in verts else 0
    vlist = list(verts)
    m = len(vlist)
    idx = {v: k for k, v in enumerate(vlist)}
    if target not in idx:
        return 0
    dp = {}
    for v in range(m):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << m):
        for v in range(m):
            if not ((mask >> v) & 1) or (mask, v) not in dp:
                continue
            for u in range(m):
                if (mask >> u) & 1:
                    continue
                if T[vlist[v]][vlist[u]]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << m) - 1
    return dp.get((full, idx[target]), 0)


def h_start_sub(T, verts, source):
    """#Ham paths through verts starting at source."""
    if len(verts) == 1:
        return 1 if source in verts else 0
    vlist = list(verts)
    m = len(vlist)
    idx = {v: k for k, v in enumerate(vlist)}
    if source not in idx:
        return 0
    dp = {(1 << idx[source], idx[source]): 1}
    for mask in range(1, 1 << m):
        for v in range(m):
            if not ((mask >> v) & 1) or (mask, v) not in dp:
                continue
            for u in range(m):
                if (mask >> u) & 1:
                    continue
                if T[vlist[v]][vlist[u]]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << m) - 1
    return sum(dp.get((full, v), 0) for v in range(m))


def analyze_stratification(T, i, j):
    """Stratify the subset convolution by |S| and compare with Delta_I."""
    n = len(T)
    fi, gj, fj, gi, V0 = compute_fi_gj(T, i, j)
    m = len(V0)

    s = {x: 1 - T[x][i] - T[j][x] for x in V0}

    # Stratify by |S|
    delta_by_size = {}
    for k in range(m + 1):
        delta_k = 0
        for mask in range(1 << m):
            if bin(mask).count('1') != k:
                continue
            rmask = ((1 << m) - 1) ^ mask
            delta_k += fi[mask] * gj.get(rmask, 0) - fj[mask] * gi.get(rmask, 0)
        delta_by_size[k] = delta_k

    total = sum(delta_by_size.values())

    # Compute actual adj counts for verification
    adj_ij = sum(1 for perm in permutations(range(n))
                 if all(T[perm[p]][perm[p+1]] for p in range(n-1))
                 for p in range(n-1) if perm[p] == i and perm[p+1] == j)

    Tp = [row[:] for row in T]
    Tp[i][j] = 0
    Tp[j][i] = 1
    adj_ji_p = sum(1 for perm in permutations(range(n))
                   if all(Tp[perm[p]][perm[p+1]] for p in range(n-1))
                   for p in range(n-1) if perm[p] == j and perm[p+1] == i)

    delta_H = adj_ij - adj_ji_p

    # Delta_I
    cycles_T = find_odd_cycles(T)
    cycles_Tp = find_odd_cycles(Tp)
    I_T = independence_poly_at_fast(cycles_T, 2)
    I_Tp = independence_poly_at_fast(cycles_Tp, 2)
    delta_I = I_T - I_Tp

    print(f"\nn={n}, arc ({i},{j})")
    print(f"s-values: {s}")
    print(f"adj(i,j)={adj_ij}, adj'(j,i)={adj_ji_p}")
    print(f"delta_H = {delta_H}, delta_I = {delta_I}, match={delta_H == delta_I}")
    print(f"total from convolution = {total}")
    print(f"\nStratification by |S|:")
    for k in range(m + 1):
        print(f"  |S|={k}: delta = {delta_by_size[k]}")

    # Now also stratify Delta_I by cycle length
    # 3-cycle contribution to Delta_I: -2*sum(s_x) for each x
    d3 = -2 * sum(s[x] for x in V0)
    print(f"\n3-cycle contribution to Delta_I: {d3}")
    print(f"|S|=0 + |S|={m} contributions: {delta_by_size[0] + delta_by_size[m]}")
    print(f"|S|=1 + |S|={m-1} contributions: {delta_by_size.get(1,0) + delta_by_size.get(m-1,0)}")


def main():
    rng = random.Random(42)
    for n in [5, 6, 7]:
        print(f"\n{'='*70}")
        T = [[0]*n for _ in range(n)]
        for a in range(n):
            for b in range(a+1, n):
                if rng.random() < 0.5:
                    T[a][b] = 1
                else:
                    T[b][a] = 1
        ii, jj = 0, 1
        if T[ii][jj] == 0:
            ii, jj = jj, ii
        analyze_stratification(T, ii, jj)


if __name__ == "__main__":
    main()
