#!/usr/bin/env python3
"""
Decompose U_T and U_T' by blocking vertex to understand the
per-vertex structure of the polynomial identity U_T' - U_T = delta_I.

Key question: does U_T decompose as sum_x f(x) where f(x) relates to s_x * H(B_x)?

Instance: opus-2026-03-05-S3
"""

import sys
sys.path.insert(0, '.')
from tournament_lib import hamiltonian_path_count, find_odd_cycles, independence_poly_at_fast
from itertools import permutations, combinations
import random


def analyze_decomposition(T, i, j):
    """Decompose U_T by which vertex blocks the swap."""
    n = len(T)
    assert T[i][j] == 1
    others = [v for v in range(n) if v != i and v != j]

    # Compute s-values
    s = {}
    for x in others:
        s[x] = 1 - T[x][i] - T[j][x]

    # Build T'
    Tp = [row[:] for row in T]
    Tp[i][j] = 0
    Tp[j][i] = 1

    # Find all Ham paths of T using arc i->j
    all_paths_T = []
    for perm in permutations(range(n)):
        valid = all(T[perm[k]][perm[k+1]] for k in range(n-1))
        if not valid:
            continue
        # Check if i->j is used consecutively
        for k in range(n-1):
            if perm[k] == i and perm[k+1] == j:
                all_paths_T.append(perm)
                break

    # Find all Ham paths of T' using arc j->i
    all_paths_Tp = []
    for perm in permutations(range(n)):
        valid = all(Tp[perm[k]][perm[k+1]] for k in range(n-1))
        if not valid:
            continue
        for k in range(n-1):
            if perm[k] == j and perm[k+1] == i:
                all_paths_Tp.append(perm)
                break

    # Classify T-paths by blocking type
    # A T-path (..., x, i, j, y, ...) is unmatched iff:
    #   pred-blocked: T[x][j]=0 (i.e. s_x = -1 when T[x][i]=1)
    #   succ-blocked: T[i][y]=0 (i.e. s_y relates to T[j][y])
    # Wait, let me be more careful about s classification.
    # s_x = 1 - T[x][i] - T[j][x]
    # pred-blocked at x: x is predecessor of i, need T[x][j]=0
    #   T[x][j]=0 means j->x, i.e. T[j][x]=1
    #   So pred is x, T[x][i]=1 (since path has x->i), T[j][x]=1
    #   s_x = 1 - 1 - 1 = -1  ✓ (s=-1 blocks T-paths via pred)
    # succ-blocked at y: y is successor of j, need T[i][y]=0
    #   T[i][y]=0 means y->i, i.e. T[y][i]=1
    #   Path has j->y so T[j][y]=1 (but wait, T[j][y] not T'[j][y])
    #   Hmm, path is a T-path so it uses T arcs. j->y means T[j][y]=1.
    #   s_y = 1 - T[y][i] - T[j][y] = 1 - 1 - 1 = -1  ✓

    # So unmatched T-paths are blocked only by s=-1 vertices. Good.

    # Decompose by blocking vertex
    U_T_by_pred = {x: 0 for x in others}  # paths blocked by pred=x
    U_T_by_succ = {x: 0 for x in others}  # paths blocked by succ=x
    U_T_pred_only = {x: 0 for x in others}
    U_T_succ_only = {x: 0 for x in others}
    U_T_both = {}  # blocked by both pred and succ
    uT_total = 0

    for path in all_paths_T:
        pos = list(path).index(i)
        assert path[pos+1] == j
        pred_blocked = False
        succ_blocked = False
        pred_v = None
        succ_v = None
        if pos > 0:
            pred_v = path[pos-1]
            if T[pred_v][j] == 0:  # j beats pred
                pred_blocked = True
        if pos+1 < n-1:
            succ_v = path[pos+2]
            if T[i][succ_v] == 0:  # succ beats i
                succ_blocked = True

        if pred_blocked or succ_blocked:
            uT_total += 1
            if pred_blocked:
                U_T_by_pred[pred_v] += 1
            if succ_blocked:
                U_T_by_succ[succ_v] += 1
            if pred_blocked and not succ_blocked:
                U_T_pred_only[pred_v] += 1
            if succ_blocked and not pred_blocked:
                U_T_succ_only[succ_v] += 1
            if pred_blocked and succ_blocked:
                key = (pred_v, succ_v)
                U_T_both[key] = U_T_both.get(key, 0) + 1

    # Same for T'-paths
    U_Tp_by_pred = {x: 0 for x in others}
    U_Tp_by_succ = {x: 0 for x in others}
    uTp_total = 0

    for path in all_paths_Tp:
        pos = list(path).index(j)
        assert path[pos+1] == i
        pred_blocked = False
        succ_blocked = False
        pred_v = None
        succ_v = None
        if pos > 0:
            pred_v = path[pos-1]
            if Tp[pred_v][i] == 0:
                pred_blocked = True
        if pos+1 < n-1:
            succ_v = path[pos+2]
            if Tp[j][succ_v] == 0:
                succ_blocked = True

        if pred_blocked or succ_blocked:
            uTp_total += 1
            if pred_blocked:
                U_Tp_by_pred[pred_v] += 1
            if succ_blocked:
                U_Tp_by_succ[succ_v] += 1

    # Compute H(B_x) for each x
    H_Bx = {}
    for x in others:
        B_x = [v for v in others if v != x]
        # Build sub-tournament
        m = len(B_x)
        sub = [[0]*m for _ in range(m)]
        for a_idx in range(m):
            for b_idx in range(m):
                sub[a_idx][b_idx] = T[B_x[a_idx]][B_x[b_idx]]
        H_Bx[x] = hamiltonian_path_count(sub)

    # Print analysis
    print(f"\nn={n}, arc ({i},{j})")
    print(f"s-values: {s}")
    print(f"U_T = {uT_total}, U_T' = {uTp_total}, diff = {uTp_total - uT_total}")

    # THM-013 formula
    delta_I_3cycle = -2 * sum(s[x] for x in others)
    print(f"-2*sum(s_x) = {delta_I_3cycle}")
    print(f"2*sum(s_x*H(B_x)) = {2 * sum(s[x]*H_Bx[x] for x in others)}")

    print(f"\nPer-vertex breakdown:")
    print(f"{'x':>4} {'s_x':>4} {'H(Bx)':>6} {'U_T pred':>8} {'U_T succ':>8} "
          f"{'U_Tp pred':>9} {'U_Tp succ':>9} {'net_pred':>8} {'net_succ':>8}")
    for x in others:
        net_pred = U_Tp_by_pred[x] - U_T_by_pred[x]
        net_succ = U_Tp_by_succ[x] - U_T_by_succ[x]
        print(f"{x:>4} {s[x]:>4} {H_Bx[x]:>6} {U_T_by_pred[x]:>8} {U_T_by_succ[x]:>8} "
              f"{U_Tp_by_pred[x]:>9} {U_Tp_by_succ[x]:>9} {net_pred:>8} {net_succ:>8}")

    total_net = uTp_total - uT_total
    # Does the net decompose per-vertex as something * s_x * H(B_x)?
    # The inclusion-exclusion double-counts, so total != sum of per-vertex.
    # Let me compute the inclusion-exclusion corrected version.
    double_counted_T = sum(v for v in U_T_both.values())
    double_counted_Tp_both = {}
    for path in all_paths_Tp:
        pos = list(path).index(j)
        pred_v = path[pos-1] if pos > 0 else None
        succ_v = path[pos+2] if pos+1 < n-1 else None
        pb = pred_v is not None and Tp[pred_v][i] == 0
        sb = succ_v is not None and Tp[j][succ_v] == 0
        if pb and sb:
            key = (pred_v, succ_v)
            double_counted_Tp_both[key] = double_counted_Tp_both.get(key, 0) + 1
    double_counted_Tp = sum(v for v in double_counted_Tp_both.values())

    print(f"\nDouble-counted (both pred+succ blocked): T={double_counted_T}, T'={double_counted_Tp}")
    sum_pred_T = sum(U_T_by_pred.values())
    sum_succ_T = sum(U_T_by_succ.values())
    sum_pred_Tp = sum(U_Tp_by_pred.values())
    sum_succ_Tp = sum(U_Tp_by_succ.values())
    print(f"sum_pred_T={sum_pred_T}, sum_succ_T={sum_succ_T}, double={double_counted_T}, "
          f"total={sum_pred_T + sum_succ_T - double_counted_T} (should be {uT_total})")
    print(f"sum_pred_Tp={sum_pred_Tp}, sum_succ_Tp={sum_succ_Tp}, double={double_counted_Tp}, "
          f"total={sum_pred_Tp + sum_succ_Tp - double_counted_Tp} (should be {uTp_total})")


def main():
    rng = random.Random(42)
    # Test at n=5, 6, 7
    for n in [5, 6, 7]:
        print(f"\n{'='*70}")
        print(f"n = {n}")
        print(f"{'='*70}")

        for trial in range(3):
            # Random tournament
            T = [[0]*n for _ in range(n)]
            for a in range(n):
                for b in range(a+1, n):
                    if rng.random() < 0.5:
                        T[a][b] = 1
                    else:
                        T[b][a] = 1

            # Pick an arc to flip
            i, j = 0, 1
            if T[i][j] == 0:
                i, j = j, i
            analyze_decomposition(T, i, j)


if __name__ == "__main__":
    main()
