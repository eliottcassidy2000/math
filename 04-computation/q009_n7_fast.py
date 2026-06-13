"""
Q-009: Test ΔI formula at n=7 using DP for Hamiltonian path counting.

Uses bitmask DP: dp[mask][v] = # Ham paths through vertices in mask ending at v.

Author: opus-2026-03-05-S2
"""

import sys
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import find_odd_cycles, conflict_graph, independence_poly_at
from itertools import combinations
import random


def hamiltonian_path_count_dp(T, verts=None):
    """Count Ham paths using bitmask DP. verts = subset to use (default all)."""
    if verts is None:
        verts = list(range(len(T)))
    n = len(verts)
    if n <= 1:
        return 1
    # Map to 0..n-1
    idx = {v: i for i, v in enumerate(verts)}
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
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


def adj_count_dp(T, verts, i, j):
    """Count Ham paths of T[verts] where i immediately precedes j."""
    n = len(verts)
    idx = {v: k for k, v in enumerate(verts)}
    if i not in idx or j not in idx:
        return 0
    ii, jj = idx[i], idx[j]
    dp = [[0] * n for _ in range(1 << n)]
    for k in range(n):
        dp[1 << k][k] = 1
    full = (1 << n) - 1
    count = 0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[verts[v]][verts[u]]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
                    # Check if this is i->j and mask covers all
                    if v == ii and u == jj and (mask | (1 << u)) == full:
                        count += dp[mask][v]
    return count


def count_cycles_using_arc_dp(T, n, i, j):
    """Count directed odd cycles using arc i→j, by length, using DP.
    A cycle using i→j: i→j→...→i (path from j back to i through subset)."""
    counts = {}
    others = [u for u in range(n) if u != i and u != j]
    for L in range(3, n + 1, 2):
        extra_needed = L - 2
        if extra_needed > len(others):
            break
        count = 0
        for extra in combinations(others, extra_needed):
            # Count paths from j to i through 'extra' vertices
            verts = list(extra)
            nv = len(verts)
            if nv == 0:
                # L=2, impossible for odd
                continue
            # DP: paths from j through verts to i
            # dp[mask][v] = # paths starting at j, going through vertices
            # indicated by mask (subset of verts), ending at v
            idx = {v: k for k, v in enumerate(verts)}
            dp = [0] * (1 << nv)
            # Initialize: j -> first vertex
            for k in range(nv):
                if T[j][verts[k]]:
                    dp[1 << k] = 1 if nv == 1 else 0
                    if nv > 1:
                        # Need to track which vertex we're at
                        pass
            # Actually need dp[mask][last_vertex]
            dp2 = [[0] * nv for _ in range(1 << nv)]
            for k in range(nv):
                if T[j][verts[k]]:
                    dp2[1 << k][k] = 1
            full = (1 << nv) - 1
            for mask in range(1, 1 << nv):
                for v in range(nv):
                    if not (mask & (1 << v)) or dp2[mask][v] == 0:
                        continue
                    for u in range(nv):
                        if mask & (1 << u):
                            continue
                        if T[verts[v]][verts[u]]:
                            dp2[mask | (1 << u)][u] += dp2[mask][v]
            # Count paths ending at vertex that beats i
            for v in range(nv):
                if dp2[full][v] > 0 and T[verts[v]][i]:
                    count += dp2[full][v]
        counts[L] = count
    return counts


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def random_tournament(n):
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1
    return T


def count_vd_pairs(cycles):
    count = 0
    for a in range(len(cycles)):
        for b in range(a + 1, len(cycles)):
            if len(cycles[a]) == 3 and len(cycles[b]) == 3:
                if not (set(cycles[a]) & set(cycles[b])):
                    count += 1
    return count


if __name__ == "__main__":
    n = 7
    print(f"Q-009: Test ΔI formula at n={n} (DP-accelerated)")
    print("=" * 70)

    # Check: I(Ω,2) = 1 + 2·|cycles| + 4·|VD pairs|
    print("Checking I(Ω,2) formula at n=7...")
    random.seed(42)
    io_ok = 0
    io_total = 0
    for _ in range(30):
        T = random_tournament(n)
        io_total += 1
        cycles = find_odd_cycles(T)
        if cycles:
            cg = conflict_graph(cycles)
            io_gt = independence_poly_at(cg, 2)
        else:
            io_gt = 1
        vd = count_vd_pairs(cycles)
        io_formula = 1 + 2 * len(cycles) + 4 * vd
        if io_gt == io_formula:
            io_ok += 1
        else:
            print(f"  MISMATCH: gt={io_gt}, formula={io_formula}, diff={io_gt-io_formula}")
            # Debug: find independent sets of size > 2
            for a in range(len(cycles)):
                for b in range(a+1, len(cycles)):
                    sa, sb = set(cycles[a]), set(cycles[b])
                    if not (sa & sb):
                        la, lb = len(cycles[a]), len(cycles[b])
                        if la != 3 or lb != 3:
                            print(f"    Non-3-3 VD pair: len {la} and {lb}")
            break

    print(f"I(Ω,2) = 1+2|cycles|+4|VD|: {io_ok}/{io_total}")

    # Test ΔI formula
    print(f"\nTesting ΔI formula...")
    random.seed(123)
    total = 0
    match = 0

    for trial in range(100):
        T = random_tournament(n)
        i, j = random.sample(range(n), 2)
        if T[i][j] == 0:
            i, j = j, i
        assert T[i][j] == 1
        total += 1

        T2 = flip_arc(T, i, j)
        others = [x for x in range(n) if x != i and x != j]

        # -2·Σ s_x·H(B_x)
        sx_hbx_sum = 0
        for x in others:
            sx = 1 - T[x][i] - T[j][x]
            bx = [u for u in others if u != x]
            hbx = hamiltonian_path_count_dp(T, bx)
            sx_hbx_sum += sx * hbx

        # Cycle terms
        dc_T = count_cycles_using_arc_dp(T, n, i, j)
        dc_T2 = count_cycles_using_arc_dp(T2, n, j, i)
        cycle_term = 0
        for L in range(3, n + 1, 2):
            cycle_term += dc_T.get(L, 0) - dc_T2.get(L, 0)

        formula_dI = -2 * sx_hbx_sum + 2 * cycle_term

        # Ground truth
        ht = hamiltonian_path_count_dp(T)
        ht2 = hamiltonian_path_count_dp(T2)
        delta_H = ht - ht2

        if formula_dI == delta_H:
            match += 1
        else:
            if total - match <= 3:
                print(f"  MISS #{total-match}: trial={trial}, formula={formula_dI}, ΔH={delta_H}")
                print(f"    sx_hbx={sx_hbx_sum}, cycle_term={cycle_term}")
                print(f"    dc_T={dc_T}, dc_T2={dc_T2}")

        if total % 20 == 0:
            print(f"  Progress: {total} done, {match} match")

    print(f"\nΔH = formula: {match}/{total}")

    if match == total:
        print("\n*** FORMULA VERIFIED AT n=7! ***")
        print("ΔI = -2·Σ_x s_x·H(B_x) + 2·Σ_L(DL-CL) generalizes to n=7!")
    else:
        print(f"\nFormula needs refinement at n=7 ({total - match} failures)")
