"""
Q-009: Test the ΔI formula at n=7.

At n=7: max independent set in Ω(T) is still 2 (pairs of VD 3-cycles,
since 3+5=8>7 and 5+5=10>7). So:
  I(Ω,2) = 1 + 2·|cycles| + 4·|{VD 3-cycle pairs}|

The formula should be:
  ΔI = -2·Σ_x s_x · H(B_x) + 2·(D5-C5) + 2·(D7-C7)

where B_x = V\{i,j,x} now has 4 vertices, and D7/C7 count destroyed/created
7-cycles (Hamilton cycles of T using arc i→j / j→i).

At n=7, 7-cycles use all vertices, so D7 = #{Ham cycles of T using i→j} and
C7 = #{Ham cycles of T' using j→i}.

Author: opus-2026-03-05-S2
"""

import sys
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (
    all_tournaments, hamiltonian_path_count, find_odd_cycles,
    conflict_graph, independence_poly_at
)
from itertools import permutations
import random


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def sub_H(T, verts):
    verts = sorted(verts)
    m = len(verts)
    if m <= 1:
        return 1
    count = 0
    for perm in permutations(verts):
        valid = True
        for k in range(m - 1):
            if T[perm[k]][perm[k + 1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count


def count_vd_pairs(cycles):
    count = 0
    for a in range(len(cycles)):
        for b in range(a + 1, len(cycles)):
            if len(cycles[a]) == 3 and len(cycles[b]) == 3:
                if not (set(cycles[a]) & set(cycles[b])):
                    count += 1
    return count


def count_odd_cycles_using_arc(T, n, i, j):
    """Count directed odd cycles using arc i→j, by cycle length."""
    counts = {}  # length -> count
    # For each possible cycle size L (odd, 3 <= L <= n):
    for L in range(3, n + 1, 2):
        if L > n:
            break
        # Choose L-2 more vertices besides i,j
        from itertools import combinations
        others = [u for u in range(n) if u != i and u != j]
        count = 0
        for extra in combinations(others, L - 2):
            # The cycle on {i, j} ∪ extra uses arc i→j
            # Cycle: i→j→p_1→...→p_{L-2}→i
            remaining = list(extra)
            for perm in permutations(remaining):
                valid = True
                # j→perm[0]
                if not T[j][perm[0]]:
                    continue
                for k in range(len(perm) - 1):
                    if not T[perm[k]][perm[k + 1]]:
                        valid = False
                        break
                if valid and T[perm[-1]][i]:
                    count += 1
        counts[L] = count
    return counts


def random_tournament(n):
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1
    return T


if __name__ == "__main__":
    n = 7
    print(f"Q-009: Test ΔI formula at n={n}")
    print("=" * 70)

    # First check: is I(Ω,2) = 1 + 2·|cycles| + 4·|VD pairs| at n=7?
    print("Checking I(Ω,2) formula at n=7...")
    io_check = 0
    io_total = 0
    random.seed(42)
    for _ in range(50):
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
            io_check += 1
        else:
            print(f"  MISMATCH: gt={io_gt}, formula={io_formula}, "
                  f"#cycles={len(cycles)}, #vd={vd}")
            # Check for larger independent sets
            for a in range(len(cycles)):
                for b in range(a+1, len(cycles)):
                    if not (set(cycles[a]) & set(cycles[b])):
                        # Check if there's a third VD cycle
                        ab_verts = set(cycles[a]) | set(cycles[b])
                        for c in range(b+1, len(cycles)):
                            if (not (set(cycles[c]) & set(cycles[a])) and
                                not (set(cycles[c]) & set(cycles[b]))):
                                print(f"    SIZE-3 IS: {cycles[a]}, {cycles[b]}, {cycles[c]}")
            break

    print(f"I(Ω,2) = 1+2|cycles|+4|VD|: {io_check}/{io_total}")

    if io_check < io_total:
        print("Formula fails at n=7! Need higher-order terms.")
        # Check: what's the max independent set size?
        # At n=7: two VD 3-cycles use 6 vertices, leaving 1 free.
        # A 3rd cycle needs at least 3 vertices. Can't fit.
        # So max IS size should still be 2!
        # Unless there's a VD pair involving a 5-cycle...
        # 3+5=8>7. So no VD 3+5 pair at n=7.
        # So max IS = 2 at n=7. Why does formula fail?
        sys.exit(1)

    # Now test ΔI formula
    print(f"\nTesting ΔI = -2·Σ s_x·H(B_x) + 2·Σ(DL-CL)...")
    total = 0
    match = 0

    random.seed(123)
    for _ in range(100):
        T = random_tournament(n)
        # Pick a random arc to flip
        i, j = random.sample(range(n), 2)
        if T[i][j] == 0:
            i, j = j, i
        assert T[i][j] == 1

        total += 1
        T2 = flip_arc(T, i, j)
        others = [x for x in range(n) if x != i and x != j]

        # Formula
        sx_hbx_sum = 0
        for x in others:
            sx = 1 - T[x][i] - T[j][x]
            bx = [u for u in others if u != x]
            hbx = sub_H(T, bx)
            sx_hbx_sum += sx * hbx

        # Odd-cycle terms by length
        dc_T = count_odd_cycles_using_arc(T, n, i, j)
        dc_T2 = count_odd_cycles_using_arc(T2, n, j, i)

        cycle_term = 0
        for L in range(3, n + 1, 2):
            dL = dc_T.get(L, 0)
            cL = dc_T2.get(L, 0)
            cycle_term += dL - cL

        formula_dI = -2 * sx_hbx_sum + 2 * cycle_term

        # Ground truth
        ht = hamiltonian_path_count(T)
        ht2 = hamiltonian_path_count(T2)
        delta_H = ht - ht2

        if formula_dI == delta_H:
            match += 1
        elif total <= 5:
            print(f"  MISS: formula={formula_dI}, gt_ΔH={delta_H}")
            print(f"    sx_hbx_sum={sx_hbx_sum}, cycle_term={cycle_term}")
            print(f"    dc_T={dc_T}, dc_T2={dc_T2}")

    print(f"\nΔH = formula: {match}/{total}")

    if match == total:
        print("\n*** FORMULA VERIFIED AT n=7! ***")
        print("ΔI = -2·Σ_x s_x·H(B_x) + 2·Σ_L(DL-CL)")
        print("generalizes correctly from n=6 to n=7!")
    else:
        print(f"\nFormula fails at n=7 ({total - match} failures)")
