"""Check: which component of ΔI formula is wrong?"""

import sys
sys.path.insert(0, '03-artifacts/code')
from tournament_lib import (
    all_tournaments, find_odd_cycles, conflict_graph, independence_poly_at
)
from itertools import permutations


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def is_cyclic_3(T, a, b, c):
    return (T[a][b] and T[b][c] and T[c][a]) or (T[a][c] and T[c][b] and T[b][a])


def count_5cycles_using_arc(T, n, i, j):
    count = 0
    for x in range(n):
        if x == i or x == j:
            continue
        others = [u for u in range(n) if u not in {i, j, x}]
        for perm in permutations(others):
            if (T[j][perm[0]] and T[perm[0]][perm[1]] and
                T[perm[1]][perm[2]] and T[perm[2]][i]):
                count += 1
    return count


n = 6
shown = 0
cycle_ok = 0
vd_ok = 0
total = 0

for t_idx, T in enumerate(all_tournaments(n)):
    if t_idx >= 200 or shown >= 5:
        break
    for i in range(n):
        for j in range(i + 1, n):
            if T[i][j] == 0 or shown >= 5:
                continue
            total += 1
            T2 = flip_arc(T, i, j)
            others = [x for x in range(n) if x != i and x != j]

            # Ground truth cycle counts
            cycles_T = find_odd_cycles(T)
            cycles_T2 = find_odd_cycles(T2)
            set_T = set(tuple(c) for c in cycles_T)
            set_T2 = set(tuple(c) for c in cycles_T2)

            gt_delta_cycles = len(cycles_T) - len(cycles_T2)

            # Formula: delta_cycles = -sum s_x + (D5 - C5)
            # where -sum_x s_x = D3 - C3
            s_sum = sum(1 - T[x][i] - T[j][x] for x in others)
            d5 = count_5cycles_using_arc(T, n, i, j)
            c5 = count_5cycles_using_arc(T2, n, j, i)
            formula_delta_cycles = -s_sum + (d5 - c5)

            if gt_delta_cycles == formula_delta_cycles:
                cycle_ok += 1
            elif shown < 5:
                shown += 1
                destroyed = set_T - set_T2
                created = set_T2 - set_T
                d3 = sum(1 for c in destroyed if len(c) == 3)
                c3 = sum(1 for c in created if len(c) == 3)
                d5_gt = sum(1 for c in destroyed if len(c) == 5)
                c5_gt = sum(1 for c in created if len(c) == 5)
                print(f"\nCYCLE MISMATCH T#{t_idx} flip {i}->{j}:")
                print(f"  GT: delta_cycles={gt_delta_cycles}, D={len(destroyed)}, C={len(created)}")
                print(f"  GT: D3={d3}, D5={d5_gt}, C3={c3}, C5={c5_gt}")
                print(f"  Formula: -s_sum={-s_sum}, D5={d5}, C5={c5}")
                print(f"  s_x values: {[1-T[x][i]-T[j][x] for x in others]}")
                # Check D3 formula specifically
                d3_formula = sum(T[j][x] * T[x][i] for x in others)
                c3_formula = sum(T[i][x] * T[x][j] for x in others)
                print(f"  D3_formula={d3_formula}, C3_formula={c3_formula}")
                print(f"  -s_sum should be D3-C3={d3_formula}-{c3_formula}={d3_formula-c3_formula}")

            # Ground truth VD pairs
            def count_vd(cycles):
                c = 0
                for a in range(len(cycles)):
                    for b in range(a + 1, len(cycles)):
                        if not (set(cycles[a]) & set(cycles[b])):
                            c += 1
                return c

            gt_delta_vd = count_vd(cycles_T) - count_vd(cycles_T2)
            formula_delta_vd = sum((1 - T[x][i] - T[j][x]) *
                                    (1 if is_cyclic_3(T, *[u for u in others if u != x]) else 0)
                                    for x in others)

            if gt_delta_vd == formula_delta_vd:
                vd_ok += 1

print(f"\nTotal: {total}")
print(f"Δ(#cycles) formula correct: {cycle_ok}/{total}")
print(f"Δ(#VD pairs) formula correct: {vd_ok}/{total}")
