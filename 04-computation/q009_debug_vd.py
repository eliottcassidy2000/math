"""Debug: why does Δ(#VD pairs) formula fail?"""

import sys
sys.path.insert(0, '03-artifacts/code')
from tournament_lib import all_tournaments, find_odd_cycles


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def is_cyclic_3(T, a, b, c):
    return (T[a][b] and T[b][c] and T[c][a]) or (T[a][c] and T[c][b] and T[b][a])


def count_vd(cycles):
    c = 0
    for a in range(len(cycles)):
        for b in range(a + 1, len(cycles)):
            if not (set(cycles[a]) & set(cycles[b])):
                c += 1
    return c


n = 6
shown = 0
for t_idx, T in enumerate(all_tournaments(n)):
    if t_idx >= 200 or shown >= 3:
        break
    for i in range(n):
        for j in range(i + 1, n):
            if T[i][j] == 0 or shown >= 3:
                continue
            T2 = flip_arc(T, i, j)
            others = [x for x in range(n) if x != i and x != j]

            cycles_T = find_odd_cycles(T)
            cycles_T2 = find_odd_cycles(T2)
            gt_delta_vd = count_vd(cycles_T) - count_vd(cycles_T2)

            formula_delta_vd = 0
            for x in others:
                sx = 1 - T[x][i] - T[j][x]
                bx = [u for u in others if u != x]
                cx = is_cyclic_3(T, *bx)
                formula_delta_vd += sx * (1 if cx else 0)

            if gt_delta_vd != formula_delta_vd:
                shown += 1
                print(f"\nVD MISMATCH T#{t_idx} flip {i}->{j}:")
                print(f"  GT: delta_VD = {gt_delta_vd}")
                print(f"  Formula: delta_VD = {formula_delta_vd}")

                # Show all VD pairs in T and T'
                vd_T = []
                for a in range(len(cycles_T)):
                    for b in range(a + 1, len(cycles_T)):
                        if not (set(cycles_T[a]) & set(cycles_T[b])):
                            vd_T.append((cycles_T[a], cycles_T[b]))
                vd_T2 = []
                for a in range(len(cycles_T2)):
                    for b in range(a + 1, len(cycles_T2)):
                        if not (set(cycles_T2[a]) & set(cycles_T2[b])):
                            vd_T2.append((cycles_T2[a], cycles_T2[b]))
                print(f"  VD pairs in T ({len(vd_T)}):")
                for p in vd_T:
                    print(f"    {p[0]} & {p[1]}, vsets={set(p[0])} & {set(p[1])}")
                print(f"  VD pairs in T' ({len(vd_T2)}):")
                for p in vd_T2:
                    print(f"    {p[0]} & {p[1]}, vsets={set(p[0])} & {set(p[1])}")

                # Show per-x analysis
                for x in others:
                    ax = {i, j, x}
                    bx = set(range(n)) - ax
                    sx = 1 - T[x][i] - T[j][x]
                    cx = is_cyclic_3(T, *sorted(bx))
                    ct = is_cyclic_3(T, i, j, x)
                    ct2 = is_cyclic_3(T2, i, j, x)
                    print(f"    x={x}: A={ax}, B={bx}")
                    print(f"      s_x={sx}, c_T(A)={ct}, c_T'(A)={ct2}, c(B)={cx}")
                    print(f"      VD_T={ct and cx}, VD_T'={ct2 and cx}, delta={int(ct2 and cx)-int(ct and cx)}")
                    print(f"      Formula contribution: s_x*c(B) = {sx * int(cx)}")
                    if (int(ct2 and cx) - int(ct and cx)) != sx * int(cx):
                        print(f"      *** MISMATCH for this x! ***")
