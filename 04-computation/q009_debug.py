"""
Debug: find where the manual decomposition disagrees with ground truth.
Compare the manual destroyed/created/persist computation against
tournament_lib's actual cycle finding + mu.

Author: opus-2026-03-05-S2
"""

import sys
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (
    all_tournaments, hamiltonian_path_count, delete_vertex,
    find_odd_cycles, mu
)
from itertools import combinations


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def is_cyclic_3(T, a, b, c):
    if T[a][b] and T[b][c] and T[c][a]:
        return True
    if T[a][c] and T[c][b] and T[b][a]:
        return True
    return False


if __name__ == "__main__":
    n = 6
    shown = 0

    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= 100 or shown >= 5:
            break
        for v in range(n):
            if shown >= 5:
                break
            others = [u for u in range(n) if u != v]
            for idx_i in range(len(others)):
                i = others[idx_i]
                for idx_j in range(idx_i + 1, len(others)):
                    j = others[idx_j]
                    if T[i][j] == 0 or shown >= 5:
                        continue
                    T2 = flip_arc(T, i, j)

                    # Ground truth
                    cycles_T = find_odd_cycles(T)
                    cycles_T_v = [c for c in cycles_T if v in set(c)]
                    cycles_T2 = find_odd_cycles(T2)
                    cycles_T2_v = [c for c in cycles_T2 if v in set(c)]

                    set_T_v = set(tuple(c) for c in cycles_T_v)
                    set_T2_v = set(tuple(c) for c in cycles_T2_v)

                    destroyed = set_T_v - set_T2_v
                    created = set_T2_v - set_T_v
                    persisted = set_T_v & set_T2_v

                    Tv_d, ol_d = delete_vertex(T, v)
                    tc_d = find_odd_cycles(Tv_d)
                    cache_T = (Tv_d, ol_d, tc_d)
                    Tv2_d, ol2_d = delete_vertex(T2, v)
                    tc2_d = find_odd_cycles(Tv2_d)
                    cache_T2 = (Tv2_d, ol2_d, tc2_d)

                    gt_sum_dest = sum(mu(T, v, c, cache_T) for c in destroyed)
                    gt_sum_creat = sum(mu(T2, v, c, cache_T2) for c in created)
                    gt_persist_delta = 0
                    for c in persisted:
                        vset = set(c)
                        if i not in vset and j not in vset:
                            m1 = mu(T, v, c, cache_T)
                            m2 = mu(T2, v, c, cache_T2)
                            gt_persist_delta += m1 - m2
                    gt_delta = gt_sum_dest - gt_sum_creat + gt_persist_delta

                    # Manual computation
                    compl = sorted(set(range(n)) - {v, i, j})
                    compl_H = 3 if is_cyclic_3(T, *compl) else 1

                    # 3-cycle (v,i,j): v->i->j->v
                    has_3_dest = T[v][i] == 1 and T[j][v] == 1
                    # 3-cycle (v,j,i) in T': v->j->i->v
                    has_3_creat = T[v][j] == 1 and T[i][v] == 1

                    man_3_dest = compl_H if has_3_dest else 0
                    man_3_creat = compl_H if has_3_creat else 0

                    # 5-cycles: count from ground truth
                    n_5_dest = sum(1 for c in destroyed if len(c) == 5)
                    n_5_creat = sum(1 for c in created if len(c) == 5)
                    n_3_dest = sum(1 for c in destroyed if len(c) == 3)
                    n_3_creat = sum(1 for c in created if len(c) == 3)

                    # Persist delta manual: 3-cycles (v, a2, b2) with a2,b2 in compl
                    man_persist = 0
                    for pair in combinations(compl, 2):
                        a2, b2 = pair
                        cyc_fwd = T[v][a2] and T[a2][b2] and T[b2][v]
                        cyc_rev = T[v][b2] and T[b2][a2] and T[a2][v]
                        if cyc_fwd or cyc_rev:
                            remaining = [x for x in compl if x not in {a2, b2}][0]
                            cc = sorted([i, j, remaining])
                            mu_T_val = 3 if is_cyclic_3(T, *cc) else 1
                            mu_T2_val = 3 if is_cyclic_3(T2, *cc) else 1
                            man_persist += mu_T_val - mu_T2_val

                    man_delta = man_3_dest - man_3_creat + n_5_dest - n_5_creat + man_persist

                    if man_delta != gt_delta:
                        shown += 1
                        print(f"\nMISMATCH T#{t_idx} v={v} flip {i}->{j}:")
                        print(f"  GT: dest_mu={gt_sum_dest} creat_mu={gt_sum_creat} "
                              f"persist={gt_persist_delta} -> delta={gt_delta}")
                        print(f"  Manual: 3_dest={man_3_dest} 3_creat={man_3_creat} "
                              f"5_dest={n_5_dest} 5_creat={n_5_creat} "
                              f"persist={man_persist} -> delta={man_delta}")
                        print(f"  Destroyed cycles: {[(c, len(c)) for c in destroyed]}")
                        print(f"  Created cycles: {[(c, len(c)) for c in created]}")
                        print(f"  GT 3-dest mu: {[mu(T, v, c, cache_T) for c in destroyed if len(c)==3]}")
                        print(f"  GT 5-dest mu: {[mu(T, v, c, cache_T) for c in destroyed if len(c)==5]}")
                        print(f"  GT 3-creat mu: {[mu(T2, v, c, cache_T2) for c in created if len(c)==3]}")
                        print(f"  GT 5-creat mu: {[mu(T2, v, c, cache_T2) for c in created if len(c)==5]}")
                        print(f"  n_3_dest={n_3_dest}, has_3_dest={has_3_dest}, "
                              f"n_3_creat={n_3_creat}, has_3_creat={has_3_creat}")
                        print(f"  compl_H={compl_H}")
                        # Show destroyed 3-cycles
                        for c in destroyed:
                            if len(c) == 3:
                                print(f"    Destroyed 3-cycle: {c}, V(C)={set(c)}, "
                                      f"contains i={i in set(c)}, j={j in set(c)}")
                                compl_c = sorted(set(range(n)) - set(c))
                                print(f"    Complement: {compl_c}")
                                h_compl = 3 if is_cyclic_3(T, *compl_c) else 1
                                print(f"    H(complement)={h_compl}, mu={mu(T, v, c, cache_T)}")

    if shown == 0:
        print("No mismatches found in first 100 tournaments!")
