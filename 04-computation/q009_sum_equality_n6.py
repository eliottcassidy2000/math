"""
OPEN-Q-009: Test the sum-equality approach at n=6.

When flipping arc i->j to j->i (i,j != v), Claim A (by induction) requires:

  delta_H - delta_Hv = 2 * delta(sum mu(C))

where delta(sum mu) decomposes into:
  (a) destroyed cycles: C through v using arc i->j (exist in T, not T')
  (b) created cycles: C' through v using arc j->i (exist in T', not T)
  (c) persisting cycles with both i,j OUTSIDE V(C): mu may change
  (d) persisting cycles with i or j IN V(C)\{v}: mu unchanged (THM-012)

The T021 sum-equality approach:
  sum_{destroyed C} mu_T(C) - sum_{created C'} mu_{T'}(C') + sum_{persist,change} [mu_T - mu_{T'}] = delta(sum mu)

By induction, mu(C) = H(T[V\V(C)]). For (c), V\V(C) contains both i,j,
so the flip changes H(T[V\V(C)]) by the adjacency formula on the complement.

Key test: Can we express delta(sum mu) purely in terms of
adjacency counts in sub-tournaments?

Author: opus-2026-03-05-S2
"""

import sys
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (
    all_tournaments, hamiltonian_path_count, delete_vertex,
    find_odd_cycles, conflict_graph, independence_poly_at, mu
)
from itertools import permutations


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def sub_tournament(T, verts):
    """Extract sub-tournament on given vertices (returns matrix + label map)."""
    verts = sorted(verts)
    n = len(verts)
    S = [[0]*n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            if a != b:
                S[a][b] = T[verts[a]][verts[b]]
    return S, verts


def adj_count_sub(T, verts, i, j):
    """Count Ham paths of T[verts] where i immediately precedes j.
    i, j are original vertex labels."""
    verts = sorted(verts)
    if i not in verts or j not in verts:
        return 0
    n = len(verts)
    count = 0
    for perm in permutations(verts):
        valid = True
        for k in range(len(perm) - 1):
            if T[perm[k]][perm[k+1]] != 1:
                valid = False
                break
        if valid:
            for k in range(len(perm) - 1):
                if perm[k] == i and perm[k+1] == j:
                    count += 1
    return count


if __name__ == "__main__":
    n = 6
    print(f"Q-009 Sum-Equality Decomposition at n={n}")
    print("="*70)

    # For each (T, v, flip i->j), decompose delta(sum mu) into components
    total = 0
    decomp_match = 0
    adj_formula_match = 0

    # Track: does the "complement adjacency" formula work?
    # For persist+change cycles C (both i,j in complement V\V(C)):
    #   mu_T(C) - mu_{T'}(C) = H(T[compl]) - H(T'[compl])
    #                        = adj_{compl}(i,j) - adj'_{compl}(j,i)
    complement_adj_ok = 0
    complement_adj_total = 0

    # Track the full decomposition
    results = []

    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= 300:  # sample
            break
        for v in range(n):
            others = [u for u in range(n) if u != v]
            for idx_i in range(len(others)):
                i = others[idx_i]
                for idx_j in range(idx_i + 1, len(others)):
                    j = others[idx_j]
                    if T[i][j] == 0:
                        continue  # only flip existing arcs
                    total += 1

                    T2 = flip_arc(T, i, j)

                    # Full delta_H and delta_Hv
                    ht = hamiltonian_path_count(T)
                    ht2 = hamiltonian_path_count(T2)
                    delta_H = ht - ht2

                    Tv_mat, _ = delete_vertex(T, v)
                    Tv2_mat, _ = delete_vertex(T2, v)
                    htv = hamiltonian_path_count(Tv_mat)
                    htv2 = hamiltonian_path_count(Tv2_mat)
                    delta_Hv = htv - htv2

                    # Cycles through v
                    cycles_T = find_odd_cycles(T)
                    cycles_T_v = [c for c in cycles_T if v in set(c)]
                    cycles_T2 = find_odd_cycles(T2)
                    cycles_T2_v = [c for c in cycles_T2 if v in set(c)]

                    # Classify: destroyed, created, persisting
                    set_T = set(tuple(c) for c in cycles_T_v)
                    set_T2 = set(tuple(c) for c in cycles_T2_v)

                    destroyed = set_T - set_T2
                    created = set_T2 - set_T
                    persisted = set_T & set_T2

                    # Compute mu for each category
                    Tv_d, ol_d = delete_vertex(T, v)
                    tc_d = find_odd_cycles(Tv_d)
                    cache_T = (Tv_d, ol_d, tc_d)

                    Tv2_d, ol2_d = delete_vertex(T2, v)
                    tc2_d = find_odd_cycles(Tv2_d)
                    cache_T2 = (Tv2_d, ol2_d, tc2_d)

                    sum_destroyed = sum(mu(T, v, c, cache_T) for c in destroyed)
                    sum_created = sum(mu(T2, v, c, cache_T2) for c in created)

                    # Persisting: split by whether both i,j outside V(C)
                    sum_persist_same = 0  # at least one of i,j in V(C)\{v}
                    sum_persist_delta = 0  # both i,j outside V(C)
                    n_persist_change = 0

                    for c in persisted:
                        vset = set(c)
                        mu_t = mu(T, v, c, cache_T)
                        mu_t2 = mu(T2, v, c, cache_T2)

                        if i not in vset and j not in vset:
                            # Both outside: mu may change
                            delta_mu_c = mu_t - mu_t2
                            sum_persist_delta += delta_mu_c

                            # Test: does delta_mu_c equal adjacency formula on complement?
                            compl_verts = [u for u in range(n) if u not in vset]
                            adj_ij_compl = adj_count_sub(T, compl_verts, i, j)
                            adj_ji_compl = adj_count_sub(T2, compl_verts, j, i)
                            adj_delta_compl = adj_ij_compl - adj_ji_compl

                            complement_adj_total += 1
                            if delta_mu_c == adj_delta_compl:
                                complement_adj_ok += 1
                            elif total <= 50:
                                print(f"  Complement adj MISMATCH: T#{t_idx} v={v} "
                                      f"flip {i}->{j} C={c}: "
                                      f"delta_mu={delta_mu_c} adj_delta={adj_delta_compl}")

                            if delta_mu_c != 0:
                                n_persist_change += 1
                        else:
                            # THM-012: mu unchanged
                            assert mu_t == mu_t2, f"THM-012 violation! {c}"
                            sum_persist_same += 0  # no contribution to delta

                    # Total delta(sum mu)
                    delta_sum_mu = sum_destroyed - sum_created + sum_persist_delta

                    # Check: delta_H - delta_Hv = 2 * delta_sum_mu
                    lhs = delta_H - delta_Hv
                    rhs = 2 * delta_sum_mu
                    if lhs == rhs:
                        decomp_match += 1

                    results.append({
                        'delta_H': delta_H, 'delta_Hv': delta_Hv,
                        'sum_destroyed': sum_destroyed,
                        'sum_created': sum_created,
                        'sum_persist_delta': sum_persist_delta,
                        'n_destroyed': len(destroyed),
                        'n_created': len(created),
                        'n_persist_change': n_persist_change,
                        'lhs': lhs, 'rhs': rhs,
                    })

    print(f"\nTotal flips analyzed: {total}")
    print(f"delta_H - delta_Hv = 2*delta(sum mu): {decomp_match}/{total}")

    print(f"\nComplement adjacency formula for persist+change cycles:")
    print(f"  mu_T(C) - mu_T'(C) = adj_compl(i,j) - adj'_compl(j,i): "
          f"{complement_adj_ok}/{complement_adj_total}")

    # Statistics
    n_with_persist_change = sum(1 for r in results if r['n_persist_change'] > 0)
    n_with_destroyed = sum(1 for r in results if r['n_destroyed'] > 0)
    n_with_created = sum(1 for r in results if r['n_created'] > 0)

    print(f"\nFlips with destroyed cycles: {n_with_destroyed}/{total}")
    print(f"Flips with created cycles: {n_with_created}/{total}")
    print(f"Flips with persisting-mu-changed cycles: {n_with_persist_change}/{total}")

    # Key question: is sum_destroyed always = sum_created
    # (with persist_delta = 0)?
    destroyed_eq_created = sum(1 for r in results
                                if r['sum_destroyed'] == r['sum_created']
                                and r['sum_persist_delta'] == 0)
    print(f"\nsum_destroyed = sum_created AND no persist change: "
          f"{destroyed_eq_created}/{total}")

    # Or is the relationship more nuanced?
    # Check: sum_destroyed - sum_created = -sum_persist_delta?
    balance = sum(1 for r in results
                  if r['sum_destroyed'] - r['sum_created'] == -r['sum_persist_delta'])
    print(f"sum_destroyed - sum_created = -sum_persist_delta: {balance}/{total}")
    # This should always hold (it's the definition of delta_sum_mu)

    # More interesting: what values does each component take?
    from collections import Counter
    d_vals = Counter(r['sum_destroyed'] for r in results)
    c_vals = Counter(r['sum_created'] for r in results)
    p_vals = Counter(r['sum_persist_delta'] for r in results)
    print(f"\nsum_destroyed distribution: {dict(sorted(d_vals.items()))}")
    print(f"sum_created distribution: {dict(sorted(c_vals.items()))}")
    print(f"sum_persist_delta distribution: {dict(sorted(p_vals.items()))}")

    # The key structural question:
    # Can we express delta(sum mu) purely in terms of adjacency counts?
    # delta(sum mu) = sum_destroyed - sum_created + persist_delta
    #
    # By induction, mu(C) = H(T[V\V(C)]).
    # For destroyed C: mu_T(C) = H(T[V\V(C)]) where V\V(C) is the complement
    # For created C': mu_{T'}(C') = H(T'[V\V(C')])
    # For persist+change: delta = H(T[compl]) - H(T'[compl]) = adj formula
    #
    # So delta(sum mu) = sum_{destroyed} H(T[compl_C])
    #                   - sum_{created} H(T'[compl_C'])
    #                   + sum_{persist,change} [adj_compl(i,j) - adj'_compl(j,i)]
    #
    # And we need this = (delta_H - delta_Hv) / 2
    # = [adj_T(i,j) - adj_{T'}(j,i) - adj_{T-v}(i,j) + adj_{T'-v}(j,i)] / 2

    print("\n" + "="*70)
    print("DETAILED EXAMPLES")
    print("="*70)
    shown = 0
    for idx, r in enumerate(results):
        if r['n_destroyed'] > 0 and r['n_persist_change'] > 0 and shown < 5:
            shown += 1
            print(f"\nExample {shown}: delta_H={r['delta_H']}, delta_Hv={r['delta_Hv']}")
            print(f"  Destroyed: {r['n_destroyed']} cycles, sum_mu={r['sum_destroyed']}")
            print(f"  Created: {r['n_created']} cycles, sum_mu={r['sum_created']}")
            print(f"  Persist+change: {r['n_persist_change']} cycles, "
                  f"sum_delta={r['sum_persist_delta']}")
            print(f"  LHS={r['lhs']}, RHS={r['rhs']}, match={r['lhs']==r['rhs']}")
