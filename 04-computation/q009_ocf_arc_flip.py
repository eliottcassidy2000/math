"""
Q-009 via OCF: Show ΔH = ΔI(Ω(T), 2) for each arc flip.

Since OCF (H(T) = I(Ω(T), 2)) holds at the transitive base case, and
any tournament is reachable by arc flips from transitive, proving
ΔH = ΔI for each flip proves OCF for all tournaments.

For a flip i→j → j→i:
  ΔH = adj_T(i,j) - adj_{T'}(j,i)

  ΔI = I(Ω(T), 2) - I(Ω(T'), 2)

Decompose ΔI using the destroyed/created/persisting cycle structure:
  Ω(T) = P ∪ D, Ω(T') = P ∪ C
where D = destroyed cycles, C = created cycles, P = persisting cycles.

Key insight: for each independent set S of P, the contribution from D vs C is:
  f_D(S_P) = sum_{S_D indep in D, no conflict with S_P} 2^|S_D|
  f_C(S_P) = sum_{S_C indep in C, no conflict with S_P} 2^|S_C|
  ΔI = sum_{S_P indep in P} 2^|S_P| * [f_D(S_P) - f_C(S_P)]

Author: opus-2026-03-05-S2
"""

import sys
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (
    all_tournaments, hamiltonian_path_count,
    find_odd_cycles, conflict_graph, independence_poly_at
)
from itertools import permutations


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def adj_count(T, i, j):
    """Count Ham paths of T where i immediately precedes j."""
    n = len(T)
    count = 0
    for perm in permutations(range(n)):
        valid = True
        for k in range(n - 1):
            if T[perm[k]][perm[k+1]] != 1:
                valid = False
                break
        if valid:
            for k in range(n - 1):
                if perm[k] == i and perm[k+1] == j:
                    count += 1
    return count


def indep_sets(cycles, cg_adj):
    """Enumerate all independent sets in the conflict graph.
    Returns list of (set of cycle indices, size)."""
    n = len(cycles)
    result = [frozenset()]
    for idx in range(n):
        new_sets = []
        for s in result:
            # Can we add idx to s?
            can_add = all(cg_adj[idx][j] == 0 for j in s)
            if can_add:
                new_sets.append(s | {idx})
        result.extend(new_sets)
    return result


def compute_delta_I_decomposed(T, T2, i, j):
    """Compute ΔI = I(Ω(T),2) - I(Ω(T'),2) with full decomposition."""
    cycles_T = find_odd_cycles(T)
    cycles_T2 = find_odd_cycles(T2)

    set_T = set(tuple(c) for c in cycles_T)
    set_T2 = set(tuple(c) for c in cycles_T2)

    destroyed = set_T - set_T2  # in T, not in T'
    created = set_T2 - set_T    # in T', not in T
    persisted = set_T & set_T2  # in both

    # Index cycles
    persist_list = sorted(persisted)
    destroy_list = sorted(destroyed)
    create_list = sorted(created)

    # Build conflict graphs
    all_T = persist_list + destroy_list
    all_T2 = persist_list + create_list

    def conflicts(c1, c2):
        return len(set(c1) & set(c2)) > 0

    n_p = len(persist_list)
    n_d = len(destroy_list)
    n_c = len(create_list)

    # For Ω(T): persist + destroy
    n_T = n_p + n_d
    adj_T = [[0]*n_T for _ in range(n_T)]
    for a in range(n_T):
        for b in range(a+1, n_T):
            if conflicts(all_T[a], all_T[b]):
                adj_T[a][b] = adj_T[b][a] = 1

    # For Ω(T'): persist + create
    n_T2 = n_p + n_c
    adj_T2 = [[0]*n_T2 for _ in range(n_T2)]
    for a in range(n_T2):
        for b in range(a+1, n_T2):
            if conflicts(all_T2[a], all_T2[b]):
                adj_T2[a][b] = adj_T2[b][a] = 1

    # Compute I(Ω(T), 2) and I(Ω(T'), 2)
    I_T = sum(2**len(s) for s in indep_sets(all_T, adj_T))
    I_T2 = sum(2**len(s) for s in indep_sets(all_T2, adj_T2))

    return {
        'I_T': I_T, 'I_T2': I_T2, 'delta_I': I_T - I_T2,
        'n_persist': n_p, 'n_destroyed': n_d, 'n_created': n_c,
    }


if __name__ == "__main__":
    print("Q-009 via OCF: ΔH = ΔI(Ω,2) for arc flips")
    print("="*70)

    for n in [4, 5]:
        print(f"\n--- n={n} ---")
        total = 0
        match = 0
        delta_patterns = {}  # (n_destroyed, n_created) -> count

        limit = 64 if n == 4 else 512
        for t_idx, T in enumerate(all_tournaments(n)):
            if t_idx >= limit:
                break
            for i in range(n):
                for j in range(i+1, n):
                    if T[i][j] == 0:
                        continue
                    total += 1
                    T2 = flip_arc(T, i, j)

                    # ΔH
                    delta_H = adj_count(T, i, j) - adj_count(T2, j, i)

                    # ΔI
                    info = compute_delta_I_decomposed(T, T2, i, j)
                    delta_I = info['delta_I']

                    if delta_H == delta_I:
                        match += 1

                    key = (info['n_destroyed'], info['n_created'])
                    delta_patterns[key] = delta_patterns.get(key, 0) + 1

        print(f"  Total flips: {total}")
        print(f"  ΔH = ΔI: {match}/{total}")
        print(f"  (destroyed, created) patterns: {dict(sorted(delta_patterns.items()))}")

    # Now at n=5: deeper analysis
    print(f"\n{'='*70}")
    print(f"n=5: Detailed decomposition of ΔI")
    print(f"{'='*70}")

    n = 5
    # For n=5: all mu=1, so I(Ω,2) = sum 2^k * alpha_k
    # When cycles are destroyed/created, the conflict graph changes.
    # Can we express ΔI in terms of just the destroyed/created cycle counts
    # and their conflicts with persisting cycles?

    results = []
    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= 200:
            break
        for i in range(n):
            for j in range(i+1, n):
                if T[i][j] == 0:
                    continue
                T2 = flip_arc(T, i, j)

                delta_H = adj_count(T, i, j) - adj_count(T2, j, i)
                info = compute_delta_I_decomposed(T, T2, i, j)

                # Count cycles through specific vertices
                cycles_T = find_odd_cycles(T)
                cycles_T2 = find_odd_cycles(T2)
                set_T = set(tuple(c) for c in cycles_T)
                set_T2 = set(tuple(c) for c in cycles_T2)
                destroyed = set_T - set_T2
                created = set_T2 - set_T

                # How many destroyed/created cycles contain both i,j?
                dest_with_ij = sum(1 for c in destroyed if i in set(c) and j in set(c))
                creat_with_ij = sum(1 for c in created if i in set(c) and j in set(c))

                # Destroyed cycle lengths
                dest_lens = sorted([len(c) for c in destroyed])
                creat_lens = sorted([len(c) for c in created])

                results.append({
                    'delta_H': delta_H, 'delta_I': info['delta_I'],
                    'n_dest': info['n_destroyed'], 'n_creat': info['n_created'],
                    'n_persist': info['n_persist'],
                    'dest_lens': dest_lens, 'creat_lens': creat_lens,
                    'dest_with_ij': dest_with_ij, 'creat_with_ij': creat_with_ij,
                })

    # At n=5, all mu=1. The conflict graph has vertices = odd cycles.
    # Two cycles conflict iff they share a vertex.
    # At n=5: 3-cycles have 3 vertices, 5-cycles have 5 vertices.
    # A 5-cycle conflicts with everything (shares vertex with any other cycle).
    # Two 3-cycles on 5 vertices: they share at least 1 vertex (pigeon), so they conflict.
    # So at n=5, Ω(T) is a COMPLETE graph on the odd cycles!
    # Therefore I(Ω,2) = 1 + 2*|cycles| (only independent sets of size 0 and 1).

    # Verify:
    all_complete = True
    for t_idx, T in enumerate(all_tournaments(5)):
        if t_idx >= 50:
            break
        cycles = find_odd_cycles(T)
        if len(cycles) <= 1:
            continue
        for a in range(len(cycles)):
            for b in range(a+1, len(cycles)):
                if not (set(cycles[a]) & set(cycles[b])):
                    all_complete = False
                    print(f"  NOT COMPLETE: T#{t_idx}, {cycles[a]} disjoint from {cycles[b]}")
    print(f"\n  At n=5: Ω(T) is always complete? {all_complete}")

    if all_complete:
        print(f"  Therefore I(Ω,2) = 1 + 2*|cycles|")
        print(f"  ΔI = 2*(|cycles_T| - |cycles_T'|) = 2*(|destroyed| - |created|)")
        print(f"  And ΔH = adj(i,j) - adj'(j,i)")
        print(f"  So: adj(i,j) - adj'(j,i) = 2*(|destroyed| - |created|)")
        print(f"\n  Verifying...")

        ok = sum(1 for r in results if r['delta_H'] == 2 * (r['n_dest'] - r['n_creat']))
        print(f"  ΔH = 2*(#destroyed - #created): {ok}/{len(results)}")

        # This is a HUGE simplification at n=5!
        # adj(i,j) - adj'(j,i) = 2 * net destroyed cycles

    # At n=6: Ω(T) is NOT always complete (can have disjoint 3-cycles).
    print(f"\n{'='*70}")
    print(f"n=6: Is Ω(T) still complete?")
    print(f"{'='*70}")

    non_complete = 0
    for t_idx, T in enumerate(all_tournaments(6)):
        if t_idx >= 200:
            break
        cycles = find_odd_cycles(T)
        if len(cycles) <= 1:
            continue
        for a in range(len(cycles)):
            for b in range(a+1, len(cycles)):
                if not (set(cycles[a]) & set(cycles[b])):
                    non_complete += 1
                    if non_complete <= 2:
                        print(f"  DISJOINT: T#{t_idx}, {cycles[a]} and {cycles[b]}")
                    break
            if non_complete > 0:
                break
    print(f"  Non-complete conflict graphs found: {non_complete > 0}")
    if non_complete:
        print(f"  At n=6, Ω(T) can have independent sets of size > 1.")
        print(f"  The simple formula ΔI = 2*(#destroyed - #created) may NOT hold.")
