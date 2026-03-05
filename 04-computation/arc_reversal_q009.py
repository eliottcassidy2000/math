"""
Investigate OPEN-Q-009: Arc-reversal invariance of D(T,v).

D(T,v) = H(T) - H(T-v) - 2 * sum_{C through v} mu(C)

Claim A says D(T,v) = 0 for all (T,v). An approach to proving this:
if D(T,v) is invariant under arc flips not involving v, then since
D(transitive, v) = 0 (no cycles), D = 0 for all T.

Define: T' = T with arc i->j flipped to j->i (where i,j != v).
Question: Is D(T,v) = D(T',v)?

Since D = 0 always (Claim A verified), this is trivially true.
But the interesting question is whether the COMPONENTS change:
- Does H(T) - H(T') have a nice formula?
- Does sum_{C through v} mu_T(C) - sum_{C through v} mu_{T'}(C) match?

If we can show these changes are equal WITHOUT assuming Claim A,
we get a proof of Claim A.

Author: opus-2026-03-05-S2
"""

import sys
sys.path.insert(0, '03-artifacts/code')
from tournament_lib import (
    all_tournaments, hamiltonian_path_count, delete_vertex,
    find_odd_cycles, conflict_graph, independence_poly_at, mu,
    tournament_from_bits
)
from itertools import permutations
import copy


def flip_arc(T, i, j):
    """Return a copy of T with arc i->j flipped to j->i."""
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def compute_rhs(T, v):
    """Compute 2 * sum_{C through v} mu(C)."""
    n = len(T)
    all_cycles = find_odd_cycles(T)
    cycles_v = [c for c in all_cycles if v in set(c)]
    if not cycles_v:
        return 0
    Tv, old_labels = delete_vertex(T, v)
    tv_cycles = find_odd_cycles(Tv)
    cache = (Tv, old_labels, tv_cycles)
    return 2 * sum(mu(T, v, c, _tv_cache=cache) for c in cycles_v)


def adj_count(T, i, j, excluded_vertices=None):
    """Count Hamiltonian paths of T (excluding some vertices) where i immediately precedes j."""
    n = len(T)
    if excluded_vertices is None:
        excluded_vertices = set()
    verts = [u for u in range(n) if u not in excluded_vertices]
    count = 0
    for perm in permutations(verts):
        valid = True
        for k in range(len(perm) - 1):
            if T[perm[k]][perm[k + 1]] != 1:
                valid = False
                break
        if valid:
            for k in range(len(perm) - 1):
                if perm[k] == i and perm[k + 1] == j:
                    count += 1
    return count


if __name__ == "__main__":
    print("="*70)
    print("OPEN-Q-009: Arc-reversal invariance investigation")
    print("="*70)

    # For each (T, v, arc i->j with i,j != v):
    # Compute delta_H = H(T) - H(T')
    # Compute delta_RHS = RHS(T,v) - RHS(T',v)
    # Check: delta_H = delta_RHS (which would follow from Claim A)
    # Also: analyze the structure of delta_H

    n = 5
    print(f"\n--- n={n}: Analyzing arc flip effects ---")

    delta_h_values = {}
    delta_rhs_values = {}
    adj_formula_works = 0
    total_flips = 0

    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= 200:  # sample
            break
        for v in range(n):
            others = [u for u in range(n) if u != v]
            for idx_i, i in enumerate(others):
                for j in others[idx_i + 1:]:
                    if T[i][j] == 0:
                        continue  # only flip existing arcs i->j
                    total_flips += 1

                    T2 = flip_arc(T, i, j)
                    ht = hamiltonian_path_count(T)
                    ht2 = hamiltonian_path_count(T2)
                    delta_h = ht - ht2

                    rhs_t = compute_rhs(T, v)
                    rhs_t2 = compute_rhs(T2, v)
                    delta_rhs = rhs_t - rhs_t2

                    # The adj formula: H(T) - H(T') = adj_T(i,j) - adj_T'(j,i)
                    # where adj_T(i,j) = #{Ham paths of T where i immediately before j}
                    adj_ij = adj_count(T, i, j)
                    adj_ji_prime = adj_count(T2, j, i)
                    adj_delta = adj_ij - adj_ji_prime

                    if delta_h == adj_delta:
                        adj_formula_works += 1

                    dh = delta_h
                    delta_h_values[dh] = delta_h_values.get(dh, 0) + 1

                    if delta_h != delta_rhs:
                        # This would be a Claim A failure
                        print(f"  CLAIM A INCONSISTENCY: T#{t_idx}, v={v}, "
                              f"flip {i}->{j}: delta_H={delta_h}, delta_RHS={delta_rhs}")

    print(f"\nTotal arc flips analyzed: {total_flips}")
    print(f"H(T)-H(T') = adj(i,j) - adj'(j,i): {adj_formula_works}/{total_flips}")
    print(f"delta_H distribution: {dict(sorted(delta_h_values.items()))}")

    # Now investigate: what's the formula for delta_RHS?
    # delta_RHS = 2[sum_{C through v in T} mu_T(C) - sum_{C through v in T'} mu_{T'}(C)]
    # = 2[sum_{C: i->j in C} mu_T(C) - sum_{C: j->i in C'} mu_{T'}(C')]
    #   + 2[sum_{C: arc(i,j) not in C} (mu_T(C) - mu_{T'}(C))]
    #
    # The first term: cycles that use the flipped arc directly.
    # The second term: cycles that DON'T use the arc but whose mu changes
    # because T-v changed (if i or j = v then T-v doesn't change, but
    # we assumed i,j != v, so T-v DOES change when we flip i<->j).

    print(f"\n--- Detailed decomposition of delta_RHS ---")

    n = 5
    examples = 0
    for t_idx, T in enumerate(all_tournaments(n)):
        if examples >= 10:
            break
        if t_idx >= 100:
            break
        for v in range(n):
            others = [u for u in range(n) if u != v]
            for i in others:
                for j in others:
                    if i == j or T[i][j] == 0:
                        continue
                    if examples >= 10:
                        break

                    T2 = flip_arc(T, i, j)

                    # Cycles through v in T
                    cycles_T = [c for c in find_odd_cycles(T) if v in set(c)]
                    cycles_T2 = [c for c in find_odd_cycles(T2) if v in set(c)]

                    # Separate: cycles using arc i->j vs not
                    uses_arc = []
                    no_arc = []
                    for c in cycles_T:
                        verts = list(c)
                        # Check if arc i->j is used in cycle c
                        found = False
                        for k in range(len(verts)):
                            if verts[k] == i and verts[(k + 1) % len(verts)] == j:
                                found = True
                                break
                        if found:
                            uses_arc.append(c)
                        else:
                            no_arc.append(c)

                    uses_arc_T2 = []
                    no_arc_T2 = []
                    for c in cycles_T2:
                        verts = list(c)
                        found = False
                        for k in range(len(verts)):
                            if verts[k] == j and verts[(k + 1) % len(verts)] == i:
                                found = True
                                break
                        if found:
                            uses_arc_T2.append(c)
                        else:
                            no_arc_T2.append(c)

                    if len(uses_arc) > 0 or len(uses_arc_T2) > 0:
                        examples += 1

                        Tv, old_labels = delete_vertex(T, v)
                        tv_cyc = find_odd_cycles(Tv)
                        cache_T = (Tv, old_labels, tv_cyc)

                        Tv2, old_labels2 = delete_vertex(T2, v)
                        tv_cyc2 = find_odd_cycles(Tv2)
                        cache_T2 = (Tv2, old_labels2, tv_cyc2)

                        mu_uses = sum(mu(T, v, c, cache_T) for c in uses_arc)
                        mu_no = sum(mu(T, v, c, cache_T) for c in no_arc)
                        mu_uses2 = sum(mu(T2, v, c, cache_T2) for c in uses_arc_T2)
                        mu_no2 = sum(mu(T2, v, c, cache_T2) for c in no_arc_T2)

                        ht = hamiltonian_path_count(T)
                        ht2 = hamiltonian_path_count(T2)

                        print(f"\nT#{t_idx}, v={v}, flip {i}->{j}:")
                        print(f"  H(T)={ht}, H(T')={ht2}, delta_H={ht-ht2}")
                        print(f"  Cycles using arc: T={len(uses_arc)}, T'={len(uses_arc_T2)}")
                        print(f"  Cycles not using arc: T={len(no_arc)}, T'={len(no_arc_T2)}")
                        print(f"  sum mu (uses arc): T={mu_uses}, T'={mu_uses2}")
                        print(f"  sum mu (no arc): T={mu_no}, T'={mu_no2}")
                        print(f"  delta_RHS = 2*({mu_uses}+{mu_no} - {mu_uses2}-{mu_no2})"
                              f" = {2*(mu_uses+mu_no-mu_uses2-mu_no2)}")
                        # Check vertex sets of cycles
                        vsets_T = sorted([frozenset(c) for c in uses_arc])
                        vsets_T2 = sorted([frozenset(c) for c in uses_arc_T2])
                        print(f"  Vertex sets (arc cycles): T={vsets_T}, T'={vsets_T2}")

    # Key question: for "no arc" cycles, does mu change?
    # These cycles exist in both T and T' (same directed cycle, arc i->j not used).
    # But mu(C) = I(Omega(T-v)|_{avoid C\{v}}, 2), and T-v CHANGED (flipped i<->j in T-v).
    # So mu CAN change.
    print(f"\n--- Does mu change for cycles NOT using the flipped arc? ---")
    mu_change_count = 0
    mu_same_count = 0
    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= 200:
            break
        for v in range(n):
            others = [u for u in range(n) if u != v]
            for i in others:
                for j in others:
                    if i == j or T[i][j] == 0:
                        continue
                    T2 = flip_arc(T, i, j)
                    cycles_T = [c for c in find_odd_cycles(T) if v in set(c)]
                    for c in cycles_T:
                        verts = list(c)
                        uses = any(verts[k] == i and verts[(k+1)%len(verts)] == j for k in range(len(verts)))
                        if uses:
                            continue
                        # This cycle exists in both T and T'. Check if mu changes.
                        # Verify cycle still valid in T'
                        still_valid = all(T2[verts[k]][verts[(k+1)%len(verts)]] == 1 for k in range(len(verts)))
                        if not still_valid:
                            continue  # cycle uses arc j->i which got flipped
                        Tv, ol = delete_vertex(T, v)
                        tc = find_odd_cycles(Tv)
                        mu_T = mu(T, v, c, (Tv, ol, tc))
                        Tv2, ol2 = delete_vertex(T2, v)
                        tc2 = find_odd_cycles(Tv2)
                        mu_T2 = mu(T2, v, c, (Tv2, ol2, tc2))
                        if mu_T != mu_T2:
                            mu_change_count += 1
                        else:
                            mu_same_count += 1
    print(f"  mu changed: {mu_change_count}, mu same: {mu_same_count}")
    if mu_change_count > 0:
        print(f"  mu DOES change for non-arc cycles ({mu_change_count} cases)")
        print(f"  This means the arc-flip approach must account for mu changes globally")
    else:
        print(f"  mu does NOT change for non-arc cycles")
        print(f"  This simplifies the arc-flip invariance proof!")
