"""
OPEN-Q-009: Investigation of the F formula for arc-reversal invariance.

Key hypothesis to test:
  F(T, v, i, j) := adj_T(i,j) - adj_{T-v}(i,j)
                 = 2 * sum_{C through v using arc i->j} mu_T(C)

If true, then:
  delta_H - delta_Hv = F(T,v,i,j) - F(T',v,j,i)
                     = 2 * [sum_{destroyed} mu_T - sum_{created} mu_{T'}]

This still leaves the "persisting cycle mu change" term:
  sum_{C persisting, i,j not in V(C)} (mu_T(C) - mu_{T'}(C))

We also need to understand this term.

Secondary hypothesis:
  sum_{C persisting, i,j not in V(C)} (mu_T(C) - mu_{T'}(C)) = 0  ???

Author: opus-2026-03-05-S1
"""

import sys
sys.path.insert(0, '03-artifacts/code')
from tournament_lib import (
    all_tournaments, hamiltonian_path_count, delete_vertex,
    find_odd_cycles, mu, tournament_from_bits
)
from itertools import permutations
import random


def flip_arc(T, i, j):
    """Return copy of T with arc i->j flipped to j->i."""
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def adj_count(T, i, j):
    """Count Ham paths of T where i immediately precedes j (i,j are indices in T)."""
    n = len(T)
    if i >= n or j >= n or T[i][j] != 1:
        return 0
    count = 0
    for perm in permutations(range(n)):
        valid = all(T[perm[k]][perm[k+1]] == 1 for k in range(n-1))
        if valid:
            for k in range(n-1):
                if perm[k] == i and perm[k+1] == j:
                    count += 1
    return count


def adj_count_relabeled(Tv, old_labels, i, j):
    """Count Ham paths of Tv where vertex i precedes vertex j, with relabeled indices."""
    if i not in old_labels or j not in old_labels:
        return 0
    ni = old_labels.index(i)
    nj = old_labels.index(j)
    return adj_count(Tv, ni, nj)


def compute_sum_mu_for_arc_cycles(T, v, i, j, arc_i_j=True):
    """
    Compute sum of mu_T(C) for cycles C through v that use arc i->j (if arc_i_j=True)
    or arc j->i (if arc_i_j=False, for use with flipped tournament).
    """
    all_cycles = find_odd_cycles(T)
    cycles_v = [c for c in all_cycles if v in set(c)]

    if arc_i_j:
        src, dst = i, j
    else:
        src, dst = j, i

    # Filter to cycles using arc src->dst
    arc_cycles = []
    for c in cycles_v:
        verts = list(c)
        L = len(verts)
        for k in range(L):
            if verts[k] == src and verts[(k+1) % L] == dst:
                arc_cycles.append(c)
                break

    if not arc_cycles:
        return 0, []

    Tv, old_labels = delete_vertex(T, v)
    tv_cycles = find_odd_cycles(Tv)
    cache = (Tv, old_labels, tv_cycles)

    total = sum(mu(T, v, c, _tv_cache=cache) for c in arc_cycles)
    return total, arc_cycles


def compute_persist_mu_change(T, T2, v, i, j):
    """
    Compute sum over persisting cycles C (present in both T and T', through v,
    not using arc i->j or j->i, and with i,j both NOT in V(C)) of (mu_T(C) - mu_{T'}(C)).

    Also reports: total sum over all persisting cycles (including those with
    endpoints in V(C)) to verify THM-012.
    """
    cycles_T = {tuple(c): c for c in find_odd_cycles(T) if v in set(c)}
    cycles_T2 = {tuple(c): c for c in find_odd_cycles(T2) if v in set(c)}

    persist_keys = set(cycles_T.keys()) & set(cycles_T2.keys())

    Tv, ol = delete_vertex(T, v)
    tc = find_odd_cycles(Tv)
    Tv2, ol2 = delete_vertex(T2, v)
    tc2 = find_odd_cycles(Tv2)

    total_both_outside = 0
    total_one_inside = 0
    total_both_inside = 0

    for key in persist_keys:
        c = cycles_T[key]
        verts = set(c)

        # Check i,j relationship to V(C)
        i_in = (i in verts)
        j_in = (j in verts)
        # v is always in verts

        mu_T = mu(T, v, c, (Tv, ol, tc))
        mu_T2 = mu(T2, v, c, (Tv2, ol2, tc2))
        diff = mu_T - mu_T2

        if not i_in and not j_in:
            total_both_outside += diff
        elif i_in and j_in:
            # THM-012 Case 1: both in C\{v}, mu should not change
            if diff != 0:
                print(f"  THM-012 VIOLATION Case 1: i,j both in C\\{{v}}, diff={diff}")
            total_both_inside += diff
        else:
            # exactly one in C\{v}, THM-012 Cases 2,3: mu should not change
            if diff != 0:
                print(f"  THM-012 VIOLATION Cases 2,3: exactly one endpoint in C\\{{v}}, diff={diff}")
            total_one_inside += diff

    return total_both_outside, total_one_inside, total_both_inside


def run_analysis(n, num_samples=None, seed=42):
    """
    For each (T, v, arc i->j with i,j != v):
    - Compute F = adj_T(i,j) - adj_{T-v}(i,j)
    - Compute 2 * sum_{C through v using i->j} mu_T(C)
    - Check if F = 2 * (arc mu sum)
    - Also compute the persist_change term
    - Check: F - F' = 2 * delta_mu (where F' uses T' and j->i)
    """
    random.seed(seed)

    print(f"n={n} Analysis: Testing F formula")
    print("="*70)

    formula_holds = 0
    formula_fails = 0

    decomp_holds = 0
    decomp_fails = 0

    persist_nonzero = 0
    persist_zero = 0

    total = 0

    if n <= 5:
        tournaments_iter = list(all_tournaments(n))
        if num_samples and len(tournaments_iter) > num_samples:
            tournaments_iter = random.sample(tournaments_iter, num_samples)
    else:
        from tournament_lib import random_tournament
        tournaments_iter = [random_tournament(n) for _ in range(num_samples or 50)]

    fail_examples = []

    for T in tournaments_iter:
        for v in range(n):
            others = [u for u in range(n) if u != v]
            for i in others:
                for j in others:
                    if i == j or T[i][j] != 1:
                        continue
                    total += 1

                    T2 = flip_arc(T, i, j)

                    # Compute F(T, v, i, j)
                    a_ij = adj_count(T, i, j)
                    Tv, ol = delete_vertex(T, v)
                    av_ij = adj_count_relabeled(Tv, ol, i, j)
                    F = a_ij - av_ij

                    # Compute 2 * sum mu for arc-cycles in T
                    mu_arc_T, arc_cycles_T = compute_sum_mu_for_arc_cycles(T, v, i, j, True)
                    two_mu_arc = 2 * mu_arc_T

                    if F == two_mu_arc:
                        formula_holds += 1
                    else:
                        formula_fails += 1
                        if len(fail_examples) < 3:
                            fail_examples.append({
                                'n': n, 'v': v, 'i': i, 'j': j,
                                'F': F, '2*mu_arc': two_mu_arc,
                                'a_ij': a_ij, 'av_ij': av_ij,
                                'mu_arc_T': mu_arc_T,
                                'num_arc_cycles': len(arc_cycles_T)
                            })

                    # Compute the full delta_H - delta_Hv = 2*delta_mu check
                    a_ji2 = adj_count(T2, j, i)
                    Tv2, ol2 = delete_vertex(T2, v)
                    av_ji2 = adj_count_relabeled(Tv2, ol2, j, i)
                    F2 = a_ji2 - av_ji2

                    # delta_H - delta_Hv
                    lhs = F - F2

                    # delta_mu = sum_destroyed - sum_created + sum_persist_change
                    mu_arc_T2, _ = compute_sum_mu_for_arc_cycles(T2, v, i, j, False)  # created
                    destroyed_minus_created = mu_arc_T - mu_arc_T2

                    persist_change, persist_one, persist_both = compute_persist_mu_change(T, T2, v, i, j)

                    delta_mu = destroyed_minus_created + persist_change
                    rhs = 2 * delta_mu

                    if lhs == rhs:
                        decomp_holds += 1
                    else:
                        decomp_fails += 1

                    if persist_change != 0:
                        persist_nonzero += 1
                    else:
                        persist_zero += 1

    print(f"Total (T, v, arc) triples: {total}")
    print()
    print(f"F = 2*mu_arc formula: {formula_holds}/{total} holds ({100*formula_holds/total:.1f}%)")
    if formula_fails > 0:
        print(f"FORMULA FAILS in {formula_fails} cases. Examples:")
        for ex in fail_examples:
            print(f"  v={ex['v']}, {ex['i']}->{ex['j']}: F={ex['F']}, 2*mu={ex['2*mu_arc']}")
            print(f"    adj_T(i,j)={ex['a_ij']}, adj_Tv(i,j)={ex['av_ij']}")
            print(f"    mu_arc_sum={ex['mu_arc_T']}, num_arc_cycles={ex['num_arc_cycles']}")
    else:
        print("  FORMULA HOLDS IN ALL CASES")
    print()
    print(f"delta_H - delta_Hv = 2*delta_mu: {decomp_holds}/{total} holds")
    print()
    print(f"Persist cycle mu change nonzero: {persist_nonzero}/{total}")
    print(f"Persist cycle mu change = 0: {persist_zero}/{total}")

    # Summary insight
    print()
    print("="*70)
    if formula_holds == total:
        print("RESULT: F formula CONFIRMED: adj_T(i,j) - adj_{T-v}(i,j) = 2 * sum_{C through v using i->j} mu(C)")
        print()
        print("COROLLARY: delta_H - delta_Hv = F - F' = 2*(sum_destroyed - sum_created)")
        if persist_nonzero > 0:
            print(f"BUT: persist_change is nonzero in {persist_nonzero} cases.")
            print("=> Need a separate argument that persist_change cancels with something.")
        else:
            print("AND: persist_change = 0 always => F-F' = 2*delta_mu exactly.")
            print("This would give a COMPLETE proof strategy!")
    else:
        print("RESULT: F formula FAILS. Need different approach.")


if __name__ == "__main__":
    # n=5: exhaustive (small)
    run_analysis(n=5)

    print()

    # n=6: sample
    print("n=6 analysis (50 random tournaments):")
    run_analysis(n=6, num_samples=50, seed=123)
