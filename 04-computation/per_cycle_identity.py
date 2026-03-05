#!/usr/bin/env python3
"""
Investigate the per-cycle identity:
For each directed 3-cycle C=(v,b,a) [v->b->a->v]:
  N_TypeII(C) + N_Orphan(C) = 2 * mu(C)

where:
  N_TypeII(C) = #{P' in Ham(T-v): (b,a) consecutive in P' [i.e., b at pos j, a at pos j+1,
                  with TypeII = s[j]=1 (v->b) and s[j+1]=0 (a->v)]}
  N_Orphan(C) = #{P in Ham(T): a->v->b appears as consecutive triple in P}
                [v between a (precedes v) and b (follows v), T[a][v]=1, T[v][b]=1, T[b][a]=1]

Key insight (from investigate_inshat.py):
  H(T) - H(T-v) = sum_{P'} #{TypeII(P')} + #{orphan Ham paths}
  Every orphan path is associated with exactly one 3-cycle (bridging)
  -> H(T)-H(T-v) = sum_{3-cycles C} [N_TypeII(C) + N_Orphan(C)]

If N_TypeII(C) + N_Orphan(C) = 2*mu(C) per cycle, then:
  Claim A follows for n where all cycles through v are 3-cycles (n <= 5).
  For n >= 6 with 5-cycles, this per-cycle identity must FAIL for some 3-cycles,
  compensating for 5-cycle mu contributions.

Session: opus-2026-03-05-S3
"""
import sys
sys.path.insert(0, '/home/e/Documents/claude/math/03-artifacts/code')
from tournament_lib import *
from itertools import permutations


def count_n_typeii(T, v, cycle_3):
    """
    Count N_TypeII(C) for 3-cycle C=(v,b,a) [v->b->a->v].
    = #{P' in Ham(T-v): (b,a) consecutive in P' at a TypeII position}
    TypeII: v->b (T[v][b]=1) and a->v (T[a][v]=1) -- which is automatic for a 3-cycle.
    """
    v_idx, b, a = cycle_3  # cycle: v->b->a->v
    assert T[v_idx][b] and T[b][a] and T[a][v_idx], f"Not a valid 3-cycle: {cycle_3}"
    Tv, old_labels = delete_vertex(T, v_idx)
    b_new = old_labels.index(b)
    a_new = old_labels.index(a)
    count = 0
    for perm in permutations(range(len(Tv))):
        valid = all(Tv[perm[i]][perm[i+1]] for i in range(len(perm)-1))
        if not valid:
            continue
        # Check if (b,a) is consecutive: b at pos j, a at pos j+1
        for j in range(len(perm)-1):
            if perm[j] == b_new and perm[j+1] == a_new:
                count += 1
    return count


def count_n_orphan(T, v, cycle_3):
    """
    Count N_Orphan(C) for 3-cycle C=(v,b,a) [v->b->a->v].
    Orphan condition: v is between a (before) and b (after) in Ham path:
       ...->a->v->b->...  where T[a][v]=1 (a->v ✓), T[v][b]=1 (v->b ✓), T[b][a]=1 (backward)
    Note: for cycle v->b->a->v, T[a][v]=1 and T[v][b]=1. ✓
    """
    v_idx, b, a = cycle_3  # cycle: v->b->a->v
    assert T[v_idx][b] and T[b][a] and T[a][v_idx], f"Not a valid 3-cycle: {cycle_3}"
    count = 0
    n = len(T)
    for perm in permutations(range(n)):
        valid = all(T[perm[i]][perm[i+1]] for i in range(n-1))
        if not valid:
            continue
        # Check if a->v->b appears consecutively: a at pos k-1, v at pos k, b at pos k+1
        pos_v = perm.index(v_idx)
        if 0 < pos_v < n-1:
            if perm[pos_v-1] == a and perm[pos_v+1] == b:
                count += 1
    return count


def check_per_cycle_identity_n4():
    """Exhaustive check of per-cycle identity at n=4."""
    print("=== Per-cycle identity check: n=4 exhaustive ===")
    total_pairs = 0
    identity_fails = 0
    for T in all_tournaments(4):
        for v in range(4):
            Tv, old_labels = delete_vertex(T, v)
            tv_cycles = find_odd_cycles(Tv)
            cache = (Tv, old_labels, tv_cycles)
            all_cyc = find_odd_cycles(T)
            cycles_v = [c for c in all_cyc if v in set(c) and len(c) == 3]
            for c in cycles_v:
                # canonical: c[0] is smallest, c = (min, ?, ?)
                # cycle c[0]->c[1]->c[2]->c[0], v is in c
                # Rewrite so v is first: rotate
                i = list(c).index(v)
                rot = tuple(c[(i+k) % 3] for k in range(3))
                # rot = (v, b, a): v->b->a->v
                v2, b, a = rot
                assert v2 == v
                assert T[v][b] and T[b][a] and T[a][v], f"Bad rotation"
                mu_c = mu(T, v, c, _tv_cache=cache)
                nt = count_n_typeii(T, v, rot)
                no = count_n_orphan(T, v, rot)
                total_pairs += 1
                if nt + no != 2 * mu_c:
                    identity_fails += 1
                    print(f"  FAIL: T={T}, v={v}, C={rot}, N_TypeII={nt}, N_Orphan={no}, mu={mu_c}, sum={nt+no} != 2*{mu_c}={2*mu_c}")
    print(f"  Total (T,v,C) triples: {total_pairs}")
    print(f"  Per-cycle identity failures: {identity_fails}")
    print()


def check_per_cycle_identity_n5():
    """Sample check at n=5."""
    import random
    random.seed(42)
    print("=== Per-cycle identity check: n=5 sample (100 random T) ===")
    total = 0
    fails = 0
    for _ in range(100):
        T = random_tournament(5)
        for v in range(5):
            Tv, old_labels = delete_vertex(T, v)
            tv_cycles = find_odd_cycles(Tv)
            cache = (Tv, old_labels, tv_cycles)
            all_cyc = find_odd_cycles(T)
            cycles_v3 = [c for c in all_cyc if v in set(c) and len(c) == 3]
            for c in cycles_v3:
                i = list(c).index(v)
                rot = tuple(c[(i+k) % 3] for k in range(3))
                v2, b, a = rot
                mu_c = mu(T, v, c, _tv_cache=cache)
                nt = count_n_typeii(T, v, rot)
                no = count_n_orphan(T, v, rot)
                total += 1
                if nt + no != 2 * mu_c:
                    fails += 1
    print(f"  Triples: {total}, Fails: {fails}")
    print()


def check_per_cycle_identity_n6_with_5cycles():
    """Check at n=6: per-cycle identity for 3-cycles when 5-cycles exist."""
    import random
    random.seed(42)
    print("=== Per-cycle identity check: n=6 with 5-cycles ===")
    total3 = 0
    fails3 = 0
    cases_mu3 = 0
    cases_mu1 = 0
    sum_excess = 0  # (nt+no) - 2*mu, should be negative when mu=1 in presence of 5-cycles
    has_5cycle_count = 0

    for _ in range(200):
        T = random_tournament(6)
        for v in range(6):
            all_cyc = find_odd_cycles(T)
            cycles_v = [c for c in all_cyc if v in set(c)]
            cycles_v3 = [c for c in cycles_v if len(c) == 3]
            cycles_v5 = [c for c in cycles_v if len(c) == 5]

            Tv, old_labels = delete_vertex(T, v)
            tv_cycles = find_odd_cycles(Tv)
            cache = (Tv, old_labels, tv_cycles)

            has_5 = len(cycles_v5) > 0
            if has_5:
                has_5cycle_count += 1

            for c in cycles_v3:
                i = list(c).index(v)
                rot = tuple(c[(i+k) % 3] for k in range(3))
                v2, b, a = rot
                mu_c = mu(T, v, c, _tv_cache=cache)
                nt = count_n_typeii(T, v, rot)
                no = count_n_orphan(T, v, rot)
                total3 += 1
                excess = (nt + no) - 2 * mu_c
                sum_excess += excess
                if nt + no != 2 * mu_c:
                    fails3 += 1
                if mu_c == 3:
                    cases_mu3 += 1
                else:
                    cases_mu1 += 1

    print(f"  (T,v) pairs: {200*6}")
    print(f"  (T,v) pairs with 5-cycles: {has_5cycle_count}")
    print(f"  3-cycle triples: {total3} (mu=1: {cases_mu1}, mu=3: {cases_mu3})")
    print(f"  Per-cycle identity failures [N_TypeII+N_Orphan != 2*mu]: {fails3}")
    print(f"  Sum of excess (N_TypeII+N_Orphan - 2*mu): {sum_excess}")
    print(f"  Mean excess per 3-cycle triple: {sum_excess/total3:.4f}")
    print()

    # For the per-cycle approach to work for n=6, we'd need sum of excess over all 3-cycles
    # to equal 2 * (sum of 5-cycle mu contributions).
    # Let's verify this decomposition.
    print("=== Verifying: Σ_{3-cycles} (N_TypeII+N_Orphan) = 2*Σ_{ALL cycles} mu ===")
    claim_a_check_fails = 0
    decomp_fails = 0
    for _ in range(50):
        T = random_tournament(6)
        for v in range(6):
            all_cyc = find_odd_cycles(T)
            cycles_v = [c for c in all_cyc if v in set(c)]
            Tv, old_labels = delete_vertex(T, v)
            tv_cycles = find_odd_cycles(Tv)
            cache = (Tv, old_labels, tv_cycles)

            # LHS of formula: Σ_{3-cycles} [N_TypeII + N_Orphan]
            # = H(T) - H(T-v)  (by the formula from investigate_inshat.py)
            ht = hamiltonian_path_count(T)
            htv = hamiltonian_path_count(Tv)
            diff = ht - htv

            # RHS: 2 * Σ_{ALL cycles through v} mu(C) [Claim A, verified]
            total_mu = sum(mu(T, v, c, _tv_cache=cache) for c in cycles_v)
            if diff != 2 * total_mu:
                claim_a_check_fails += 1

    print(f"  Claim A failures: {claim_a_check_fails}/300")
    print()


if __name__ == '__main__':
    check_per_cycle_identity_n4()
    check_per_cycle_identity_n5()
    check_per_cycle_identity_n6_with_5cycles()
