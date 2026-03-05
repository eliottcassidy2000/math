"""
Q-009: Attempt to express delta(sum mu) in a closed form at n=6.

The full equation (Claim A under arc flip i->j, i,j != v):
  adj_T^v(i,j) - adj_{T'}^v(j,i) = 2 * delta(sum mu)

where adj_T^v(i,j) = # of n-vertex Ham paths of T with i immediately before j.

delta(sum mu) decomposes as:
  D = sum_{destroyed C} mu(C) - sum_{created C'} mu(C') + persist_delta

For n=6:
- Destroyed/created 3-cycles: C = (v,i,j), mu = H(T[{a,b,c}]) where {a,b,c} = V\{v,i,j}
- Destroyed/created 5-cycles: C = (v,i,j,a,b), mu = 1
- Persist-change: 3-cycles C not containing i or j, mu change = adj(i,j) on 3-vertex complement

Key question: can we express D purely in terms of:
  - The 3-vertex complement T[{a,b,c}] structure (cyclic vs transitive)
  - The arc directions involving v
  - The destroyed/created 5-cycle counts

Author: opus-2026-03-05-S2
"""

import sys
sys.path.insert(0, '03-artifacts/code')
from tournament_lib import (
    all_tournaments, hamiltonian_path_count, delete_vertex,
    find_odd_cycles, mu
)
from itertools import permutations, combinations
from collections import Counter


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def is_cyclic_3(T, a, b, c):
    """Check if {a,b,c} form a directed 3-cycle in T."""
    # 3-cycle a->b->c->a or a->c->b->a
    if T[a][b] and T[b][c] and T[c][a]:
        return True
    if T[a][c] and T[c][b] and T[b][a]:
        return True
    return False


def count_5cycles_through_vij(T, v, i, j, n):
    """Count directed 5-cycles through v that contain both i and j
    and use the arc i->j."""
    other_verts = [u for u in range(n) if u not in {v, i, j}]
    count = 0
    for pair in combinations(other_verts, 2):
        a, b = pair
        verts = [v, i, j, a, b]
        # Enumerate all directed 5-cycles on these 5 vertices
        # that go through v and use arc i->j
        for perm in permutations(verts):
            if perm[0] != v:  # fix starting vertex to v
                continue
            # Check if it's a valid directed cycle using arc i->j
            valid = True
            uses_ij = False
            for k in range(5):
                u1 = perm[k]
                u2 = perm[(k + 1) % 5]
                if T[u1][u2] != 1:
                    valid = False
                    break
                if u1 == i and u2 == j:
                    uses_ij = True
            if valid and uses_ij:
                count += 1
    # Each 5-cycle is counted 5 times (one for each rotation starting at v)
    # Actually no: we fixed start = v, so each cycle is counted once
    # But wait, the cycle v->a->i->j->b->v is the same as going around
    # We fixed the start to v and enumerate orderings, so each cycle appears once
    return count


def adj_v(T, v, i, j, n):
    """Count n-vertex Hamiltonian paths of T with i immediately before j."""
    count = 0
    for perm in permutations(range(n)):
        valid = True
        for k in range(n - 1):
            if T[perm[k]][perm[k + 1]] != 1:
                valid = False
                break
        if valid:
            for k in range(n - 1):
                if perm[k] == i and perm[k + 1] == j:
                    count += 1
    return count


if __name__ == "__main__":
    n = 6
    print(f"Q-009 Closed-Form Analysis at n={n}")
    print("=" * 70)

    # For each (T, v, i->j flip), compute all terms
    results = []
    total = 0

    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= 200:
            break
        for v in range(n):
            others = [u for u in range(n) if u != v]
            for idx_i in range(len(others)):
                i = others[idx_i]
                for idx_j in range(idx_i + 1, len(others)):
                    j = others[idx_j]
                    if T[i][j] == 0:
                        continue
                    total += 1
                    T2 = flip_arc(T, i, j)

                    # Complement {a,b,c} = V\{v,i,j}
                    compl = sorted(set(range(n)) - {v, i, j})
                    a, b, c = compl

                    # Is complement cyclic?
                    compl_cyclic = is_cyclic_3(T, a, b, c)
                    compl_H = 3 if compl_cyclic else 1

                    # 3-cycle (v,i,j) exists? v->i->j->v
                    has_3cycle_vij = (T[v][i] == 1 and T[j][v] == 1)
                    # Note: T[i][j] == 1 is given

                    # Created 3-cycle (v,j,i) in T'? v->j->i->v
                    # Needs T'[v][j]=T[v][j]=1 and T'[i][v]=T[i][v]=1
                    has_3cycle_vji = (T[v][j] == 1 and T[i][v] == 1)

                    # Destroyed/created 5-cycle counts
                    # Destroyed: 5-cycles in T through v using i->j
                    n_destroyed_5 = count_5cycles_through_vij(T, v, i, j, n)
                    # Created: 5-cycles in T' through v using j->i
                    n_created_5 = count_5cycles_through_vij(T2, v, j, i, n)

                    # Persist-change: 3-cycles through v with {i,j} both outside
                    # These are 3-cycles (v, a', b') where a', b' ∈ compl
                    # complement of such cycle = {i, j, remaining vertex}
                    persist_delta = 0
                    for pair in combinations(compl, 2):
                        a2, b2 = pair
                        # Check if (v, a2, b2) is a 3-cycle in T
                        cycle_exists_T = (T[v][a2] and T[a2][b2] and T[b2][v]) or \
                                         (T[v][b2] and T[b2][a2] and T[a2][v])
                        cycle_exists_T2 = (T2[v][a2] and T2[a2][b2] and T2[b2][v]) or \
                                          (T2[v][b2] and T2[b2][a2] and T2[a2][v])

                        if cycle_exists_T and cycle_exists_T2:
                            # Persisting cycle. Its complement contains i and j.
                            remaining = [x for x in compl if x not in {a2, b2}][0]
                            compl_cycle = sorted([i, j, remaining])
                            # Check if complement is cyclic in T vs T'
                            cyc_T = is_cyclic_3(T, *compl_cycle)
                            cyc_T2 = is_cyclic_3(T2, *compl_cycle)
                            mu_T = 3 if cyc_T else 1
                            mu_T2 = 3 if cyc_T2 else 1
                            persist_delta += mu_T - mu_T2

                    # Also check 5-cycles through v not using i->j arc
                    # At n=6, 5-cycle through v with both i,j outside V(C):
                    # impossible since V(C) has 5 vertices and |V|=6,
                    # so only 1 vertex outside. Can't have both i,j outside.
                    # So NO persist-change 5-cycles at n=6.

                    # Compute delta(sum mu)
                    sum_destroyed_mu = 0
                    if has_3cycle_vij:
                        sum_destroyed_mu += compl_H
                    sum_destroyed_mu += n_destroyed_5  # each has mu=1

                    sum_created_mu = 0
                    if has_3cycle_vji:
                        # Complement is same {a,b,c}, check in T'
                        # But T and T' agree on {a,b,c} (no i,j there)
                        sum_created_mu += compl_H
                    sum_created_mu += n_created_5  # each has mu=1

                    delta_sum_mu = sum_destroyed_mu - sum_created_mu + persist_delta

                    # LHS: adj_T^v(i,j) - adj_{T'}^v(j,i)
                    adj_T_v = adj_v(T, v, i, j, n)
                    adj_T2_v = adj_v(T2, v, j, i, n)
                    lhs = adj_T_v - adj_T2_v

                    match = (lhs == 2 * delta_sum_mu)

                    results.append({
                        'lhs': lhs,
                        'delta_sum_mu': delta_sum_mu,
                        'match': match,
                        'compl_H': compl_H,
                        'has_3_dest': has_3cycle_vij,
                        'has_3_creat': has_3cycle_vji,
                        'n_5_dest': n_destroyed_5,
                        'n_5_creat': n_created_5,
                        'persist_delta': persist_delta,
                        'adj_T_v': adj_T_v,
                        'adj_T2_v': adj_T2_v,
                        # Arc pattern: v->i? v->j?
                        'v_beats_i': T[v][i] == 1,
                        'v_beats_j': T[v][j] == 1,
                        'j_beats_v': T[j][v] == 1,
                        'i_beats_v': T[i][v] == 1,
                    })

    matches = sum(1 for r in results if r['match'])
    print(f"Total: {total}, Matches: {matches}/{total}")

    if matches < total:
        print("\nFAILURES:")
        for r in results:
            if not r['match']:
                print(f"  LHS={r['lhs']}, 2*delta={2*r['delta_sum_mu']}")
                break

    # Analyze: what determines the destroyed/created 3-cycle pattern?
    # has_3_dest = T[v][i]=1 AND T[j][v]=1
    # has_3_creat = T[v][j]=1 AND T[i][v]=1
    # Note: T[v][i]+T[i][v]=1, T[v][j]+T[j][v]=1
    # So (v_beats_i, v_beats_j) determines the pattern

    print("\n--- Pattern by arc directions (v->i?, v->j?) ---")
    for vi in [True, False]:
        for vj in [True, False]:
            subset = [r for r in results
                      if r['v_beats_i'] == vi and r['v_beats_j'] == vj]
            if not subset:
                continue
            n_sub = len(subset)
            n_dest = sum(1 for r in subset if r['has_3_dest'])
            n_creat = sum(1 for r in subset if r['has_3_creat'])
            avg_5_dest = sum(r['n_5_dest'] for r in subset) / n_sub
            avg_5_creat = sum(r['n_5_creat'] for r in subset) / n_sub
            avg_persist = sum(r['persist_delta'] for r in subset) / n_sub
            avg_lhs = sum(r['lhs'] for r in subset) / n_sub

            dest_str = "3-cycle destroyed" if vi and not vj else ""
            creat_str = "3-cycle created" if not vi and vj else ""
            both_str = ""
            if vi and vj:
                # v->i, v->j: j->v=0, so no destroyed; i->v=0, so no created
                both_str = "no 3-cycle destroyed or created"
            if not vi and not vj:
                # i->v, j->v: v->i=0 so no dest needs v->i; v->j=0 so no creat needs v->j
                both_str = "both impossible"

            # Wait, let me recalculate:
            # Destroyed (v,i,j): v->i AND j->v. So vi=True AND vj=False (j->v iff v->j=False)
            # Created (v,j,i): v->j AND i->v. So vj=True AND vi=False
            desc = f"v->i={vi}, v->j={vj}"
            has_dest = vi and (not vj)  # j->v
            has_creat = (not vi) and vj  # i->v

            print(f"\n  {desc} (n={n_sub}):")
            print(f"    3-cycle dest possible: {has_dest}, actual: {n_dest}")
            print(f"    3-cycle creat possible: {has_creat}, actual: {n_creat}")
            print(f"    Avg 5-destroyed: {avg_5_dest:.2f}, 5-created: {avg_5_creat:.2f}")
            print(f"    Avg persist_delta: {avg_persist:.3f}")
            print(f"    Avg LHS: {avg_lhs:.2f}")

    # Key: when v->i and j->v (i.e., v->i=T, v->j=F):
    # - 3-cycle (v,i,j) is destroyed with mu = H(complement)
    # - No 3-cycle created
    # - persist_delta from complement changes
    # So: adj_T^v(i,j) - adj_{T'}^v(j,i) = 2*(H(compl) + n_5_dest - n_5_creat + persist_delta)

    print("\n--- Testing simplified formula ---")
    # For each case, test:
    # LHS = 2 * (delta_3cycle_mu + delta_5cycles + persist_delta)
    all_ok = sum(1 for r in results if r['match'])
    print(f"All match: {all_ok}/{total}")

    # Decompose: what fraction of delta_sum_mu comes from each component?
    from_3 = Counter()
    from_5 = Counter()
    from_p = Counter()
    for r in results:
        d3 = (r['compl_H'] if r['has_3_dest'] else 0) - \
             (r['compl_H'] if r['has_3_creat'] else 0)
        d5 = r['n_5_dest'] - r['n_5_creat']
        dp = r['persist_delta']
        from_3[d3] = from_3.get(d3, 0) + 1
        from_5[d5] = from_5.get(d5, 0) + 1
        from_p[dp] = from_p.get(dp, 0) + 1

    print(f"\n3-cycle contribution distribution: {dict(sorted(from_3.items()))}")
    print(f"5-cycle contribution distribution: {dict(sorted(from_5.items()))}")
    print(f"Persist contribution distribution: {dict(sorted(from_p.items()))}")
