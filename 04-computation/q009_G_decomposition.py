"""
OPEN-Q-009: Deeper decomposition of F = adj_T(i,j) - adj_{T-v}(i,j).

From the formula:
  adj_T(i,j) = sum_{P': Ham of T-v, i->j consecutive} (insact(v,P') - b(P',i,j))

where b(P',i,j) = T[i][v] * T[v][j]  (1 if v can be inserted between i and j)

So:
  F = adj_T(i,j) - adj_{T-v}(i,j)
    = sum_{P': i->j} (insact(v,P') - b(P',i,j)) - #{P': i->j}
    = sum_{P': i->j} (insact(v,P') - b(P',i,j) - 1)

Define:
  G(T,v,i,j) = sum_{P': i->j} (insact(v,P') - 1)
  b_sum(T,v,i,j) = sum_{P': i->j} b(P',i,j) = #{P': i->j, T[i][v]=1, T[v][j]=1}

Then F = G - b_sum.

Now sum_{P' ALL} (insact(v,P') - 1) = H(T) - H(T-v) = 2 * sum_C mu(C)  [Claim A]

QUESTION: What is G(T,v,i,j) in terms of cycles?

Also investigate: what is the CORRECT formula for F?

If F = 2 * sum_{C through v using i->j} mu(C) - b_sum, then
  G = 2 * mu_arc_sum

Let's test this.

Author: opus-2026-03-05-S1
"""

import sys
sys.path.insert(0, '03-artifacts/code')
from tournament_lib import (
    all_tournaments, hamiltonian_path_count, delete_vertex,
    find_odd_cycles, mu
)
from itertools import permutations
import random


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def insact_count(T, v, P_prime):
    """
    Count valid insertion positions for v into P_prime (a Ham path of T-v).
    A position between P_prime[k-1] and P_prime[k] is valid if
    T[P_prime[k-1]][v]=1 and T[v][P_prime[k]]=1.
    Position 0 (before P_prime[0]) is valid if T[v][P_prime[0]]=1.
    Position len (after P_prime[-1]) is valid if T[P_prime[-1]][v]=1.
    """
    m = len(P_prime)
    count = 0
    # Position 0: v before P_prime[0]
    if T[v][P_prime[0]] == 1:
        count += 1
    # Interior positions
    for k in range(m - 1):
        if T[P_prime[k]][v] == 1 and T[v][P_prime[k+1]] == 1:
            count += 1
    # Position m: v after P_prime[-1]
    if T[P_prime[-1]][v] == 1:
        count += 1
    return count


def compute_G_and_bsum(T, v, i, j):
    """
    Compute G(T,v,i,j) = sum_{P': i->j} (insact(v,P') - 1)
    and b_sum = #{P': i->j, T[i][v]=1 and T[v][j]=1}
    """
    n = len(T)
    others = [u for u in range(n) if u != v]

    b_ok = (T[i][v] == 1 and T[v][j] == 1)  # can v go between i and j?

    G = 0
    b_sum = 0
    count_paths = 0

    for perm in permutations(others):
        # Check if perm is a valid Ham path of T-v
        valid = all(T[perm[k]][perm[k+1]] == 1 for k in range(len(perm)-1))
        if not valid:
            continue
        # Check if i->j is consecutive in perm
        has_ij = False
        for k in range(len(perm)-1):
            if perm[k] == i and perm[k+1] == j:
                has_ij = True
                break
        if not has_ij:
            continue
        count_paths += 1
        ins = insact_count(T, v, perm)
        G += ins - 1
        if b_ok:
            b_sum += 1

    return G, b_sum, count_paths


def compute_mu_arc_sum(T, v, i, j):
    """Sum of mu_T(C) for cycles C through v using arc i->j."""
    all_cycles = find_odd_cycles(T)
    cycles_v = [c for c in all_cycles if v in set(c)]
    arc_cycles = []
    for c in cycles_v:
        verts = list(c)
        L = len(verts)
        for k in range(L):
            if verts[k] == i and verts[(k+1) % L] == j:
                arc_cycles.append(c)
                break
    if not arc_cycles:
        return 0
    Tv, old_labels = delete_vertex(T, v)
    tv_cycles = find_odd_cycles(Tv)
    cache = (Tv, old_labels, tv_cycles)
    return sum(mu(T, v, c, _tv_cache=cache) for c in arc_cycles)


def run_G_analysis(n, num_samples=None, seed=42):
    """Test: G(T,v,i,j) = 2 * mu_arc_sum(T,v,i,j)?"""
    random.seed(seed)
    print(f"n={n}: G decomposition analysis")
    print("="*60)

    if n <= 5:
        tournaments_iter = list(all_tournaments(n))
        if num_samples:
            random.shuffle(tournaments_iter)
            tournaments_iter = tournaments_iter[:num_samples]
    else:
        from tournament_lib import random_tournament
        tournaments_iter = [random_tournament(n) for _ in range(num_samples or 30)]

    G_eq_2mu = 0
    G_ne_2mu = 0
    total = 0
    fail_examples = []

    for T in tournaments_iter:
        for v in range(n):
            others = [u for u in range(n) if u != v]
            for ii, i in enumerate(others):
                for j in others[ii+1:]:  # avoid double-counting arcs
                    # Check which direction the arc goes
                    if T[i][j] == 1:
                        src, dst = i, j
                    else:
                        src, dst = j, i

                    # Only count arc src->dst (not both directions)
                    G, b_sum, count_paths = compute_G_and_bsum(T, v, src, dst)
                    mu_arc = compute_mu_arc_sum(T, v, src, dst)

                    if count_paths == 0:
                        continue

                    total += 1

                    if G == 2 * mu_arc:
                        G_eq_2mu += 1
                    else:
                        G_ne_2mu += 1
                        if len(fail_examples) < 5:
                            fail_examples.append({
                                'v': v, 'src': src, 'dst': dst,
                                'G': G, '2mu': 2*mu_arc, 'mu': mu_arc,
                                'b_sum': b_sum, 'count_paths': count_paths,
                                'F': G - b_sum
                            })

    print(f"Total (T,v,arc) with paths: {total}")
    print(f"G = 2*mu_arc: {G_eq_2mu}/{total} ({100*G_eq_2mu/max(total,1):.1f}%)")
    if fail_examples:
        print(f"FAILS: {G_ne_2mu}. Examples:")
        for ex in fail_examples:
            print(f"  v={ex['v']}, arc {ex['src']}->{ex['dst']}: G={ex['G']}, "
                  f"2*mu={ex['2mu']}, b_sum={ex['b_sum']}, "
                  f"F=G-b={ex['F']}, #paths={ex['count_paths']}")
    else:
        print("G = 2*mu_arc HOLDS IN ALL CASES!")
        print()
        print("THEOREM CANDIDATE:")
        print("  G(T,v,i,j) = sum_{P': i->j in T-v} (insact(v,P') - 1)")
        print("             = 2 * sum_{C through v using arc i->j} mu_T(C)")
        print()
        print("This means:")
        print("  F = G - b_sum = 2*mu_arc_sum - b_sum")
        print("  adj_T(i,j) - adj_{T-v}(i,j) = 2*mu_arc_sum - b_sum")
    print()


if __name__ == "__main__":
    run_G_analysis(n=5)
    run_G_analysis(n=6, num_samples=30)
    run_G_analysis(n=7, num_samples=10)
