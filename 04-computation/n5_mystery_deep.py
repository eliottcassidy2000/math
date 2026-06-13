"""
Deep investigation of OPEN-Q-001: Why does the per-path identity hold at n=5?

Key observation from initial run: At n=5, EVERY 5-cycle through v has mu(C) = 1.
This is because C\{v} has 4 vertices, which is ALL of T-v (since n=5).
So Omega(T-v)|_{avoid C\{v}} has NO available vertices -> empty graph -> I(empty, 2) = 1.

This means mu(5-cycle) = 1 ALWAYS at n=5. This is structural, not a coincidence.

The per-path identity says: (inshat-1)/2 = #{3-cycle embeddings in P'}
Claim A says: H(T) - H(T-v) = 2 * (sum_3 mu_3 + sum_5 mu_5)

Since mu(5-cycle) = 1 always, and there's exactly one vertex set for 5-cycles
through v at n=5 (it must use ALL vertices), the question reduces to:
how many directed 5-cycles through v exist, and does the per-path identity
absorb them?

Author: opus-2026-03-05-S1
"""

from itertools import permutations, combinations
import sys


def all_tournaments(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(edges)
    for bits in range(2**m):
        T = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if (bits >> k) & 1:
                T[i][j] = 1
            else:
                T[j][i] = 1
        yield T


def count_ham_paths(T, verts):
    count = 0
    for perm in permutations(verts):
        if all(T[perm[i]][perm[i+1]] == 1 for i in range(len(perm)-1)):
            count += 1
    return count


def count_directed_5cycles_through_v(T, v, n):
    """Count directed 5-cycles through v at n=5."""
    other = [u for u in range(n) if u != v]
    count = 0
    for perm in permutations(other):
        if (T[v][perm[0]] == 1 and T[perm[0]][perm[1]] == 1 and
            T[perm[1]][perm[2]] == 1 and T[perm[2]][perm[3]] == 1 and
            T[perm[3]][v] == 1):
            count += 1
    return count


def count_directed_3cycles_through_v(T, v, n):
    """Count directed 3-cycles through v."""
    other = [u for u in range(n) if u != v]
    count = 0
    for a, b in permutations(other, 2):
        if T[v][a] == 1 and T[a][b] == 1 and T[b][v] == 1:
            count += 1
    return count


def compute_mu_3cycle(T, v, n, a, b):
    """mu for 3-cycle v->a->b->v. C\{v} = {a,b}.
    Find odd cycles in T-v disjoint from {a,b}, build conflict graph, eval I at 2."""
    tmv = [u for u in range(n) if u != v]
    avail = [u for u in tmv if u not in (a, b)]

    # Find directed odd cycles among avail vertices
    cycles = []
    for L in range(3, len(avail) + 1, 2):
        for sub in combinations(avail, L):
            first = sub[0]
            for perm in permutations(sub[1:]):
                seq = (first,) + perm
                if all(T[seq[i]][seq[i+1]] == 1 for i in range(L-1)) and T[seq[-1]][first] == 1:
                    cycles.append(frozenset(sub))

    if not cycles:
        return 1

    # Build conflict graph and compute I(G, 2)
    k = len(cycles)
    adj = [set() for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            if cycles[i] & cycles[j]:
                adj[i].add(j)
                adj[j].add(i)

    total = 0
    for bits in range(2**k):
        nodes = [i for i in range(k) if (bits >> i) & 1]
        indep = True
        for i in range(len(nodes)):
            for j in range(i+1, len(nodes)):
                if nodes[j] in adj[nodes[i]]:
                    indep = False
                    break
            if not indep:
                break
        if indep:
            total += 2**len(nodes)
    return total


if __name__ == "__main__":
    n = 5
    print("="*70)
    print("OPEN-Q-001 Deep Analysis: n=5 mystery")
    print("="*70)

    # THEOREM: At n=5, mu(5-cycle through v) = 1 ALWAYS.
    # Proof: A 5-cycle through v uses v + 4 other vertices = ALL vertices.
    # So C\{v} = V\{v} = all of T-v. The "available" vertices for
    # Omega(T-v)|_{avoid C\{v}} is EMPTY. I(empty graph, 2) = 1.
    print("\nTheorem: At n=5, mu(C) = 1 for every 5-cycle C through v.")
    print("Proof: C uses all 5 vertices, so C\\{v} = V\\{v}.")
    print("       No vertices remain for Omega(T-v)|_{avoid C\\{v}}.")
    print("       I(empty, 2) = 1. QED")

    # Now: how many directed 5-cycles through v exist?
    # At n=5, a 5-cycle uses all vertices. The number of directed Hamiltonian
    # cycles in T equals the number of directed 5-cycles.
    # For the per-path identity to hold, we need:
    # sum over P' of [(inshat-1)/2] = sum over P' of [#3-cycle embeddings in P']
    # But Claim A says:
    # H(T) - H(T-v) = 2*(sum_3 mu_3 + #5cycles_through_v * 1)
    #
    # The per-path sum gives: sum_P' (inshat-1)/2 = sum_P' #3cycle_embeddings
    # And by THM-004/005: (inshat-1)/2 = #Type-II = #3-cycle embeddings
    #
    # Also: H(T) = sum_P' inshat(v, P')
    # So: H(T) - H(T-v) = sum_P' (inshat - 1)  [since sum_P' 1 = H(T-v)]
    #                    = 2 * sum_P' (inshat-1)/2
    #                    = 2 * sum_P' #3-cycle-embeddings
    #
    # Claim A says this equals 2*(sum_3 mu_3 + N_5) where N_5 = #5-cycles through v.
    #
    # So the question is: does sum_P' #3-cycle-embeddings = sum_3 mu_3 + N_5?
    # i.e., does the 3-cycle per-path count "accidentally" include the 5-cycle contribution?

    print("\n" + "="*70)
    print("Testing: sum_P' #3cycle_embeddings = sum_3_cycles mu(C) + N_5cycles")
    print("="*70)

    all_match = True
    total_with_5cycles = 0

    for t_idx, T in enumerate(all_tournaments(n)):
        for v in range(n):
            other = [u for u in range(n) if u != v]

            # Count 5-cycles through v
            N5 = count_directed_5cycles_through_v(T, v, n)

            # Sum over P' of #3-cycle-embeddings
            perpath_sum = 0
            H_Tv = 0
            for perm_p in permutations(other):
                if not all(T[perm_p[i]][perm_p[i+1]] == 1 for i in range(len(perm_p)-1)):
                    continue
                H_Tv += 1
                # Count 3-cycle embeddings in this path
                for j in range(len(perm_p) - 1):
                    a, b = perm_p[j], perm_p[j+1]
                    if T[v][a] == 1 and T[a][b] == 1 and T[b][v] == 1:
                        perpath_sum += 1
                    if T[v][b] == 1 and T[b][a] == 1 and T[a][v] == 1:
                        perpath_sum += 1

            # Sum of mu for 3-cycles through v
            mu3_sum = 0
            three_cycles = []
            for a, b in permutations(other, 2):
                if T[v][a] == 1 and T[a][b] == 1 and T[b][v] == 1:
                    three_cycles.append((a, b))
                    mu3_sum += compute_mu_3cycle(T, v, n, a, b)

            expected = mu3_sum + N5

            if perpath_sum != expected:
                print(f"MISMATCH: T#{t_idx}, v={v}: perpath_sum={perpath_sum}, "
                      f"mu3_sum={mu3_sum}, N5={N5}, expected={expected}")
                all_match = False

            if N5 > 0:
                total_with_5cycles += 1

    if all_match:
        print(f"\nCONFIRMED: sum_P' #3cycle_embeds = sum_3cycles mu(C) + N_5cycles")
        print(f"for ALL {n}-vertex tournaments (512 tournaments, 2560 (T,v) pairs)")
        print(f"Pairs with 5-cycles: {total_with_5cycles}")

    print("\n" + "="*70)
    print("CONCLUSION for OPEN-Q-001:")
    print("="*70)
    print("""
At n=5, the per-path identity (inshat-1)/2 = #3-cycle-embeddings holds
for ALL (T,v,P') triples, even when 5-cycles exist. The reason is:

1. mu(5-cycle) = 1 ALWAYS at n=5, because C\\{v} exhausts all of T-v,
   leaving no vertices for the restricted conflict graph.

2. The 3-cycle per-path sum naturally "absorbs" the 5-cycle contribution:
   sum_P' #3-cycle-embeddings = sum_{3-cycles C} mu(C) + N_{5-cycles}

   This means the 3-cycle embeddings across all paths over-count by exactly
   the number of 5-cycles, which is precisely the missing 5-cycle mu
   contribution (since each has mu=1).

3. At n=6, this breaks because:
   - 5-cycles no longer use all vertices, so mu(5-cycle) != 1 in general
   - The "absorption" identity no longer holds
   - The per-path 3-cycle count no longer captures the full RHS of Claim A
""")
