"""
Q-009: Decompose adj(i,j) to match the ΔI formula.

The identity to prove at n=6:
  adj(i,j) - adj'(j,i) = -2·Σ_x s_x·H(B_x) + 2·(D5-C5)

adj(i,j) = # Ham paths of T with i immediately before j.
Split by position a (number of vertices before i):

  adj(i,j) = Σ_{a=0}^{4} Σ_{L:|L|=a, R=others\L}
              h_end(L∪{i}, i) · h_start({j}∪R, j)

Similarly, adj'(j,i) = adj_{T'}(j,i) = Σ same with j before i.
Since T and T' agree on all arcs except i↔j:

  adj'(j,i) = Σ_{a=0}^{4} Σ_{L:|L|=a, R=others\L}
               h_end(L∪{j}, j) · h_start({i}∪R, i)

(using T arcs, since L∪{j} and {i}∪R never contain both i and j)

Can we match terms with the RHS by grouping by x (excluded vertex)?

Author: opus-2026-03-05-S2
"""

import sys
sys.path.insert(0, '03-artifacts/code')
from tournament_lib import all_tournaments, hamiltonian_path_count
from itertools import permutations, combinations


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def sub_H(T, verts):
    verts = sorted(verts)
    m = len(verts)
    if m <= 1:
        return 1
    count = 0
    for perm in permutations(verts):
        valid = True
        for k in range(m - 1):
            if T[perm[k]][perm[k + 1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count


def h_end(T, verts, target):
    """# Ham paths of T[verts] ending at target."""
    verts = sorted(verts)
    assert target in verts
    count = 0
    for perm in permutations(verts):
        if perm[-1] != target:
            continue
        valid = True
        for k in range(len(perm) - 1):
            if T[perm[k]][perm[k + 1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count


def h_start(T, verts, source):
    """# Ham paths of T[verts] starting at source."""
    verts = sorted(verts)
    assert source in verts
    count = 0
    for perm in permutations(verts):
        if perm[0] != source:
            continue
        valid = True
        for k in range(len(perm) - 1):
            if T[perm[k]][perm[k + 1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count


if __name__ == "__main__":
    n = 6
    print(f"Q-009: Decompose adj(i,j) at n={n}")
    print("=" * 70)

    # For each (T, i->j), decompose adj(i,j) by split position a
    # and check relationship to Σ_x s_x·H(B_x)

    # Key insight: group terms by which vertex x is NOT on the same
    # side as both i and j.
    # Actually, let me try a different grouping: by excluded vertex x.
    # Split V\{i,j} = {a,b,c,d} into L (before i) and R (after j).
    # For each split (L,R), there's no natural "x" unless we can
    # relate the sum to B_x terms.

    # Let me compute adj(i,j) decomposed by a (= |L|) and also by
    # which endpoint of B_x gets connected to j (before) and i (after).

    # First, verify the basic decomposition
    total = 0
    match = 0

    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= 100:
            break
        for i in range(n):
            for j in range(i + 1, n):
                if T[i][j] == 0:
                    continue
                total += 1
                others = [x for x in range(n) if x != i and x != j]

                # Compute adj(i,j) by decomposition
                adj_decomp = 0
                for a in range(5):  # |L| = 0,...,4
                    for L in combinations(others, a):
                        L_set = set(L)
                        R = [x for x in others if x not in L_set]
                        he = h_end(T, list(L) + [i], i)
                        hs = h_start(T, [j] + R, j)
                        adj_decomp += he * hs

                # Direct computation
                adj_direct = 0
                for perm in permutations(range(n)):
                    valid = True
                    for k in range(n - 1):
                        if T[perm[k]][perm[k + 1]] != 1:
                            valid = False
                            break
                    if valid:
                        for k in range(n - 1):
                            if perm[k] == i and perm[k + 1] == j:
                                adj_direct += 1

                if adj_decomp == adj_direct:
                    match += 1

    print(f"Decomposition check: {match}/{total}")

    # Now, the key question: what is adj(i,j) - adj'(j,i)?
    # Let's express it differently.
    # adj(i,j) - adj'(j,i) = Σ_{L,R} [h_end(L∪{i},i)·h_start({j}∪R,j)
    #                                  - h_end(L∪{j},j)·h_start({i}∪R,i)]
    #
    # For a=0: R = all others = {a,b,c,d}
    #   1 · h_start({j,a,b,c,d}, j) - 1 · h_start({i,a,b,c,d}, i)
    #   = H_5_from_j - H_5_from_i (paths starting from j vs i in 5-vertex tournaments)
    #
    # For a=4: L = all others
    #   h_end({a,b,c,d,i}, i) · 1 - h_end({a,b,c,d,j}, j) · 1
    #   = H_5_to_i - H_5_to_j
    #
    # For a=1: L={x}, R=others\{x}
    #   Σ_x T[x][i] · h_start({j}∪R_x, j) - T[x][j] · h_start({i}∪R_x, i)
    #   where R_x = others\{x}
    #
    # Note: {j}∪R_x has 4 vertices ({j} + 3 others). This is {j}∪B_x where B_x = others\{x}.
    # Similarly: {i}∪R_x = {i}∪B_x.
    #
    # For a=3: L has 3 vertices = B_x for each x not in L
    #   Σ over 3-element subsets L of others:
    #     h_end(L∪{i}, i) · T[j][r] - h_end(L∪{j}, j) · T[i][r]
    #   where {r} = others\L

    # This is getting complex. Let me instead look at the SYMMETRY between
    # adj(i,j) and adj'(j,i).
    #
    # adj(i,j) = Σ paths with ...→i→j→...
    # adj'(j,i) = Σ paths with ...→j→i→...
    #
    # The key difference: in adj(i,j), i beats all vertices before it, j is
    # beaten by all vertices after it, and i→j. In adj'(j,i), j beats all
    # before it, i is beaten by all after, and j→i.

    # Let me try a completely different approach: express everything as a sum
    # over the 4 "other" vertices grouped by pairs.
    #
    # There are C(4,2) = 6 ways to split the 4 others into (L,R) with |L|=|R|=2.
    # Plus degenerate cases.

    # Actually, let me try the simplest possible test: for a=1 and a=3,
    # the terms involve B_x directly. Let me see if the a=1 and a=3 terms
    # alone match the -2·Σ s_x·H(B_x) term.

    print(f"\n--- Decompose adj(i,j) - adj'(j,i) by position a ---")
    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= 5:
            break
        for i in range(n):
            for j in range(i + 1, n):
                if T[i][j] == 0:
                    continue
                T2 = flip_arc(T, i, j)
                others = [x for x in range(n) if x != i and x != j]

                terms = {}
                for a in range(5):
                    term = 0
                    for L in combinations(others, a):
                        L_set = set(L)
                        R = [x for x in others if x not in L_set]
                        he_i = h_end(T, list(L) + [i], i)
                        hs_j = h_start(T, [j] + R, j)
                        he_j = h_end(T, list(L) + [j], j)
                        hs_i = h_start(T, [i] + R, i)
                        term += he_i * hs_j - he_j * hs_i
                    terms[a] = term

                total_diff = sum(terms.values())
                sx_hbx = sum((1 - T[x][i] - T[j][x]) * sub_H(T, [u for u in others if u != x])
                             for x in others)

                # 5-cycle term
                d5c5 = 0
                for x in others:
                    bx = [u for u in others if u != x]
                    for perm in permutations(bx):
                        internal = T[perm[0]][perm[1]] * T[perm[1]][perm[2]]
                        if internal:
                            d5c5 += (T[j][perm[0]] * T[perm[2]][i] -
                                     T[i][perm[0]] * T[perm[2]][j])

                rhs = -2 * sx_hbx + 2 * d5c5

                print(f"\nT#{t_idx} flip {i}->{j}:")
                print(f"  By position: {terms}")
                print(f"  Total: adj-adj' = {total_diff}")
                print(f"  RHS = -2*{sx_hbx} + 2*{d5c5} = {rhs}")
                print(f"  Match: {total_diff == rhs}")

                # Check: a=1 term vs s_x relationship
                a1_per_x = {}
                for x in others:
                    R_x = [u for u in others if u != x]
                    t1 = T[x][i] * h_start(T, [j] + R_x, j) - T[x][j] * h_start(T, [i] + R_x, i)
                    a1_per_x[x] = t1

                a3_per_x = {}
                for x in others:
                    L_x = [u for u in others if u != x]  # = B_x
                    t3 = h_end(T, L_x + [i], i) * T[j][x] - h_end(T, L_x + [j], j) * T[i][x]
                    a3_per_x[x] = t3

                print(f"  a=1 per x: {a1_per_x}")
                print(f"  a=3 per x: {a3_per_x}")
                print(f"  s_x: {dict((x, 1-T[x][i]-T[j][x]) for x in others)}")
                print(f"  H(B_x): {dict((x, sub_H(T, [u for u in others if u!=x])) for x in others)}")

                break
            else:
                continue
            break
