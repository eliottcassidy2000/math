"""
Even-Odd Split Lemma: deeper analysis.

The alternating sum sum_{S subset V\{i,j}} (-1)^|S| Delta(S, R) = 0
where Delta(S,R) = h_end(S+{i},i)*h_start({j}+R,j) - h_end(S+{j},j)*h_start({i}+R,i).

Let f(S) = h_end(S+{i},i) * h_start({j}+R,j)  [R = complement]
Let g(S) = h_end(S+{j},j) * h_start({i}+R,i)

The lemma says: sum (-1)^|S| [f(S) - g(S)] = 0
i.e., sum (-1)^|S| f(S) = sum (-1)^|S| g(S).

Question: does each side vanish individually?
Or is it only the difference that vanishes?

Also: what's the connection to the subset convolution / Mobius function perspective?
The alternating sum sum (-1)^|S| F(S) over S ⊆ [m] is the Mobius function evaluation
of F at the empty set, if F is viewed as a function on the Boolean lattice.

Instance: opus-2026-03-05-S4
"""

import sys
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import all_tournaments
from itertools import permutations, combinations


def h_end(T, verts, target):
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


def analyze_split(n, max_tournaments=20):
    print(f"=== Even-Odd Split Analysis at n={n} ===")

    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= max_tournaments:
            break
        for i in range(n):
            for j in range(n):
                if i == j or T[i][j] != 1:
                    continue
                others = [x for x in range(n) if x != i and x != j]
                m = len(others)

                alt_f = 0  # sum (-1)^|S| f(S)
                alt_g = 0  # sum (-1)^|S| g(S)

                for size in range(m + 1):
                    for S in combinations(others, size):
                        S_set = set(S)
                        R = [x for x in others if x not in S_set]
                        S_list = list(S)
                        sign = (-1) ** size

                        f_S = h_end(T, S_list + [i], i) * h_start(T, [j] + R, j)
                        g_S = h_end(T, S_list + [j], j) * h_start(T, [i] + R, i)

                        alt_f += sign * f_S
                        alt_g += sign * g_S

                if alt_f != alt_g or (t_idx < 3 and j == (i+1) % n):
                    print(f"  T#{t_idx} flip {i}->{j}: alt_f={alt_f}, alt_g={alt_g}, diff={alt_f-alt_g}")

                if alt_f != alt_g:
                    print(f"  *** INDIVIDUAL SIDES DIFFER! ***")
                    return False

    print(f"  Do f and g sides individually vanish?")
    # Check if alt_f = 0 individually
    any_nonzero_f = False
    any_nonzero_g = False
    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= max_tournaments:
            break
        for i in range(n):
            for j in range(n):
                if i == j or T[i][j] != 1:
                    continue
                others = [x for x in range(n) if x != i and x != j]
                m = len(others)

                alt_f = 0
                alt_g = 0
                for size in range(m + 1):
                    for S in combinations(others, size):
                        S_set = set(S)
                        R = [x for x in others if x not in S_set]
                        S_list = list(S)
                        sign = (-1) ** size
                        alt_f += sign * h_end(T, S_list + [i], i) * h_start(T, [j] + R, j)
                        alt_g += sign * h_end(T, S_list + [j], j) * h_start(T, [i] + R, i)

                if alt_f != 0:
                    any_nonzero_f = True
                if alt_g != 0:
                    any_nonzero_g = True

    print(f"    alt_f ever nonzero: {any_nonzero_f}")
    print(f"    alt_g ever nonzero: {any_nonzero_g}")

    if any_nonzero_f:
        print("  => Only the DIFFERENCE vanishes, not each side individually.")
    else:
        print("  => Both sides vanish individually!")

    return True


def factored_form_analysis(n, max_tournaments=10):
    """
    Try to understand WHY the alternating sum vanishes.

    Key insight: h_end(S+{i}, i) depends on arcs within S and from S to i.
    h_start({j}+R, j) depends on arcs within R and from j to R.
    These are independent sets of arcs (no overlap).

    So f(S) = h_end(S+{i},i) * h_start({j}+R,j) is a product of functions
    that depend on disjoint parts of the tournament (except for arcs between S and R).

    Wait -- arcs between S and R are NOT used by either factor!
    h_end only uses arcs within S+{i}, h_start only uses arcs within {j}+R.

    So f(S) depends on: arcs within S, arcs S->i (and i->S), arcs within R, arcs j->R (and R->j).
    NOT on arcs between S and R.
    """
    print(f"\n=== Factored Form Analysis at n={n} ===")
    print("Testing whether f(S) depends on arcs between S and R...")

    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= 2:
            break
        i, j = 0, 1
        if T[i][j] != 1:
            continue
        others = [x for x in range(n) if x != i and x != j]
        m = len(others)

        # For a fixed S, flip an arc between S and R and see if f(S) changes
        if m >= 2:
            S = [others[0]]
            R = others[1:]
            S_list = list(S)
            R_list = list(R)

            # Original f(S)
            f_orig = h_end(T, S_list + [i], i) * h_start(T, [j] + R_list, j)

            # Flip an arc between S[0] and R[0]
            a, b = S[0], R[0]
            T2 = [row[:] for row in T]
            T2[a][b], T2[b][a] = T2[b][a], T2[a][b]
            f_flip = h_end(T2, S_list + [i], i) * h_start(T2, [j] + R_list, j)

            print(f"  T#{t_idx}, S={S}, R={R}: f_orig={f_orig}, f_flip={f_flip}, changed={f_orig!=f_flip}")

    print("\nKey observation: f(S) = h_end(S+{i},i) * h_start({j}+R,j)")
    print("h_end(S+{i},i) uses arcs among S+{i} ONLY")
    print("h_start({j}+R,j) uses arcs among {j}+R ONLY")
    print("=> f(S) does NOT depend on arcs between S and R!")
    print("=> f is a 'product' function on the Boolean lattice")


def mobius_interpretation(n, max_tournaments=10):
    """
    If f(S) = A(S) * B(complement(S)) where A depends only on arcs within S+{i}
    and B depends only on arcs within {j}+R, then:

    sum (-1)^|S| f(S) = sum (-1)^|S| A(S) * B(complement(S))

    This is the "boolean convolution" of A and B evaluated via Mobius inversion.

    Can we show this vanishes by analyzing A and B separately?

    Actually: sum (-1)^|S| A(S)*B(V\S) is the Hadamard product of A*(-1)^|.| with B_rev.

    More concretely, define:
    a(S) = h_end(S+{i}, i)
    b(R) = h_start({j}+R, j)

    sum (-1)^|S| a(S)*b(V\S) = ?

    For this to vanish for ALL tournaments, it must be a structural identity.
    """
    print(f"\n=== Mobius / Boolean Convolution Analysis at n={n} ===")

    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= 5:
            break
        i, j = 0, 1
        if T[i][j] != 1:
            continue
        others = [x for x in range(n) if x != i and x != j]
        m = len(others)

        # Compute a(S) and b(S) tables
        a_table = {}  # a(S) = h_end(S+{i}, i)
        b_table = {}  # b(R) = h_start({j}+R, j)
        c_table = {}  # c(S) = h_end(S+{j}, j)
        d_table = {}  # d(R) = h_start({i}+R, i)

        for size in range(m + 1):
            for S in combinations(others, size):
                S_set = frozenset(S)
                R = tuple(x for x in others if x not in S_set)
                R_set = frozenset(R)

                a_table[S_set] = h_end(T, list(S) + [i], i)
                b_table[R_set] = h_start(T, [j] + list(R), j)
                c_table[S_set] = h_end(T, list(S) + [j], j)
                d_table[R_set] = h_start(T, [i] + list(R), i)

        # Check: sum (-1)^|S| a(S)*b(comp(S)) vs sum (-1)^|S| c(S)*d(comp(S))
        alt_ab = 0
        alt_cd = 0
        for size in range(m + 1):
            for S in combinations(others, size):
                S_set = frozenset(S)
                R_set = frozenset(others) - S_set
                sign = (-1) ** size
                alt_ab += sign * a_table[S_set] * b_table[R_set]
                alt_cd += sign * c_table[S_set] * d_table[R_set]

        print(f"  T#{t_idx}: alt(a*b) = {alt_ab}, alt(c*d) = {alt_cd}")

        # What if we look at a and c Mobius transforms separately?
        # a_hat(S) = sum_{T subset S} (-1)^{|S|-|T|} a(T)  [Mobius transform]
        a_hat = {}
        c_hat = {}
        for size in range(m + 1):
            for S in combinations(others, size):
                S_set = frozenset(S)
                val_a = 0
                val_c = 0
                for size2 in range(size + 1):
                    for T_sub in combinations(list(S_set), size2):
                        T_set = frozenset(T_sub)
                        sign = (-1) ** (size - size2)
                        val_a += sign * a_table[T_set]
                        val_c += sign * c_table[T_set]
                a_hat[S_set] = val_a
                c_hat[S_set] = val_c

        # Similarly for b and d (but complement direction)
        print(f"    a_hat values: {dict((tuple(sorted(k)), v) for k, v in a_hat.items() if v != 0)}")
        print(f"    c_hat values: {dict((tuple(sorted(k)), v) for k, v in c_hat.items() if v != 0)}")


if __name__ == "__main__":
    for n in [4, 5]:
        analyze_split(n, max_tournaments=10000)

    factored_form_analysis(5)
    mobius_interpretation(5, max_tournaments=5)
