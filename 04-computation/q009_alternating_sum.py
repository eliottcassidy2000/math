"""
Q-009: Prove the alternating sum identity algebraically.

The key lemma: sum_S (-1)^|S| Delta(S, others\S) = 0
where Delta(S,R) = L_j(S)*R_i(R) - L_i(S)*R_j(R).

Approach 1: Expand as a polynomial and check coefficient-by-coefficient.
Approach 2: Find a combinatorial proof (sign-reversing involution).
Approach 3: Matrix identity.

Let me first understand the structure better by examining small cases.

Author: opus-2026-03-05-S4
"""

import random
from itertools import permutations, combinations


def random_tournament(n):
    T = [[0]*n for _ in range(n)]
    for a in range(n):
        for b in range(a+1, n):
            if random.random() < 0.5:
                T[a][b] = 1
            else:
                T[b][a] = 1
    return T


def h_end_dp(T, verts, target):
    n = len(verts)
    if n == 1: return 1
    idx = {v: k for k, v in enumerate(verts)}
    ti = idx[target]
    dp = [[0]*n for _ in range(1<<n)]
    for k in range(n): dp[1<<k][k] = 1
    full = (1<<n)-1
    for mask in range(1, 1<<n):
        for v in range(n):
            if not (mask & (1<<v)) or dp[mask][v]==0: continue
            for u in range(n):
                if mask & (1<<u): continue
                if T[verts[v]][verts[u]]: dp[mask|(1<<u)][u] += dp[mask][v]
    return dp[full][ti]


def h_start_dp(T, verts, source):
    n = len(verts)
    if n == 1: return 1
    idx = {v: k for k, v in enumerate(verts)}
    si = idx[source]
    dp = [[0]*n for _ in range(1<<n)]
    dp[1<<si][si] = 1
    full = (1<<n)-1
    for mask in range(1, 1<<n):
        for v in range(n):
            if not (mask & (1<<v)) or dp[mask][v]==0: continue
            for u in range(n):
                if mask & (1<<u): continue
                if T[verts[v]][verts[u]]: dp[mask|(1<<u)][u] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


def analyze_alternating_sum(T, i, j, others):
    """Decompose the alternating sum into per-vertex contributions."""
    m = len(others)

    # Compute the alternating sum in detail
    alt_sum = 0
    terms = {}

    for size in range(m + 1):
        for S in combinations(others, size):
            R = tuple(x for x in others if x not in S)
            S_list = sorted(S)
            R_list = sorted(R)

            if not S_list:
                L_j = 1; L_i = 1
            else:
                L_j = sum(h_end_dp(T, S_list, a) * T[a][j] for a in S_list)
                L_i = sum(h_end_dp(T, S_list, a) * T[a][i] for a in S_list)

            if not R_list:
                R_i = 1; R_j = 1
            else:
                R_i = sum(h_start_dp(T, R_list, b) * T[i][b] for b in R_list)
                R_j = sum(h_start_dp(T, R_list, b) * T[j][b] for b in R_list)

            delta = L_j * R_i - L_i * R_j
            sign = (-1) ** size
            terms[S] = (sign, delta, L_j, L_i, R_i, R_j)
            alt_sum += sign * delta

    return alt_sum, terms


def test_factored_form():
    """Test if the alternating sum factors nicely."""
    random.seed(42)

    for n in [5, 6, 7]:
        print(f"\n=== n={n} ===")
        for trial in range(5):
            T = random_tournament(n)
            i, j = random.sample(range(n), 2)
            if T[i][j] == 0: i, j = j, i
            others = [x for x in range(n) if x != i and x != j]

            alt_sum, terms = analyze_alternating_sum(T, i, j, others)

            # Compute per-vertex s values
            s = {x: 1 - T[x][i] - T[j][x] for x in others}

            # Compute Alt_L = sum_S (-1)^|S| L_j(S) and Alt_L' = sum_S (-1)^|S| L_i(S)
            # Similarly for R side
            Alt_Lj = sum(terms[S][0] * terms[S][2] for S in terms)
            Alt_Li = sum(terms[S][0] * terms[S][3] for S in terms)

            # Actually, since Delta = L_j*R_i - L_i*R_j, the alternating sum is:
            # sum (-1)^|S| [L_j(S)*R_i(R) - L_i(S)*R_j(R)]
            # This is NOT simply factored because R depends on S.

            # Test: does sum_S (-1)^|S| L_j(S) * R_i(others\S) = 0?
            cross1 = sum(terms[S][0] * terms[S][2] * terms[S][4] for S in terms)
            cross2 = sum(terms[S][0] * terms[S][3] * terms[S][5] for S in terms)

            print(f"  Trial {trial}: alt_sum={alt_sum}, "
                  f"cross1={cross1}, cross2={cross2}, "
                  f"cross1==cross2: {cross1==cross2}")


def test_vertex_swap_symmetry():
    """Test if the alternating sum has a sign-reversing involution.

    Key idea: for each subset S, consider the "complement swap" S -> others\S.
    This sends |S| -> m-|S|, and if m is even, parity flips.
    """
    random.seed(42)

    for n in [5, 6, 7]:
        print(f"\n=== n={n}, m={n-2} ===")
        for trial in range(5):
            T = random_tournament(n)
            i, j = random.sample(range(n), 2)
            if T[i][j] == 0: i, j = j, i
            others = [x for x in range(n) if x != i and x != j]
            m = len(others)

            _, terms = analyze_alternating_sum(T, i, j, others)

            # Check: Delta(S, R) vs Delta(R, S) ?
            for S in terms:
                R = tuple(x for x in others if x not in S)
                if S <= R:  # avoid double counting
                    d_SR = terms[S][1]
                    d_RS = terms[R][1] if R in terms else None
                    if d_RS is not None and trial == 0:
                        print(f"  S={S}, R={R}: Delta(S,R)={d_SR}, Delta(R,S)={d_RS}, "
                              f"sum={d_SR+d_RS}, neg_sum={d_SR-d_RS}")

            # Key test: is Delta(S,R) = -Delta(R,S) when |S|+|R| = m?
            # If so, (-1)^|S| Delta(S,R) + (-1)^|R| Delta(R,S)
            # = (-1)^|S| * [Delta(S,R) + (-1)^m * Delta(R,S)]
            # If m is even: = (-1)^|S| [Delta(S,R) + Delta(R,S)]
            # If m is odd: = (-1)^|S| [Delta(S,R) - Delta(R,S)]
            # For the alternating sum to be 0 with m even, we need Delta(S,R) = -Delta(R,S).
            # For m odd, we need Delta(S,R) = Delta(R,S).

            symmetry_holds = True
            for S in terms:
                R = tuple(x for x in others if x not in S)
                d_SR = terms[S][1]
                d_RS = terms[R][1]
                if m % 2 == 0:
                    if d_SR != -d_RS:
                        symmetry_holds = False
                        break
                else:
                    if d_SR != d_RS:
                        symmetry_holds = False
                        break

            if trial < 3 or not symmetry_holds:
                expected = "antisymmetric" if m % 2 == 0 else "symmetric"
                print(f"  Trial {trial}: m={m}, {expected}: {symmetry_holds}")


if __name__ == "__main__":
    print("=== Testing factored form ===")
    test_factored_form()

    print("\n\n=== Testing vertex swap symmetry ===")
    test_vertex_swap_symmetry()
