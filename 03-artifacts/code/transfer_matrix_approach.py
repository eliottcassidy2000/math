#!/usr/bin/env python3
"""
Explore transfer matrix / permanent representations of H(T) and I(Omega(T), 2).

Key question: is there a matrix M(T) such that H(T) = per(M) or H(T) = tr(M^k)?

Known: H(T) = sum over Ham paths = permanent of a certain (n-1)! x (n-1)! matrix,
but we want something simpler.

Also: I(Omega(T), 2) = sum_{k>=0} alpha_k * 2^k where alpha_k = #independent sets of size k.

The independence polynomial is related to the matching polynomial of the complement graph.
For the conflict graph Omega(T), we have I(Omega, x) = sum_k alpha_k x^k.

Key identity to test: H(T) = per(J + A) where A is the adjacency matrix?
Or: H(T) = sum of permanents of submatrices?

Instance: opus-2026-03-05-S3
"""

import sys
sys.path.insert(0, '.')
from tournament_lib import (
    hamiltonian_path_count, find_odd_cycles, independence_poly_at_fast,
    all_tournaments
)
from itertools import permutations
import numpy as np


def permanent(M):
    """Compute permanent of a square matrix (Ryser's formula)."""
    n = len(M)
    if n == 0:
        return 1
    # Ryser's formula: per(M) = (-1)^n sum_{S subset {1..n}} (-1)^|S| prod_i sum_{j in S} M[i][j]
    total = 0
    for mask in range(1, 1 << n):
        # S = set of columns in mask
        bits = bin(mask).count('1')
        prod = 1
        for i in range(n):
            row_sum = 0
            for j in range(n):
                if (mask >> j) & 1:
                    row_sum += M[i][j]
            prod *= row_sum
        total += ((-1) ** (n - bits)) * prod
    return ((-1) ** n) * total


def test_permanent_formulas(T):
    """Test various permanent-based formulas for H(T)."""
    n = len(T)
    H = hamiltonian_path_count(T)

    # Test 1: per(A) where A is adjacency matrix
    A = [[T[i][j] for j in range(n)] for i in range(n)]
    perA = permanent(A)

    # Test 2: per(A) with diagonal zeroed
    A_nodiag = [[T[i][j] if i != j else 0 for j in range(n)] for i in range(n)]
    perA_nodiag = permanent(A_nodiag)

    # Test 3: H(T) should equal the number of Ham paths = sum over permutations
    # of prod T[pi(k)][pi(k+1)]. This is NOT the permanent of A.
    # permanent(A) = sum_pi prod_i A[i][pi(i)] = sum_pi prod T[i][pi(i)]
    # This counts permutation matrices, not paths.

    # H(T) = sum over linear orderings pi: prod_{k=0}^{n-2} T[pi(k)][pi(k+1)]
    # This can be written as: sum_pi prod_{k} T[pi(k)][pi(k+1)]
    # = tr(A^{n-1}) in some sense? No, that overcounts.

    # Actually H(T) = sum_pi prod A[pi(k)][pi(k+1)] where sum is over all n! permutations
    # This is related to the permanent of a DIFFERENT matrix.

    # Let B[i][j] = T[i][j] for i != j, B[i][i] = 0.
    # Then H(T) = sum_{v_1,...,v_n distinct} prod B[v_k][v_{k+1}]
    # This equals the (1,1,...,1) coefficient of prod (x_{v_k} B[v_k][v_{k+1}])
    # which is related to the hafnian or the permanent of the path matrix.

    # Actually, the standard result: H(T) can be computed via the inclusion-exclusion
    # formula for the number of Hamiltonian paths, using the transfer matrix method:
    # H(T) = sum_{S subset V, |S|>=2} (-1)^{|S|+...} ... (complicated)

    # Or using the DP: H(T) = sum_v dp[V][v] where dp[S][v] = #{paths through S ending at v}

    # A simpler observation: H(T) relates to the TRACE of the path algebra.
    # Define M_k = sum over k-vertex paths. Then H(T) = sum_v (M_{n-1} * e_v).

    return {
        'H': H,
        'per(A)': perA,
        'per(A_nodiag)': perA_nodiag,
    }


def explore_skew_permanent(T):
    """
    Test: H(T) = per(A - A^T + I) / something?

    For a tournament, A - A^T = 2A - J + I (skew-symmetric part).
    The pfaffian of the skew-adjacency matrix is related to perfect matchings.
    """
    n = len(T)
    H = hamiltonian_path_count(T)

    # Skew adjacency: S[i][j] = T[i][j] - T[j][i] for i < j
    S = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            S[i][j] = T[i][j] - T[j][i]  # +1 if i->j, -1 if j->i, 0 on diagonal

    return {'H': H, 'skew': S}


def test_ocf_structure(n_max=6):
    """
    For small tournaments, compute H(T) and I(Omega, 2) and look for
    a common matrix expression.
    """
    print(f"Testing matrix representations for n <= {n_max}")
    print("="*60)

    for n in range(3, n_max + 1):
        print(f"\nn = {n}:")
        count = 0
        for T in all_tournaments(n):
            count += 1
            if count > 5:
                break
            result = test_permanent_formulas(T)
            H = result['H']
            cycles = find_odd_cycles(T)
            I = independence_poly_at_fast(cycles, 2)
            print(f"  H={H}, I={I}, per(A)={result['per(A)']}, "
                  f"per(A_nodiag)={result['per(A_nodiag)']}, "
                  f"H==I: {H==I}")


def test_hadamard_product_idea():
    """
    New idea: H(T) = sum over orderings pi of prod T[pi(k)][pi(k+1)].
    Let's write this as a bilinear form.

    Define for each pair (u,v), the "transition" T[u][v].
    A Ham path is a sequence of n-1 transitions.
    H(T) = sum_pi prod_{k=1}^{n-1} T[pi(k-1)][pi(k)]

    This is the permanent of the "path matrix" P where P[k][v] encodes
    "vertex v at position k". But this matrix is not well-defined without
    fixing the start.

    Alternative: H(T) = sum_{v} h_end(T, V, v) where h_end is the
    #paths through all vertices ending at v. This is the DP.

    Key insight: the DP satisfies dp[S][v] = sum_{u in S, u!=v, T[u][v]=1} dp[S\{v}][u].
    In matrix form: dp[S] = A_S * dp[S\{v}] for each v.

    This is a product of transfer matrices! Specifically:
    dp[V] = product over vertices added in some order.
    But the order matters (it's not commutative).
    """
    pass


def test_chromatic_independence_relation():
    """
    I(G, x) at x=2 counts 2-colorings of independent sets.
    Equivalently: I(G, 2) = sum over vertex subsets U: [U is independent in G] * 2^|U|
    = sum_U [U independent] * 2^|U|
    = number of pairs (U, f) where U is independent and f: U -> {red, blue}

    So I(Omega(T), 2) counts "colored independent cycle sets" in Omega(T).

    H(T) should equal this count. This means each Ham path of T
    corresponds to a colored independent cycle set in Omega(T).

    Can we find this bijection?
    """
    for n in [5, 6]:
        print(f"\nn = {n}:")
        count = 0
        for T in all_tournaments(n):
            count += 1
            if count > 3:
                break
            H = hamiltonian_path_count(T)
            cycles = find_odd_cycles(T)
            I = independence_poly_at_fast(cycles, 2)

            # Decompose I by size
            # alpha_0 = 1 (empty set), alpha_1 = |cycles|, etc.
            num_cycles = len(cycles)
            # Count VD pairs
            vd_pairs = 0
            for ci in range(num_cycles):
                for cj in range(ci+1, num_cycles):
                    vi = set(cycles[ci])
                    vj = set(cycles[cj])
                    if vi.isdisjoint(vj):
                        vd_pairs += 1

            I_check = 1 + 2*num_cycles + 4*vd_pairs
            print(f"  H={H}, #cycles={num_cycles}, #VD_pairs={vd_pairs}, "
                  f"I=1+2*{num_cycles}+4*{vd_pairs}={I_check}, match={H==I_check}")

            # So H = 1 + 2*(#odd_cycles) + 4*(#VD_pairs) + ...
            # The "1" corresponds to... what? The transitive Ham path?
            # The "2*#cycles" corresponds to... 2 Ham paths per odd cycle?
            # The "4*#VD_pairs" corresponds to... 4 paths per VD pair?


def main():
    test_ocf_structure(6)
    print("\n" + "="*60)
    test_chromatic_independence_relation()


if __name__ == "__main__":
    main()
