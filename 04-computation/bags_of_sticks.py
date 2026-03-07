#!/usr/bin/env python3
"""
Explore the bags-of-sticks decomposition of U_T for tournaments.

A "bag of sticks" is a disjoint union of directed paths.
U_{P_lambda} is the Redei-Berge function of a bag of sticks with parts lambda.

For a tournament T, U_T can be expressed as a linear combination of U_{P_lambda}.

Key question: what does this decomposition look like in terms of our
independence polynomial / OCF framework?

Recall: ps_1(U_T) = H(T) (principal specialization = number of listings/paths).

For a bag of sticks P_lambda with lambda = (lambda_1, ..., lambda_k):
  ps_1(U_{P_lambda})(m) = number of listings of [n] such that within each stick,
    elements appear in the correct order.

Actually, U_{P_lambda} for a directed path P_n on [n] with edges 1->2->...->n:
  U_{P_n} = sum_{sigma in S_n} F_{P_n-Des(sigma)}

  where P_n-Des(sigma) = {i : (sigma_i, sigma_{i+1}) is an edge of P_n}
       = {i : sigma_{i+1} = sigma_i + 1}

For the COMPLETE directed path (edges i->i+1 for all i), this counts
permutations by their "ascent-to-successor" set.

Let me compute ps_1 for small bags of sticks and check formulas.

opus-2026-03-07-S36
"""
from itertools import permutations
from collections import defaultdict
from math import factorial

def compute_UT_ps1(n, edges, m):
    """Compute ps_1(U_X)(m) for digraph X on [n] with given edges.

    ps_1(U_X)(m) = #{(sigma, f): sigma is a listing of [n],
                     f: [n] -> [m] weakly increasing along sigma,
                     strictly increasing at X-descents}
    """
    vertices = list(range(n))
    edge_set = set(edges)

    count = 0
    for sigma in permutations(vertices):
        # X-descent set
        descents = set()
        for i in range(n-1):
            if (sigma[i], sigma[i+1]) in edge_set:
                descents.add(i)

        # Count number of weakly increasing f: [n] -> [m]
        # with f(sigma_i) < f(sigma_{i+1}) at descents
        # and f(sigma_i) <= f(sigma_{i+1}) at non-descents
        # This equals C(m + n - 1 - |descents|, n) ... actually it's
        # the number of weakly increasing sequences of length n from [m]
        # with strict increases at positions in descents.

        # By stars and bars:
        # Free positions = n - 1 - |descents| (non-descent boundaries)
        # Strict positions = |descents|
        # Total = C(m - |descents| + n - 1 - |descents|, n - 1)
        # Wait, this isn't right for general m.

        # Actually ps_1(F_I)(m) = C(m, n) for I = empty
        # and more generally C(m - |I|, n - |I|) ... no.

        # The correct formula: ps_1(F_I)(m) = number of sequences
        # 1 <= a_1 <= a_2 <= ... <= a_n <= m with a_j < a_{j+1} for j in I
        # = C(m, n) * (product over blocks) ...
        # Actually it's C(m + n - 1 - |I|, n) wait no...
        pass

    # Just compute directly for small m
    return None

def ham_paths_of_digraph(n, edges):
    """Count Hamiltonian paths in digraph (= listings with all consecutive
    pairs being edges). This is ps_1(U_X)(1) when X is a tournament."""
    edge_set = set(edges)
    count = 0
    for sigma in permutations(range(n)):
        if all((sigma[i], sigma[i+1]) in edge_set for i in range(n-1)):
            count += 1
    return count

def redei_berge_power_sum(n, edges):
    """Compute [p_lambda] U_X using Theorem 2.2 from Mitrovic-Stojadinovic.

    U_X = sum_{pi in S_V(X, X_bar)} (-1)^{phi(pi)} p_{type(pi)}

    where phi(pi) = sum_{gamma cycle of pi that is a cycle of X} (len(gamma) - 1)
    """
    # X_bar for tournament = opposite tournament (complement edges)
    # For a tournament, edge_set union opp_edge_set = all pairs
    edge_set = set(edges)
    opp_edges = set((j, i) for (i, j) in edge_set)  # T^op = reverse all arcs

    coeffs = defaultdict(int)  # lambda -> coefficient

    for sigma in permutations(range(n)):
        # Decompose sigma into cycles
        visited = [False] * n
        cycles = []
        for start in range(n):
            if visited[start]:
                continue
            cycle = []
            curr = start
            while not visited[curr]:
                visited[curr] = True
                cycle.append(curr)
                curr = sigma[curr]
            if len(cycle) > 0:
                cycles.append(tuple(cycle))

        # Check: each nontrivial cycle must be a directed cycle of X or X_bar
        valid = True
        phi = 0
        for cyc in cycles:
            if len(cyc) == 1:
                continue  # trivial cycle, always OK
            # Check if cyc is a directed cycle of X
            is_X_cycle = all((cyc[i], cyc[(i+1) % len(cyc)]) in edge_set for i in range(len(cyc)))
            is_Xbar_cycle = all((cyc[i], cyc[(i+1) % len(cyc)]) in opp_edges for i in range(len(cyc)))

            if not is_X_cycle and not is_Xbar_cycle:
                valid = False
                break

            if is_X_cycle:
                phi += len(cyc) - 1

        if valid:
            cycle_type = tuple(sorted([len(c) for c in cycles], reverse=True))
            coeffs[cycle_type] += (-1)**phi

    return dict(coeffs)

def main():
    print("=== Rédei-Berge Power Sum Coefficients ===\n")

    # Transitive tournament on 5 vertices: 0->1->2->3->4 (plus all i->j for i<j)
    n = 5
    trans_edges = [(i, j) for i in range(n) for j in range(i+1, n)]

    print(f"Transitive tournament T_5 (edges i->j for i<j):")
    H = ham_paths_of_digraph(n, trans_edges)
    print(f"  H(T) = {H}")

    coeffs = redei_berge_power_sum(n, trans_edges)
    print(f"  [p_lambda] U_T:")
    for lam in sorted(coeffs.keys()):
        if coeffs[lam] != 0:
            print(f"    p_{lam}: {coeffs[lam]}")

    # Check: specializing p_k -> 1 for odd k, p_k -> 0 for even k,
    # should give... wait, for tournaments U_T specializes differently.
    # ps_1(U_T)(1) counts the number of listings = n! / 1 = ... no.
    # Actually for a tournament, ps_1(U_T)(1) = H(T) by Corollary 20.

    # The OCF specialization: p_k -> 2 for odd k >= 3, p_1 -> 1, others 0?
    # From GS-OCF bridge: U_T = sum_{S indep} 2^|S| prod p_{len(c)} p_1^{n-sum_len}
    # Specializing p_k -> 1 for all k: get sum_{S indep} 2^|S| = I(Omega, 2) = H(T)

    print(f"\n  Specialization check (p_k -> 1 for all k):")
    H_from_coeffs = sum(coeffs[lam] for lam in coeffs)
    print(f"    sum of coeffs = {H_from_coeffs} (should be {H})")

    # Now a cyclic tournament
    print(f"\nCyclic tournament C_5 (edges i->(i+1 mod 5), i->(i+2 mod 5)):")
    cyc_edges = [(i, (i+1) % n) for i in range(n)] + [(i, (i+2) % n) for i in range(n)]
    H_cyc = ham_paths_of_digraph(n, cyc_edges)
    print(f"  H(T) = {H_cyc}")

    coeffs_cyc = redei_berge_power_sum(n, cyc_edges)
    print(f"  [p_lambda] U_T:")
    for lam in sorted(coeffs_cyc.keys()):
        if coeffs_cyc[lam] != 0:
            print(f"    p_{lam}: {coeffs_cyc[lam]}")

    H_from_coeffs_cyc = sum(coeffs_cyc[lam] for lam in coeffs_cyc)
    print(f"  sum of coeffs = {H_from_coeffs_cyc} (should be {H_cyc})")

    # Paley T_7
    print(f"\nPaley tournament T_7:")
    n = 7
    QR = {1, 2, 4}  # QR mod 7
    paley_edges = []
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in QR:
                paley_edges.append((i, j))
    H_p7 = ham_paths_of_digraph(n, paley_edges)
    print(f"  H(T_7) = {H_p7}")

    coeffs_p7 = redei_berge_power_sum(n, paley_edges)
    print(f"  [p_lambda] U_T_7:")
    for lam in sorted(coeffs_p7.keys()):
        if coeffs_p7[lam] != 0:
            print(f"    p_{lam}: {coeffs_p7[lam]}")

    # Check partition into odd vs even parts
    odd_only_sum = sum(coeffs_p7[lam] for lam in coeffs_p7
                       if all(l % 2 == 1 for l in lam))
    all_sum = sum(coeffs_p7[lam] for lam in coeffs_p7)
    print(f"  sum (all): {all_sum}")
    print(f"  sum (odd parts only): {odd_only_sum}")
    print(f"  sum (even parts): {all_sum - odd_only_sum}")
    print(f"  H(T_7) = {H_p7} (should match sum(all))")

if __name__ == "__main__":
    main()
