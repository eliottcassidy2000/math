#!/usr/bin/env python3
"""
Grinberg-Stanley U_T computation for small tournaments.

Computes:
1. E_T(t) = sum_H t^{des(H)} for Hamiltonian paths H
2. U_T as a symmetric function in the power-sum basis
3. Verifies principal specialization of U_T matches E_T(t)
4. Determines exact coefficients in U_T = f(p_1, 2p_3, 2p_5, ...)

Reference: Grinberg-Stanley, arXiv:2307.05569
  Theorem 1.30: U_D = sum_{sigma in S_V(D,Dbar)} (-1)^{phi(sigma)} p_{type(sigma)}
  Theorem 1.38: For tournaments, U_T is a polynomial in p_1, 2p_3, 2p_5, ...
                 with nonneg integer coefficients.

Key question: Are the coefficients in the power-sum expansion of U_T
exactly the independence numbers alpha_k of Omega(T)?
"""

from itertools import permutations, combinations
from collections import defaultdict, Counter
from fractions import Fraction
import sys


# ============================================================
# Tournament utilities
# ============================================================

def make_tournament(n, adj_dict):
    """Create adjacency matrix from dict of edges {(i,j): True means i->j}."""
    A = [[0]*n for _ in range(n)]
    for (i,j), val in adj_dict.items():
        if val:
            A[i][j] = 1
        else:
            A[j][i] = 1
    return A

def transitive_tournament(n):
    """T_n: transitive tournament, i->j iff i<j."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    return A

def cyclic_tournament(n):
    """C_n: cyclic tournament on n vertices (n must be odd).
    i->j iff (j-i) mod n in {1, 2, ..., (n-1)/2}."""
    assert n % 2 == 1
    A = [[0]*n for _ in range(n)]
    half = (n-1) // 2
    for i in range(n):
        for d in range(1, half+1):
            j = (i + d) % n
            A[i][j] = 1
    return A

def all_tournaments(n):
    """Generate all tournaments on n vertices."""
    edges = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(edges):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A

def tournament_name(A, n):
    """Try to identify a tournament."""
    T_trans = transitive_tournament(n)
    if A == T_trans:
        return f"T_{n} (transitive)"
    if n % 2 == 1:
        C = cyclic_tournament(n)
        if A == C:
            return f"C_{n} (cyclic)"
    return "unnamed"


# ============================================================
# Hamiltonian path enumeration with descent tracking
# ============================================================

def hamiltonian_paths(A, n):
    """Return list of all directed Hamiltonian paths as tuples."""
    paths = []
    for perm in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if A[perm[i]][perm[i+1]] != 1:
                valid = False
                break
        if valid:
            paths.append(perm)
    return paths

def descent_count(perm):
    """Count descents: positions i where perm[i] > perm[i+1]."""
    return sum(1 for i in range(len(perm)-1) if perm[i] > perm[i+1])

def eulerian_polynomial(A, n):
    """E_T(t) = sum_H t^{des(H)} as a dict {k: count}."""
    paths = hamiltonian_paths(A, n)
    poly = defaultdict(int)
    for p in paths:
        d = descent_count(p)
        poly[d] += 1
    return dict(poly)


# ============================================================
# Odd cycles and conflict graph
# ============================================================

def find_all_odd_cycles(A, n):
    """Find all directed odd cycles, returned as tuples (canonical form)."""
    cycles = []
    seen = set()
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts):
                valid = True
                for i in range(length):
                    if A[perm[i]][perm[(i+1) % length]] != 1:
                        valid = False
                        break
                if valid:
                    # Canonical: rotate so smallest vertex is first
                    min_idx = list(perm).index(min(perm))
                    canonical = perm[min_idx:] + perm[:min_idx]
                    if canonical not in seen:
                        seen.add(canonical)
                        cycles.append(canonical)
    return cycles

def conflict_graph_adj(cycles):
    """Build adjacency matrix of conflict graph Omega(T)."""
    m = len(cycles)
    adj = [[False]*m for _ in range(m)]
    cycle_sets = [set(c) for c in cycles]
    for i in range(m):
        for j in range(i+1, m):
            if cycle_sets[i] & cycle_sets[j]:
                adj[i][j] = adj[j][i] = True
    return adj

def independence_polynomial(cycles, return_coeffs=False):
    """Compute I(Omega(T), x) as polynomial coefficients [alpha_0, alpha_1, ...]."""
    m = len(cycles)
    adj = conflict_graph_adj(cycles)

    alpha = defaultdict(int)
    for mask in range(2**m):
        verts = [i for i in range(m) if (mask >> i) & 1]
        indep = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj[verts[i]][verts[j]]:
                    indep = False
                    break
            if not indep:
                break
        if indep:
            alpha[len(verts)] += 1

    max_k = max(alpha.keys()) if alpha else 0
    coeffs = [alpha.get(k, 0) for k in range(max_k + 1)]
    return coeffs


# ============================================================
# Grinberg-Stanley U_T computation
# ============================================================

def compute_phi(sigma, A, n):
    """
    Compute phi(sigma) = number of "descents" in the permutation sigma
    with respect to the tournament T.

    In Grinberg-Stanley, for a tournament T on V = {1,...,n}:
    A permutation sigma is in S_V(T, Tbar) means it's a permutation of V.

    phi(sigma) = number of i in {1,...,n-1} such that sigma(i) -> sigma(i+1)
    is NOT an arc of T (i.e., it's an arc of Tbar).

    Wait -- let me reconsider. In Grinberg-Stanley:
    - D is a digraph, Dbar is its complement
    - S_V(D, Dbar) = set of permutations sigma of V such that for each i,
      the arc (sigma(i), sigma(i+1)) is in D or Dbar (always true for tournaments)
    - phi(sigma) = #{i : (sigma(i), sigma(i+1)) in Dbar}

    For a tournament T, every arc is in T or Tbar, so S_V = S_n (all perms).
    phi(sigma) = #{i in {0,...,n-2} : T[sigma[i]][sigma[i+1]] = 0}
               = #{i : sigma[i+1] -> sigma[i] in T}
    """
    count = 0
    for i in range(n-1):
        if A[sigma[i]][sigma[i+1]] == 0:  # arc is in Tbar
            count += 1
    return count

def cycle_type(sigma):
    """Return the cycle type of permutation sigma as a sorted tuple (descending)."""
    n = len(sigma)
    visited = [False]*n
    lengths = []
    for i in range(n):
        if not visited[i]:
            length = 0
            j = i
            while not visited[j]:
                visited[j] = True
                j = sigma[j]
                length += 1
            lengths.append(length)
    return tuple(sorted(lengths, reverse=True))

def power_sum_monomial(partition):
    """Return a hashable representation of p_{lambda} = p_{l1} * p_{l2} * ..."""
    return tuple(sorted(partition, reverse=True))

def compute_UT_power_sum(A, n):
    """
    Compute U_T = sum_{sigma in S_n} (-1)^{phi(sigma)} p_{type(sigma)}

    Returns a dict: partition -> coefficient
    where partition is a tuple like (3,1,1) meaning p_3 * p_1^2
    """
    UT = defaultdict(int)
    for sigma in permutations(range(n)):
        phi = compute_phi(sigma, A, n)
        ct = cycle_type(sigma)
        sign = (-1)**phi
        UT[ct] += sign
    return dict(UT)

def principal_specialization_pk(k, n_vars, t):
    """
    p_k(1, t, t^2, ..., t^{n_vars-1}) = (1 - t^{n_vars*k}) / (1 - t^k)

    For t=1, this equals n_vars.
    We work with exact rational arithmetic using Fraction.
    Actually, let's just work with polynomial evaluation.
    """
    if t == 1:
        return n_vars
    return sum(t**(i*k) for i in range(n_vars))

def evaluate_UT_principal(UT, n_vars, t):
    """
    Evaluate U_T under the principal specialization x_i = t^{i-1} for i=1,...,n_vars.

    p_k -> p_k(1, t, ..., t^{n_vars-1}) = sum_{i=0}^{n_vars-1} t^{ik}
    """
    result = Fraction(0)
    for partition, coeff in UT.items():
        if coeff == 0:
            continue
        term = Fraction(coeff)
        for k in partition:
            pk_val = sum(Fraction(t)**( i*k) for i in range(n_vars))
            term *= pk_val
        result += term
    return result

def evaluate_UT_principal_poly(UT, n_vars):
    """
    Evaluate U_T under principal specialization, returning polynomial in t.

    p_k -> sum_{i=0}^{n_vars-1} t^{ik}

    Returns dict {power: coefficient}.
    """
    # Represent polynomials as dicts {exponent: coeff}
    def poly_mul(p1, p2):
        result = defaultdict(Fraction)
        for e1, c1 in p1.items():
            for e2, c2 in p2.items():
                result[e1+e2] += c1 * c2
        return dict(result)

    def poly_add(p1, p2):
        result = defaultdict(Fraction)
        for e, c in p1.items():
            result[e] += c
        for e, c in p2.items():
            result[e] += c
        return dict(result)

    def poly_scale(p, s):
        return {e: c*s for e, c in p.items()}

    def pk_poly(k, nv):
        """p_k(1, t, ..., t^{nv-1}) as polynomial."""
        return {i*k: Fraction(1) for i in range(nv)}

    total = defaultdict(Fraction)
    for partition, coeff in UT.items():
        if coeff == 0:
            continue
        # Product of p_{k_i} polynomials
        term = {0: Fraction(1)}
        for k in partition:
            term = poly_mul(term, pk_poly(k, n_vars))
        # Scale by coefficient
        for e, c in term.items():
            total[e] += Fraction(coeff) * c

    return dict(total)


def rewrite_in_odd_power_sums(UT, n):
    """
    Rewrite U_T in terms of p_1, 2p_3, 2p_5, ...

    By Theorem 1.38, for tournaments, U_T is a polynomial in p_1, 2p_3, 2p_5, ...
    with nonneg integer coefficients.

    Let q_k = 2*p_k for odd k >= 3, and q_1 = p_1.
    Then U_T = sum c_{lambda} prod q_{lambda_i} where c_{lambda} >= 0.

    This means: coefficient of p_{lambda} in U_T, where lambda has parts
    l_1 >= l_2 >= ... >= l_r (all odd), equals c_{lambda} * 2^{#{parts > 1}}.

    Actually, let's just display the p-expansion and check which partitions appear.
    """
    print(f"\n  Power-sum expansion of U_T:")
    for partition in sorted(UT.keys()):
        coeff = UT[partition]
        if coeff != 0:
            parts_str = " * ".join(f"p_{k}" for k in partition)
            print(f"    {coeff:+d} * {parts_str}")

    # Check: only odd parts should appear (Theorem 1.38 for tournaments)
    print(f"\n  Checking only odd-part partitions have nonzero coefficients...")
    all_odd = True
    for partition, coeff in UT.items():
        if coeff != 0:
            if any(k % 2 == 0 for k in partition):
                print(f"    VIOLATION: partition {partition} has even part, coeff = {coeff}")
                all_odd = False
    if all_odd:
        print(f"    PASS: all nonzero partitions have only odd parts")

    # Rewrite in terms of q_1 = p_1, q_k = 2*p_k (k odd, k>=3)
    print(f"\n  Rewriting in terms of q_1=p_1, q_k=2*p_k (odd k>=3):")
    for partition in sorted(UT.keys()):
        coeff = UT[partition]
        if coeff == 0:
            continue
        if any(k % 2 == 0 for k in partition):
            continue
        # p_{l1}*...*p_{lr} where each l_i is odd
        # = (q_1)^{a} * prod_{k>=3 odd} (q_k/2)^{b_k}
        # = prod q_i / 2^{number of parts >= 3}
        num_big = sum(1 for k in partition if k >= 3)
        q_coeff = Fraction(coeff, 2**num_big)
        parts_str = " * ".join(f"q_{k}" for k in partition)
        print(f"    {q_coeff} * {parts_str}")


# ============================================================
# Main analysis
# ============================================================

def analyze_tournament(A, n, name=""):
    print(f"\n{'='*70}")
    print(f"Tournament: {name} (n={n})")
    print(f"{'='*70}")

    # 1. Adjacency matrix
    print(f"\nAdjacency matrix:")
    for row in A:
        print(f"  {row}")

    # 2. Hamiltonian paths and E_T(t)
    paths = hamiltonian_paths(A, n)
    ET = eulerian_polynomial(A, n)
    print(f"\nH(T) = {len(paths)} Hamiltonian paths")
    print(f"E_T(t) = {' + '.join(f'{c}*t^{k}' for k,c in sorted(ET.items()))}")

    # 3. Odd cycles and independence polynomial
    cycles = find_all_odd_cycles(A, n)
    print(f"\nOdd cycles in T: {len(cycles)}")
    for c in cycles:
        print(f"  {c} (length {len(c)})")

    ip_coeffs = independence_polynomial(cycles)
    print(f"\nI(Omega(T), x) coefficients: {ip_coeffs}")
    ip_str = " + ".join(f"{c}*x^{k}" for k,c in enumerate(ip_coeffs) if c != 0)
    print(f"I(Omega(T), x) = {ip_str}")
    ip_at_2 = sum(c * 2**k for k, c in enumerate(ip_coeffs))
    print(f"I(Omega(T), 2) = {ip_at_2}")
    print(f"H(T) = {len(paths)}, match: {ip_at_2 == len(paths)}")

    # 4. Grinberg-Stanley U_T
    print(f"\n--- Grinberg-Stanley U_T computation ---")
    UT = compute_UT_power_sum(A, n)
    rewrite_in_odd_power_sums(UT, n)

    # 5. Principal specialization check
    print(f"\n--- Principal specialization verification ---")
    # With n variables: p_k(1,t,...,t^{n-1}) = (1-t^{nk})/(1-t^k)
    # E_T(t) should equal (1/n!) * U_T evaluated at this specialization?
    # No -- Grinberg-Stanley Theorem 1.30 gives U_T directly as a symmetric function
    # whose principal specialization should give... let me think.
    #
    # Actually, the relationship is:
    # U_T = sum_{sigma} (-1)^{phi(sigma)} p_{type(sigma)}
    #
    # The key identity (Corollary 1.36 of Grinberg-Stanley):
    # If we apply the "standard" specialization to n variables x_1,...,x_n,
    # the coefficient of x_1 * x_2 * ... * x_n in U_T equals
    # sum_{sigma} (-1)^{phi(sigma)} = sum_H (-1)^{des(H)} ...
    #
    # Wait, let me reconsider. Their U_T is a symmetric function.
    # Corollary 1.36: The "monomial extraction" [x_1 x_2 ... x_n] U_T
    # = sum_{sigma in S_n} (-1)^{phi(sigma)} = ?
    #
    # Actually, the simpler relationship for tournaments:
    # Under x_i = t^{i-1}, the symmetric function ps(U_T)(t) should give
    # a polynomial related to E_T(t).
    #
    # Let's just compute and see what happens.

    UT_poly = evaluate_UT_principal_poly(UT, n)

    # Clean up near-zero coefficients
    UT_poly_clean = {e: c for e, c in UT_poly.items() if c != 0}
    max_deg = max(UT_poly_clean.keys()) if UT_poly_clean else 0

    print(f"\n  U_T under x_i = t^{{i-1}} (n={n} variables):")
    for e in range(max_deg + 1):
        c = UT_poly_clean.get(e, Fraction(0))
        if c != 0:
            print(f"    t^{e}: {c}")

    # Compare with E_T(t)
    print(f"\n  E_T(t) = sum_H t^{{des(H)}}:")
    for k in sorted(ET.keys()):
        print(f"    t^{k}: {ET[k]}")

    # Check if they match up to a scalar or normalization
    # The principal specialization of p_lambda with n variables gives
    # the evaluation of the symmetric function at (1, t, t^2, ..., t^{n-1})
    # For the identity permutation (type (1,1,...,1)): p_{(1^n)} = (sum x_i)^n
    # but that's not right either.
    #
    # Let me check: does ps(U_T)/n! = E_T(t)?  Or some other normalization?

    print(f"\n  Checking ratios ps(U_T)(t^k) / E_T(t^k):")
    # Evaluate at a few numerical points
    for t_val in [Fraction(1,3), Fraction(1,2), Fraction(2,1)]:
        ps_val = Fraction(0)
        for e, c in UT_poly_clean.items():
            ps_val += c * t_val**e
        et_val = sum(Fraction(c) * t_val**k for k, c in ET.items())
        if et_val != 0:
            ratio = ps_val / et_val
            print(f"    t={t_val}: ps(U_T)={ps_val}, E_T={et_val}, ratio={ratio} = {float(ratio):.6f}")
        else:
            print(f"    t={t_val}: ps(U_T)={ps_val}, E_T={et_val}, E_T=0!")

    return UT, ET, ip_coeffs, UT_poly_clean


def main():
    print("=" * 70)
    print("GRINBERG-STANLEY U_T COMPUTATION")
    print("=" * 70)

    # ---- n=3 ----
    print("\n\n" + "#"*70)
    print("# n = 3")
    print("#"*70)

    # Transitive tournament T_3: 0->1, 0->2, 1->2
    T3 = transitive_tournament(3)
    UT3, ET3, IP3, PS3 = analyze_tournament(T3, 3, "T_3 (transitive)")

    # 3-cycle C_3: 0->1, 1->2, 2->0
    C3 = cyclic_tournament(3)
    UTC3, ETC3, IPC3, PSC3 = analyze_tournament(C3, 3, "C_3 (3-cycle)")

    # ---- n=5 ----
    print("\n\n" + "#"*70)
    print("# n = 5")
    print("#"*70)

    T5 = transitive_tournament(5)
    UT5, ET5, IP5, PS5 = analyze_tournament(T5, 5, "T_5 (transitive)")

    C5 = cyclic_tournament(5)
    UTC5, ETC5, IPC5, PSC5 = analyze_tournament(C5, 5, "C_5 (cyclic)")

    # ---- All n=3 tournaments ----
    print("\n\n" + "#"*70)
    print("# ALL TOURNAMENTS n=3 (exhaustive)")
    print("#"*70)

    seen = set()
    for A in all_tournaments(3):
        key = tuple(tuple(row) for row in A)
        if key not in seen:
            seen.add(key)
            name = tournament_name(A, 3)
            analyze_tournament(A, 3, name)

    # ---- Key comparison ----
    print("\n\n" + "#"*70)
    print("# SUMMARY: U_T coefficients vs I(Omega, x) coefficients")
    print("#"*70)

    print("""
For the central question: are the coefficients of U_T in the
{p_1, 2p_3, 2p_5, ...} basis related to the independence polynomial
I(Omega(T), x)?

Specifically, if U_T = sum_{lambda} c_lambda * q_{lambda_1} * ... * q_{lambda_r}
where q_1 = p_1, q_k = 2*p_k (odd k>=3), and c_lambda >= 0,

does I(Omega(T), x) = sum_k alpha_k * x^k relate to these c_lambda?
""")

    # Additional analysis: compute U_T for ALL n=3 and n=5 tournaments
    # and compare the q-expansion coefficients with independence polynomial

    print("\n--- Exhaustive n=3 comparison ---")
    for A in all_tournaments(3):
        n = 3
        cycles = find_all_odd_cycles(A, n)
        ip = independence_polynomial(cycles)
        UT = compute_UT_power_sum(A, n)

        # Extract the coefficient of p_(3) = q_3/2
        p3_coeff = UT.get((3,), 0)
        p111_coeff = UT.get((1,1,1), 0)

        # q_3 coeff = p3_coeff / 2  (since p_3 = q_3/2)
        # Wait: q_3 = 2*p_3, so p_3 = q_3/2, meaning coeff of q_3 = p3_coeff/2

        ip_at_2 = sum(c * 2**k for k, c in enumerate(ip))
        paths = hamiltonian_paths(A, n)

        print(f"  A={[row for row in A]}")
        print(f"    H(T)={len(paths)}, I(Omega,2)={ip_at_2}")
        print(f"    I(Omega,x) coeffs = {ip}")
        print(f"    U_T: p_(1,1,1) coeff = {p111_coeff}, p_(3) coeff = {p3_coeff}")
        # alpha_1 = number of odd cycles = ip[1] if len(ip)>1 else 0
        alpha_1 = ip[1] if len(ip) > 1 else 0
        print(f"    alpha_1 (# odd cycles) = {alpha_1}")
        print(f"    p_(3) coeff / 2 = {Fraction(p3_coeff, 2)} (should be alpha_1 = {alpha_1}?)")
        print()


if __name__ == "__main__":
    main()
