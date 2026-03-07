#!/usr/bin/env python3
"""
Bridge between Grinberg-Stanley's U_T and our OCF / G_T(t,x).

Grinberg-Stanley (2307.05569) Theorem 1.30:
  U_T = sum_{sigma in S_V(T,T-bar)} (-1)^{phi(sigma)} p_{type(sigma)}

For tournaments (Thm 1.38): U_T is a polynomial in p_1, 2p_3, 2p_5, ...
with nonneg integer coefficients.

OUR OCF: H(T) = I(Omega(T), 2) = sum_S 2^{|S|} over independent sets S in Omega(T).

QUESTION: Does principal specialization of U_T recover E_T(t)?
QUESTION: Are GS coefficients exactly the alpha_k of Omega(T)?

We compute U_T directly from the definition for small tournaments and compare.

opus-2026-03-07-S34
"""
from itertools import permutations, combinations
from collections import Counter, defaultdict
from math import factorial, comb
from fractions import Fraction
import sympy
from sympy import symbols, Poly, expand, Rational

def tournament_from_bits(n, bits):
    """Build adjacency matrix from bit encoding."""
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits[idx] == 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def all_tournaments(n):
    """Generate all tournaments on n vertices."""
    m = n*(n-1)//2
    for bits in range(2**m):
        b = [(bits >> k) & 1 for k in range(m)]
        yield tournament_from_bits(n, b)

def hamiltonian_paths_with_descents(A, n):
    """Return list of (path, descent_count) for all Hamiltonian paths."""
    results = []
    for perm in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if A[perm[i]][perm[i+1]] != 1:
                valid = False
                break
        if valid:
            desc = sum(1 for i in range(n-1) if perm[i] > perm[i+1])
            results.append((perm, desc))
    return results

def descent_polynomial(A, n):
    """Compute E_T(t) = sum_H t^{des(H)}."""
    t = symbols('t')
    paths = hamiltonian_paths_with_descents(A, n)
    poly = sum(t**d for _, d in paths)
    return expand(poly)

def find_directed_odd_cycles(A, n):
    """Find all directed odd cycles in tournament A."""
    cycles = []
    for length in range(3, n+1, 2):
        for subset in combinations(range(n), length):
            for perm in permutations(subset):
                # Check if this is a directed cycle
                is_cycle = True
                for i in range(length):
                    if A[perm[i]][perm[(i+1)%length]] != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    # Canonical form: start from smallest vertex, check direction
                    min_idx = perm.index(min(perm))
                    canonical = tuple(perm[min_idx:] + perm[:min_idx])
                    if canonical not in [c for c in cycles]:
                        cycles.append(canonical)
                        break  # Each ordered cycle found once per rotation
    # Deduplicate
    unique = set()
    result = []
    for c in cycles:
        min_idx = list(c).index(min(c))
        canon = c[min_idx:] + c[:min_idx]
        if canon not in unique:
            unique.add(canon)
            result.append(canon)
    return result

def conflict_graph(cycles, n):
    """Build conflict graph: edge if cycles share a vertex."""
    nc = len(cycles)
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if set(cycles[i]) & set(cycles[j]):
                adj[i][j] = adj[j][i] = True
    return adj

def independent_sets(adj, nc):
    """Find all independent sets in graph with adjacency matrix adj."""
    result = []
    for mask in range(2**nc):
        nodes = [i for i in range(nc) if (mask >> i) & 1]
        is_indep = True
        for i in range(len(nodes)):
            for j in range(i+1, len(nodes)):
                if adj[nodes[i]][nodes[j]]:
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            result.append(nodes)
    return result

def independence_polynomial(cycles, n):
    """Compute I(Omega(T), x)."""
    if not cycles:
        return {0: 1}
    adj = conflict_graph(cycles, n)
    nc = len(cycles)
    indep = independent_sets(adj, nc)
    poly = defaultdict(int)
    for s in indep:
        poly[len(s)] += 1
    return dict(poly)

def compute_UT_permutation_expansion(A, n):
    """
    Compute U_T via Grinberg-Stanley's Theorem 1.30.

    U_T = sum_{sigma in S_V(T, T-bar)} (-1)^{phi(sigma)} p_{type(sigma)}

    S_V(T, T-bar) = permutations whose cycles are all either:
      - directed cycles of T, or
      - directed cycles of the complement T-bar

    phi(sigma) = sum over T-cycles gamma of (len(gamma) - 1)

    For a tournament T, the complement T-bar is T^op (opposite tournament).
    """
    T_bar = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                T_bar[i][j] = 1 - A[i][j]

    # For each permutation, check if every cycle is either a T-cycle or T-bar-cycle
    result = defaultdict(int)  # partition -> coefficient

    for perm in permutations(range(n)):
        # Decompose into cycles
        visited = [False]*n
        cycles = []
        for i in range(n):
            if not visited[i]:
                cycle = []
                j = i
                while not visited[j]:
                    visited[j] = True
                    cycle.append(j)
                    j = perm[j]
                cycles.append(tuple(cycle))

        # For each cycle, check if it's a T-cycle or T-bar-cycle
        valid = True
        phi = 0
        for cyc in cycles:
            k = len(cyc)
            if k == 1:
                # Fixed points are always valid (trivial 1-cycles)
                # They contribute k-1 = 0 to phi regardless
                continue
            # Check if directed cycle in T
            is_T_cycle = all(A[cyc[i]][cyc[(i+1)%k]] == 1 for i in range(k))
            # Check if directed cycle in T-bar
            is_Tbar_cycle = all(T_bar[cyc[i]][cyc[(i+1)%k]] == 1 for i in range(k))

            if is_T_cycle:
                phi += k - 1
            elif is_Tbar_cycle:
                pass  # T-bar cycle contributes 0 to phi
            else:
                valid = False
                break

        if valid:
            # Get cycle type (partition)
            cycle_type = tuple(sorted([len(c) for c in cycles], reverse=True))
            sign = (-1)**phi
            result[cycle_type] += sign

    return dict(result)

def partition_to_pstring(partition):
    """Convert partition to p-monomial string like p_3^2*p_1."""
    cnt = Counter(partition)
    parts = []
    for k in sorted(cnt.keys(), reverse=True):
        if cnt[k] == 1:
            parts.append(f"p_{k}")
        else:
            parts.append(f"p_{k}^{cnt[k]}")
    return "*".join(parts)

def main():
    print("=" * 70)
    print("GRINBERG-STANLEY U_T <-> OCF BRIDGE")
    print("=" * 70)

    # n=3: all tournaments
    print("\n--- n=3 ---")
    n = 3
    for bits in range(2**(n*(n-1)//2)):
        b = [(bits >> k) & 1 for k in range(n*(n-1)//2)]
        A = tournament_from_bits(n, b)

        # Compute U_T
        ut = compute_UT_permutation_expansion(A, n)

        # Compute OCF
        cycles = find_directed_odd_cycles(A, n)
        ip = independence_polynomial(cycles, n)

        # Compute E_T(t)
        et = descent_polynomial(A, n)

        # Compute H(T)
        H = sum(2**k * ip.get(k, 0) for k in ip)

        print(f"\nTournament bits={b}, H={H}, #3-cycles={len(cycles)}")
        print(f"  I(Omega, x) = {ip}")
        print(f"  E_T(t) = {et}")
        print(f"  U_T expansion:")
        for part, coeff in sorted(ut.items()):
            if coeff != 0:
                print(f"    {coeff:+d} * {partition_to_pstring(part)}")

    # n=5: sample tournaments
    print("\n\n--- n=5 (first 16 tournaments) ---")
    n = 5
    m = n*(n-1)//2
    for bits_int in range(16):
        b = [(bits_int >> k) & 1 for k in range(m)]
        A = tournament_from_bits(n, b)

        ut = compute_UT_permutation_expansion(A, n)
        cycles = find_directed_odd_cycles(A, n)
        ip = independence_polynomial(cycles, n)
        H = sum(2**k * ip.get(k, 0) for k in ip)

        print(f"\nbits={bits_int:05b}, H={H}")
        print(f"  I(Omega,x) = {ip}")
        print(f"  U_T (nonzero terms):")
        for part, coeff in sorted(ut.items()):
            if coeff != 0:
                print(f"    {coeff:+d} * {partition_to_pstring(part)}")

    # Key test: convert U_T to p_1, 2p_3, 2p_5 basis
    print("\n\n--- KEY TEST: U_T in {p_1, 2p_3, 2p_5} basis ---")
    print("For tournaments, Thm 1.38 says U_T = poly in p_1, 2p_3, 2p_5,... with nonneg coeffs")
    print("We check: are these coefficients = independence polynomial alpha_k?")

    n = 3
    for bits in range(2**(n*(n-1)//2)):
        b = [(bits >> k) & 1 for k in range(n*(n-1)//2)]
        A = tournament_from_bits(n, b)
        ut = compute_UT_permutation_expansion(A, n)
        cycles = find_directed_odd_cycles(A, n)
        ip = independence_polynomial(cycles, n)

        # For n=3: possible cycle types are (1,1,1), (3,)
        # p_{(1,1,1)} = p_1^3, p_{(3)} = p_3
        # In the {p_1, 2p_3} basis: U_T = a * p_1^3 + b * (2p_3)
        # So coeff of p_3 should be b = coeff_p3 / 2? No...
        # Actually: U_T = a*p_1^3 + c*(2p_3), so c = coeff_p3 / 2

        print(f"\nbits={b}, I(Omega,x)={ip}")
        coeff_p1_cubed = ut.get((1,1,1), 0)
        coeff_p3 = ut.get((3,), 0)
        print(f"  U_T = {coeff_p1_cubed}*p_1^3 + {coeff_p3}*p_3")
        print(f"  In (p_1, 2p_3) basis: {coeff_p1_cubed}*p_1^3 + {coeff_p3//2}*(2p_3)")
        print(f"  alpha_0={ip.get(0,0)}, alpha_1={ip.get(1,0)}")

    # n=5
    print("\n\n--- n=5 KEY TEST ---")
    n = 5
    m = n*(n-1)//2
    for bits_int in [0, 1, 5, 15, 31, 100, 200, 500]:
        if bits_int >= 2**m:
            continue
        b = [(bits_int >> k) & 1 for k in range(m)]
        A = tournament_from_bits(n, b)
        ut = compute_UT_permutation_expansion(A, n)
        cycles = find_directed_odd_cycles(A, n)
        ip = independence_polynomial(cycles, n)
        H = sum(2**k * ip.get(k, 0) for k in ip)

        print(f"\nbits={bits_int}, H={H}, I(Omega,x)={ip}")
        # For n=5, cycle types involving only odd parts:
        # (5), (3,1,1), (1,1,1,1,1)
        # Could also have (3,1,1) from mixed

        # Extract key coefficients
        for part in sorted(ut.keys()):
            c = ut[part]
            if c != 0:
                # Check if all cycle lengths are odd
                all_odd = all(p % 2 == 1 for p in part)
                print(f"  {c:+d} * {partition_to_pstring(part)} {'[all odd]' if all_odd else '[has even]'}")

    # Principal specialization test
    print("\n\n--- PRINCIPAL SPECIALIZATION ---")
    print("Under x_i = t^{i-1}: p_k = sum t^{k*i} = (1-t^{nk})/(1-t^k)")
    print("Does specializing U_T give E_T(t)?")

    t = symbols('t')
    n = 3
    for bits in range(2**(n*(n-1)//2)):
        b = [(bits >> k) & 1 for k in range(n*(n-1)//2)]
        A = tournament_from_bits(n, b)
        ut = compute_UT_permutation_expansion(A, n)
        et = descent_polynomial(A, n)

        # Specialize: p_k -> sum_{i=0}^{n-1} t^{ki} = (1-t^{nk})/(1-t^k) for k>1, p_1 -> n for k=1... no
        # Actually for n variables: p_k(1, t, t^2, ..., t^{n-1}) = 1 + t^k + t^{2k} + ... + t^{(n-1)k}
        spec = 0
        for part, coeff in ut.items():
            term = coeff
            for k in part:
                pk_val = sum(t**(k*i) for i in range(n))
                term *= pk_val
            spec += term
        spec = expand(spec)

        # BUT: U_T uses quasisymmetric functions and the connection to E_T might need
        # a different specialization. Let's check.
        print(f"\nbits={b}")
        print(f"  E_T(t) = {et}")
        print(f"  spec(U_T) = {spec}")
        print(f"  Match: {expand(et - spec) == 0}")

if __name__ == "__main__":
    main()
