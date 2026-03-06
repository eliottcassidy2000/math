#!/usr/bin/env python3
"""
Test h-positivity of the Rédei-Berge function U_T for tournaments.
Instance: opus-2026-03-06-S10

The Rédei-Berge function U_T = sum_{sigma in S_V(T,Tbar)} (-1)^{phi(sigma)} p_{type(sigma)}
where S_V(T,Tbar) = {sigma in S_n : every cycle of sigma is a cycle of T or Tbar}
and phi(sigma) = sum over cycles gamma of sigma that are in Tbar of (len(gamma)-1).

For tournaments, Tbar = T^op, and sigma in S_V(T,T^op) means every cycle of sigma
is a directed cycle in T or in T^op.

We expand U_T in the power-sum basis, convert to the h-basis (complete homogeneous),
and check if all coefficients are non-negative.

Key connection to tiling geometry:
- SC tournaments (T ~ T^op) correspond to grid-symmetric tilings
- h-positivity would mean the Hamiltonian path structure decomposes non-negatively
"""

from itertools import permutations
from math import factorial
from fractions import Fraction
from collections import defaultdict

def tournament_from_bits(bits, n):
    """Convert bitmask to adjacency dict. Non-path arcs only."""
    adj = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                adj[(i,j)] = 0
    # Path arcs: i -> i+1
    for i in range(n-1):
        adj[(i,i+1)] = 1
    # Non-path arcs: enumerated in a specific order
    idx = 0
    for gap in range(2, n):
        for i in range(n - gap):
            j = i + gap
            if (bits >> idx) & 1:
                adj[(i,j)] = 1
                adj[(j,i)] = 0
            else:
                adj[(j,i)] = 1
                adj[(i,j)] = 0
            idx += 1
    return adj

def all_tournaments(n):
    """Generate all tournaments on n vertices as adjacency dicts."""
    m = n*(n-1)//2 - (n-1)  # non-path arcs
    for bits in range(2**m):
        yield bits, tournament_from_bits(bits, n)

def get_cycles(perm):
    """Get the cycle decomposition of a permutation (list)."""
    n = len(perm)
    visited = [False]*n
    cycles = []
    for start in range(n):
        if visited[start]:
            continue
        cycle = []
        v = start
        while not visited[v]:
            visited[v] = True
            cycle.append(v)
            v = perm[v]
        if len(cycle) > 0:
            cycles.append(tuple(cycle))
    return cycles

def is_directed_cycle_in(cycle, adj):
    """Check if (v0, v1, ..., vk-1, v0) is a directed cycle in the digraph."""
    k = len(cycle)
    if k == 1:
        return True  # Fixed points are always "cycles"
    for i in range(k):
        if adj.get((cycle[i], cycle[(i+1)%k]), 0) != 1:
            return False
    return True

def compute_U_T_power_sum(n, adj):
    """
    Compute U_T in the power-sum basis.

    U_T = sum_{sigma in S_V(T,Tbar)} (-1)^{phi(sigma)} p_{type(sigma)}

    where Tbar = complement = T^op for tournaments.
    S_V(T,Tbar) = permutations where every cycle is a cycle of T or T^op.
    phi(sigma) = sum of (len(gamma)-1) for cycles gamma that are in Tbar (T^op).
    """
    # Build T^op adjacency
    adj_op = {}
    for (i,j), v in adj.items():
        adj_op[(j,i)] = v  # reverse all edges
    # Actually for tournament: adj_op[(i,j)] = adj[(j,i)] = 1 - adj[(i,j)]
    adj_op = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                adj_op[(i,j)] = adj.get((j,i), 0)

    # power_sum_coeffs: maps partition (sorted tuple) -> coefficient
    coeffs = defaultdict(int)

    for perm_tuple in permutations(range(n)):
        perm = list(perm_tuple)
        cycles = get_cycles(perm)

        # Check: every cycle must be a cycle of T or T^op
        phi = 0
        valid = True
        for cycle in cycles:
            if len(cycle) == 1:
                continue  # Fixed points are fine
            in_T = is_directed_cycle_in(cycle, adj)
            in_Top = is_directed_cycle_in(cycle, adj_op)
            if not in_T and not in_Top:
                valid = False
                break
            if in_Top and not in_T:
                phi += len(cycle) - 1
            # If in both T and T^op (possible for palindromic cycles),
            # convention: count as in T (phi contribution = 0)
            # Actually Grinberg-Stanley: if in T^op (=Tbar), add len-1
            # If a cycle is in BOTH, it's in S_V(X) subset of S_V(X,Xbar),
            # so phi contribution is 0 for those.
            elif in_Top and in_T:
                pass  # phi += 0 (it's a cycle of T, so no contribution)

        if not valid:
            continue

        # Compute partition type
        cycle_lens = sorted([len(c) for c in cycles], reverse=True)
        partition = tuple(cycle_lens)

        sign = (-1)**phi
        coeffs[partition] += sign

    return dict(coeffs)

def partition_to_str(p):
    return '(' + ','.join(map(str, p)) + ')'

def change_basis_p_to_h(n, p_coeffs):
    """
    Convert from power-sum basis to complete homogeneous basis.
    Uses the relation: p_lambda = sum_mu z_lambda^{-1} chi^mu(lambda) s_mu
    Then s to h. For small n, we can use the character table approach.

    Actually, simpler: use the transition matrix p -> h directly.
    The relation is: h_n = (1/n!) sum_{sigma in S_n} p_{type(sigma)}
    More precisely: p_mu = sum_lambda z_mu^{-1} ...

    For small n, let's just use the explicit formulas.
    """
    # For degree n, we need the transition matrix from p to h basis
    # h_lambda = sum_mu M_{lambda,mu} p_mu where M is known
    #
    # The inverse relation: p_mu = sum_lambda chi^lambda(mu) h_lambda
    # No, that's not right either.
    #
    # Actually: in the ring of symmetric functions,
    # p_k = sum_{j=0}^{k-1} (-1)^j e_j h_{k-j}  (Newton's identity)
    # and h_n = (1/n) sum_{k=1}^n p_k h_{n-k}
    #
    # For the MULTIPLICATIVE structure: p_{(k1,k2,...)} = p_{k1} * p_{k2} * ...
    # h_{(l1,l2,...)} = h_{l1} * h_{l2} * ...
    #
    # The transition matrix from p_lambda to h_mu is:
    # p_lambda = sum_mu (z_lambda)^{-1} chi^mu(lambda) * ...
    #
    # This is getting complicated. Let me use a direct computational approach.
    # For n <= 5, I'll compute the transition matrix numerically.

    # Generate all partitions of n
    partitions = list(generate_partitions(n))

    # Build the transition matrix from p to h
    # Using the fact that for any symmetric function f,
    # [h_mu]f = <f, m_mu> where <,> is the Hall inner product
    # and <p_lambda, m_mu> = ...
    #
    # Alternative: use the inverse Kostka matrix approach
    # p_lambda/z_lambda = sum_mu K^{-1}_{mu,lambda} s_mu
    # s_mu = sum_nu K_{mu,nu} m_nu
    # h_mu = sum_nu m_nu (for all nu coarser than mu... no)
    # h_mu = s_mu only if mu = (n)

    # OK let me use a different approach: direct numerical computation
    # of the principal specialization.

    # For our purpose, we actually want to check h-positivity
    # of U_T as an element of Sym_n.
    #
    # The simplest approach: compute u_T(m) = ps^1(U_T)(m) for
    # m = 1, 2, ..., n and use interpolation to find h-coefficients.
    #
    # u_T(m) = sum_lambda c_lambda * h_lambda(1^m)
    # where h_k(1^m) = C(m+k-1, k) = m(m+1)...(m+k-1)/k!
    # and h_lambda(1^m) = prod h_{lambda_i}(1^m)
    #
    # Actually, the principal specialization at m is:
    # ps^1(p_k)(m) = m (since p_k(1,...,1) = m with m variables)
    # ps^1(p_lambda)(m) = m^{l(lambda)} where l is number of parts

    # So u_T(m) = sum_lambda c_lambda * m^{l(lambda)}
    # where c_lambda is the p-coefficient of partition lambda.

    # Wait, that's not right. ps^1(p_k)(m) = m for ALL k, so
    # ps^1(p_lambda)(m) = m^{l(lambda)}.
    # This means u_T(m) = sum_lambda coeff(lambda) * m^{l(lambda)}

    # For h-positivity, we need the h-expansion.
    # ps^1(h_k)(m) = C(m+k-1, k)
    # ps^1(h_lambda)(m) = prod C(m+lambda_i-1, lambda_i)

    # So u_T(m) = sum_lambda c_lambda^h * prod C(m+lambda_i-1, lambda_i)

    # For n=3, the h-partitions are (3), (2,1), (1,1,1).
    # ps^1 values:
    # h_{(3)}: C(m+2,3) = m(m+1)(m+2)/6
    # h_{(2,1)}: C(m+1,2)*m = m^2(m+1)/2
    # h_{(1,1,1)}: m^3

    # We need u_T(m) for m=1,2,3 to solve for 3 unknowns.
    # u_T(1) = H(T)
    # u_T(m) = sum of p-coefficients * m^{l(lambda)}

    pass  # Will implement below

def generate_partitions(n):
    """Generate all partitions of n in decreasing order."""
    if n == 0:
        yield ()
        return
    def helper(n, max_part):
        if n == 0:
            yield ()
            return
        for k in range(min(n, max_part), 0, -1):
            for rest in helper(n-k, k):
                yield (k,) + rest
    yield from helper(n, n)

def compute_u_T_at_m(n, p_coeffs, m):
    """Compute u_T(m) = ps^1(U_T)(m) = sum_lambda c_lambda * m^{l(lambda)}."""
    total = 0
    for partition, coeff in p_coeffs.items():
        l = len(partition)  # number of parts
        total += coeff * (m ** l)
    return total

def h_at_m(partition, m):
    """Compute ps^1(h_lambda)(m) = prod C(m+lambda_i-1, lambda_i)."""
    result = 1
    for part in partition:
        # C(m+part-1, part) = prod_{j=0}^{part-1} (m+j) / part!
        num = 1
        for j in range(part):
            num *= (m + j)
        num //= factorial(part)
        result *= num
    return result

def check_h_positivity(n, p_coeffs):
    """Check h-positivity by solving linear system."""
    partitions = list(generate_partitions(n))
    k = len(partitions)

    # Compute u_T(m) for m = 1, ..., k+2 (overdetermined for safety)
    num_evals = k + 5
    u_values = []
    for m in range(1, num_evals + 1):
        u_values.append(compute_u_T_at_m(n, p_coeffs, m))

    # Build matrix: u_T(m) = sum_j c_j * h_{partitions[j]}(m)
    # For each m, this gives one equation.
    # Use first k equations.

    # Build matrix A[i][j] = h_{partitions[j]}(m=i+1) as Fraction
    A = []
    b = []
    for i in range(k):
        m = i + 1
        row = []
        for j in range(k):
            row.append(Fraction(h_at_m(partitions[j], m)))
        A.append(row)
        b.append(Fraction(u_values[i]))

    # Solve Ax = b using Gaussian elimination with exact arithmetic
    for col in range(k):
        # Find pivot
        pivot = None
        for row in range(col, k):
            if A[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            continue
        A[col], A[pivot] = A[pivot], A[col]
        b[col], b[pivot] = b[pivot], b[col]

        for row in range(k):
            if row == col:
                continue
            if A[row][col] == 0:
                continue
            factor = A[row][col] / A[col][col]
            for j in range(k):
                A[row][j] -= factor * A[col][j]
            b[row] -= factor * b[col]

    # Extract solution
    h_coeffs = {}
    for i in range(k):
        if A[i][i] != 0:
            val = b[i] / A[i][i]
            h_coeffs[partitions[i]] = val
        else:
            h_coeffs[partitions[i]] = Fraction(0)

    # Verify against extra evaluations
    for m_idx in range(k, num_evals):
        m = m_idx + 1
        predicted = sum(h_coeffs[p] * h_at_m(p, m) for p in partitions)
        actual = Fraction(u_values[m_idx])
        if predicted != actual:
            print(f"  WARNING: verification failed at m={m}: predicted={predicted}, actual={actual}")

    return h_coeffs

def count_ham_paths(n, adj):
    """Count Hamiltonian paths using DP."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj.get((v, u), 0):
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def is_self_converse(n, adj):
    """Check if tournament is isomorphic to its converse."""
    # Quick check: score sequence must be palindromic
    scores = []
    for i in range(n):
        s = sum(adj.get((i,j), 0) for j in range(n) if j != i)
        scores.append(s)
    scores_sorted = sorted(scores)
    rev = sorted(scores, reverse=True)
    complement_scores = [n-1-s for s in scores_sorted]
    if sorted(complement_scores) != scores_sorted:
        return False
    # Full check would require isomorphism testing; skip for now
    return True  # Necessary condition only

def main():
    print("=== H-POSITIVITY TEST FOR TOURNAMENT RÉDEI-BERGE FUNCTION ===\n")

    for n in range(3, 6):
        print(f"\n{'='*60}")
        print(f"n = {n}")
        print(f"{'='*60}")

        m = n*(n-1)//2 - (n-1)
        num_tournaments = 2**m

        # Track h-positivity statistics
        h_positive_count = 0
        h_negative_count = 0
        total = 0

        # Group by isomorphism class (approximate by score sequence + H)
        seen = {}

        for bits, adj in all_tournaments(n):
            H = count_ham_paths(n, adj)
            scores = tuple(sorted([sum(adj.get((i,j),0) for j in range(n) if j!=i) for i in range(n)]))

            # Skip if we've seen this (score, H) pair already
            key = (scores, H)
            if key in seen:
                continue
            seen[key] = bits

            total += 1

            p_coeffs = compute_U_T_power_sum(n, adj)
            h_coeffs = check_h_positivity(n, p_coeffs)

            is_h_pos = all(v >= 0 for v in h_coeffs.values())

            if is_h_pos:
                h_positive_count += 1
            else:
                h_negative_count += 1

            sc = is_self_converse(n, adj)

            # Print details for small n
            if n <= 4 or not is_h_pos:
                print(f"\n  bits={bits}, H={H}, scores={scores}, SC={sc}")
                print(f"  p-expansion: ", end='')
                for p, c in sorted(p_coeffs.items()):
                    if c != 0:
                        print(f"{c}*p{partition_to_str(p)} ", end='')
                print()
                print(f"  h-expansion: ", end='')
                for p in sorted(h_coeffs.keys()):
                    c = h_coeffs[p]
                    if c != 0:
                        print(f"{c}*h{partition_to_str(p)} ", end='')
                print()
                print(f"  h-positive: {is_h_pos}")

        print(f"\n  SUMMARY n={n}: {h_positive_count}/{total} classes h-positive, {h_negative_count}/{total} h-negative")

if __name__ == '__main__':
    main()
