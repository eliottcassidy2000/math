"""
Discriminant approach to real-rootedness of I(Omega(T), x).

For n <= 8, alpha(Omega(T)) <= floor(n/3) <= 2, so
  I(Omega(T), x) = 1 + a1*x + a2*x^2
has real roots iff a1^2 >= 4*a2.

This script verifies the discriminant condition exhaustively at n=6
and by sampling at n=7,8, providing an ELEMENTARY alternative to
Chudnovsky-Seymour for proving real-rootedness at n <= 8.

For n=9 (alpha <= 3), checks degree-3 discriminant conditions.

Author: opus-2026-03-06-S15
"""

import itertools
import numpy as np
from collections import defaultdict

def tournament_from_bits(n, bits):
    """Build adjacency matrix from upper-triangular bits."""
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits[idx]:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def find_all_directed_cycles(A, n, length):
    """Find all directed cycles of given length in tournament with adjacency matrix A.
    Fix minimum vertex to avoid counting each cycle multiple times."""
    cycles = []
    for subset in itertools.combinations(range(n), length):
        min_v = subset[0]  # fix first vertex
        rest = list(subset[1:])
        for perm in itertools.permutations(rest):
            path = [min_v] + list(perm)
            valid = True
            for k in range(length):
                if A[path[k]][path[(k+1) % length]] != 1:
                    valid = False
                    break
            if valid:
                cycles.append(frozenset(subset))
    return cycles

def find_all_odd_cycles(A, n):
    """Find all odd directed cycles (length 3, 5, 7, ...) in tournament."""
    all_cycles = []
    for length in range(3, n+1, 2):
        cycles = find_all_directed_cycles(A, n, length)
        all_cycles.extend(cycles)
    return all_cycles

def compute_independence_poly(cycles):
    """Compute I(Omega, x) coefficients: alpha_k = #{independent sets of size k}.
    Omega has cycles as vertices, edges between cycles sharing a vertex.
    Independent sets = collections of pairwise vertex-disjoint cycles."""
    m = len(cycles)
    if m == 0:
        return [1]  # I = 1

    # Build adjacency (conflict): cycles[i] and cycles[j] conflict if they share a vertex
    conflict = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cycles[i] & cycles[j]:  # shared vertex
                conflict[i][j] = conflict[j][i] = True

    # Count independent sets by size using inclusion
    # For small alpha (<=3), we can enumerate efficiently
    alpha_0 = 1
    alpha_1 = m

    # alpha_2: pairs of non-conflicting cycles
    alpha_2 = 0
    for i in range(m):
        for j in range(i+1, m):
            if not conflict[i][j]:
                alpha_2 += 1

    # alpha_3: triples of pairwise non-conflicting cycles
    alpha_3 = 0
    for i in range(m):
        for j in range(i+1, m):
            if conflict[i][j]:
                continue
            for k in range(j+1, m):
                if not conflict[i][k] and not conflict[j][k]:
                    alpha_3 += 1

    coeffs = [1, alpha_1, alpha_2, alpha_3]
    # Trim trailing zeros
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

def check_real_roots(coeffs):
    """Check if polynomial with given coefficients has all real roots."""
    if len(coeffs) <= 1:
        return True, float('inf')
    if len(coeffs) == 2:
        return True, float('inf')  # 1 + a1*x always has real root
    if len(coeffs) == 3:
        a1, a2 = coeffs[1], coeffs[2]
        disc = a1*a1 - 4*a2
        return disc >= 0, disc
    if len(coeffs) == 4:
        # Cubic: use numpy to find roots
        poly_coeffs = list(reversed(coeffs))  # highest degree first
        roots = np.roots(poly_coeffs)
        all_real = all(abs(r.imag) < 1e-10 for r in roots)
        return all_real, None
    return None, None

def run_exhaustive(n):
    """Run exhaustive check at given n."""
    num_edges = n * (n - 1) // 2
    total = 1 << num_edges

    min_ratio = float('inf')
    min_ratio_tournament = None
    stats = defaultdict(int)
    failures = 0

    step = max(1, total // 20)

    for mask in range(total):
        bits = [(mask >> i) & 1 for i in range(num_edges)]
        A = tournament_from_bits(n, bits)
        cycles = find_all_odd_cycles(A, n)
        coeffs = compute_independence_poly(cycles)

        degree = len(coeffs) - 1
        stats[degree] += 1

        real, disc = check_real_roots(coeffs)

        if not real:
            failures += 1
            print(f"  FAILURE at mask={mask}: coeffs={coeffs}")

        if degree == 2:
            a1, a2 = coeffs[1], coeffs[2]
            ratio = a1*a1 / (4*a2) if a2 > 0 else float('inf')
            if ratio < min_ratio:
                min_ratio = ratio
                min_ratio_tournament = (mask, coeffs, a1, a2)

        if mask % step == 0:
            print(f"  n={n}: {mask}/{total} ({100*mask/total:.0f}%)", end='\r')

    print(f"\nn={n}: {total} tournaments checked, {failures} failures")
    print(f"  Degree distribution: {dict(stats)}")
    if min_ratio_tournament:
        mask, coeffs, a1, a2 = min_ratio_tournament
        print(f"  Tightest discriminant: a1={a1}, a2={a2}, a1^2/(4*a2) = {min_ratio:.4f}")
        print(f"  (discriminant = {a1*a1 - 4*a2})")
    return failures, min_ratio

def run_sampled(n, num_samples=5000):
    """Run sampled check at given n."""
    import random
    num_edges = n * (n - 1) // 2

    min_ratio = float('inf')
    stats = defaultdict(int)
    failures = 0
    ratios = []

    for trial in range(num_samples):
        bits = [random.randint(0, 1) for _ in range(num_edges)]
        A = tournament_from_bits(n, bits)
        cycles = find_all_odd_cycles(A, n)
        coeffs = compute_independence_poly(cycles)

        degree = len(coeffs) - 1
        stats[degree] += 1

        real, disc = check_real_roots(coeffs)

        if not real:
            failures += 1
            print(f"  FAILURE at trial={trial}: coeffs={coeffs}")

        if degree == 2 and coeffs[2] > 0:
            a1, a2 = coeffs[1], coeffs[2]
            ratio = a1*a1 / (4*a2)
            ratios.append(ratio)
            if ratio < min_ratio:
                min_ratio = ratio
        elif degree == 3 and coeffs[3] > 0:
            a1, a2, a3 = coeffs[1], coeffs[2], coeffs[3]
            # Newton's inequalities for real-rooted degree-3:
            # a1^2 >= 3*a2 and a2^2 >= 3*a1*a3
            n1 = a1*a1 / (3*a2) if a2 > 0 else float('inf')
            n2 = a2*a2 / (3*a1*a3) if a1*a3 > 0 else float('inf')
            ratios.append(min(n1, n2))

    print(f"n={n}: {num_samples} sampled, {failures} failures")
    print(f"  Degree distribution: {dict(stats)}")
    if ratios:
        print(f"  Min ratio (discriminant/Newton): {min(ratios):.4f}")
        print(f"  Mean ratio: {np.mean(ratios):.4f}")
    return failures, min_ratio

if __name__ == '__main__':
    print("="*60)
    print("DISCRIMINANT APPROACH TO REAL-ROOTEDNESS")
    print("="*60)

    # n=5: exhaustive (quick)
    print("\n--- n=5 (exhaustive, 1024 tournaments) ---")
    run_exhaustive(5)

    # n=6: exhaustive (32768 tournaments)
    print("\n--- n=6 (exhaustive, 32768 tournaments) ---")
    run_exhaustive(6)

    # n=7: sampled
    print("\n--- n=7 (sampled, 3000 tournaments) ---")
    run_sampled(7, 3000)

    # n=8: sampled
    print("\n--- n=8 (sampled, 1000 tournaments) ---")
    run_sampled(8, 1000)

    # n=9: sampled (first case where degree can be 3)
    print("\n--- n=9 (sampled, 500 tournaments) ---")
    run_sampled(9, 500)
