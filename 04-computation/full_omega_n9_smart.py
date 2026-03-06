#!/usr/bin/env python3
"""
Compute I(Omega(T), x) for the n=9 counterexample using a smarter algorithm.

The counterexample has 62 odd cycles. Brute force is 2^62, impossible.
Use the deletion-contraction recurrence for independence polynomials.

I(G, x) = I(G-v, x) + x * I(G-N[v], x)

where G-v removes vertex v, and G-N[v] removes v and all its neighbors.
This has worst-case 2^n but works well for dense graphs (many neighbors => fast reduction).

For our graph with 62 vertices but very dense (1882 edges), many vertices
have high degree, so N[v] is large and the recursion shrinks fast.

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
import numpy as np
from itertools import combinations, permutations
from functools import lru_cache

# The counterexample tournament
T = [
    [0, 1, 0, 1, 0, 0, 1, 1, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0],
    [1, 1, 0, 0, 1, 1, 1, 1, 0],
    [0, 0, 1, 0, 0, 1, 0, 1, 0],
    [1, 1, 0, 1, 0, 0, 0, 1, 0],
    [1, 1, 0, 0, 1, 0, 1, 1, 1],
    [0, 1, 0, 1, 1, 0, 0, 1, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 1, 1, 1, 0, 1, 1, 0],
]
n = 9

def find_odd_cycles(T, n, max_len=None):
    if max_len is None:
        max_len = n
    cycles = set()
    for length in range(3, max_len + 1, 2):
        for subset in combinations(range(n), length):
            for perm in permutations(subset[1:]):
                cycle = (subset[0],) + perm
                is_cycle = True
                for i in range(length):
                    if not T[cycle[i]][cycle[(i+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    cycles.add(frozenset(subset))
                    break
    return cycles

def independence_poly_dc(adj_matrix):
    """
    Compute independence polynomial using deletion-contraction.
    adj_matrix: dict mapping vertex -> frozenset of neighbors
    Returns: list of coefficients [a0, a1, a2, ...]
    """
    # Use memoization on the vertex set
    memo = {}

    def solve(vertices):
        key = vertices  # frozenset
        if key in memo:
            return memo[key]

        if not vertices:
            memo[key] = [1]
            return [1]

        # Pick vertex with highest degree (greedy for speed)
        v = max(vertices, key=lambda u: len(adj_matrix[u] & vertices))
        remaining = vertices - {v}

        # I(G, x) = I(G-v, x) + x * I(G-N[v], x)
        p1 = solve(remaining)
        closed_nbrs = (adj_matrix[v] & vertices) | {v}
        p2 = solve(vertices - closed_nbrs)

        # Add polynomials: p1 + x * p2
        result = list(p1) + [0] * max(0, len(p2) + 1 - len(p1))
        for i in range(len(result)):
            if i < len(p1):
                pass  # already there
            # now add x * p2
        # Redo properly
        maxlen = max(len(p1), len(p2) + 1)
        result = [0] * maxlen
        for i in range(len(p1)):
            result[i] += p1[i]
        for i in range(len(p2)):
            result[i + 1] += p2[i]

        memo[key] = result
        return result

    all_verts = frozenset(adj_matrix.keys())
    return solve(all_verts)

print("=" * 70)
print("FULL OMEGA(T) INDEPENDENCE POLYNOMIAL — COUNTEREXAMPLE")
print("=" * 70)

# Find all odd cycles
print("\nFinding odd cycles...")
all_cycles = find_odd_cycles(T, n, 9)
cycle_list = list(all_cycles)
m = len(cycle_list)
print(f"  {m} odd cycles found")

# Build adjacency
print("Building adjacency...")
adj = {}
for i in range(m):
    adj[i] = frozenset(j for j in range(m) if j != i and cycle_list[i] & cycle_list[j])

degrees = [len(adj[i]) for i in range(m)]
print(f"  Degree distribution: min={min(degrees)}, max={max(degrees)}, mean={sum(degrees)/len(degrees):.1f}")

# Compute independence polynomial
print("Computing independence polynomial (deletion-contraction)...")
import time
t0 = time.time()
coeffs = independence_poly_dc(adj)
t1 = time.time()
print(f"  Time: {t1-t0:.2f}s")
print(f"  I(Omega, x) coefficients: {coeffs}")

deg = len(coeffs) - 1
while deg > 0 and coeffs[deg] == 0:
    deg -= 1

print(f"  Degree: {deg}")
poly_str = " + ".join(f"{coeffs[i]}*x^{i}" for i in range(deg+1) if coeffs[i] > 0)
print(f"  I(Omega, x) = {poly_str}")

# Check real-rootedness
if deg > 0:
    poly_c = [coeffs[deg - i] for i in range(deg + 1)]
    roots = np.roots(poly_c)
    all_real = all(abs(r.imag) < 1e-8 for r in roots)
    print(f"\n  Roots: {roots}")
    print(f"  All real: {all_real}")

    # Check Newton's inequalities
    print(f"\n  Newton's inequalities:")
    for k in range(1, deg):
        lhs = coeffs[k]**2
        rhs_exact = coeffs[k-1] * coeffs[k+1] * (k+1) / k
        ok = lhs >= rhs_exact - 0.001
        print(f"    k={k}: a_{k}^2={lhs} vs a_{k-1}*a_{k+1}*(k+1)/k = {rhs_exact:.1f}: {'OK' if ok else 'FAIL'}")

# Also evaluate at x=2
val_at_2 = sum(coeffs[i] * 2**i for i in range(len(coeffs)))
print(f"\n  I(Omega, 2) = {val_at_2}")

# Compare with H(T) from Hamiltonian path counting
print("\n  Computing H(T) by direct enumeration...")
from tournament_lib import count_hamiltonian_paths_by_parity
# Actually let's just count
def hamiltonian_paths(T, n):
    count = 0
    from itertools import permutations
    for perm in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if not T[perm[i]][perm[i+1]]:
                valid = False
                break
        if valid:
            count += 1
    return count

h = hamiltonian_paths(T, n)
print(f"  Total Hamiltonian paths: {h}")
print(f"  H(T) should equal I(Omega, 2) = {val_at_2}")
print(f"  Match: {h == val_at_2}")
