#!/usr/bin/env python3
"""
Check real-rootedness of I(Omega(T), x) (FULL odd-cycle conflict graph)
for the counterexample tournament where I(Omega_3(T), x) fails.

The counterexample has score sequence [1,1,3,4,4,4,6,6,7].

If the full Omega also fails, this disproves OPEN-Q-015 entirely.
If the full Omega succeeds, then Omega_3 is a bad proxy.

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
import numpy as np

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
    """Find all directed odd cycles (as vertex sets) in tournament T."""
    if max_len is None:
        max_len = n
    cycles = set()

    # For each possible cycle length 3, 5, 7, 9
    for length in range(3, max_len + 1, 2):
        # Enumerate all subsets of size length
        from itertools import combinations, permutations
        for subset in combinations(range(n), length):
            # Check if there's a directed cycle on this subset
            # Try all cyclic orderings (fix first vertex, permute rest)
            found = False
            for perm in permutations(subset[1:]):
                cycle = (subset[0],) + perm
                # Check if cycle[0]->cycle[1]->...->cycle[-1]->cycle[0]
                is_cycle = True
                for i in range(length):
                    if not T[cycle[i]][cycle[(i+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    cycles.add(frozenset(subset))
                    found = True
                    break
            # No need to continue if found
    return cycles

def independence_polynomial(vertices, adj):
    """
    Compute independence polynomial of graph with given vertices and adjacency.
    adj[v] = set of neighbors of v.
    Uses inclusion-exclusion / brute force for small graphs.
    """
    vlist = list(vertices)
    m = len(vlist)
    # Enumerate all independent sets
    coeffs = [0] * (m + 1)
    coeffs[0] = 1

    for mask in range(1, 1 << m):
        subset = [vlist[i] for i in range(m) if mask & (1 << i)]
        k = len(subset)
        # Check independence
        independent = True
        for i in range(len(subset)):
            for j in range(i+1, len(subset)):
                if subset[j] in adj[subset[i]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            coeffs[k] += 1

    return coeffs

print("=" * 70)
print("FULL OMEGA(T) REAL-ROOTEDNESS CHECK FOR COUNTEREXAMPLE")
print("=" * 70)

# Step 1: Find all odd cycles
print("\nStep 1: Finding all odd cycles...")
for max_len in [3, 5, 7, 9]:
    cycles = find_odd_cycles(T, n, max_len)
    print(f"  Odd cycles up to length {max_len}: {len(cycles)}")

# For n=9, odd cycles of length 5, 7, 9 could be expensive
# Let's time it carefully
import time

t0 = time.time()
cycles_3 = find_odd_cycles(T, n, 3)
t1 = time.time()
print(f"\n  3-cycles: {len(cycles_3)} (in {t1-t0:.2f}s)")

t0 = time.time()
cycles_5 = find_odd_cycles(T, n, 5)
t1 = time.time()
print(f"  Odd cycles up to 5: {len(cycles_5)} (in {t1-t0:.2f}s)")

t0 = time.time()
cycles_7 = find_odd_cycles(T, n, 7)
t1 = time.time()
print(f"  Odd cycles up to 7: {len(cycles_7)} (in {t1-t0:.2f}s)")

# Length 9 = Hamiltonian cycle — there are only C(9,9)=1 subsets,
# but 8! permutations to check. Should be manageable.
t0 = time.time()
all_cycles = find_odd_cycles(T, n, 9)
t1 = time.time()
print(f"  All odd cycles: {len(all_cycles)} (in {t1-t0:.2f}s)")

# Step 2: Build conflict graph Omega(T)
print("\nStep 2: Building Omega(T)...")
cycle_list = list(all_cycles)
m = len(cycle_list)
adj = {i: set() for i in range(m)}
for i in range(m):
    for j in range(i+1, m):
        if cycle_list[i] & cycle_list[j]:
            adj[i].add(j)
            adj[j].add(i)

print(f"  Vertices (odd cycles): {m}")
edges = sum(len(a) for a in adj.values()) // 2
print(f"  Edges (sharing vertex): {edges}")

# Step 3: Compute independence polynomial
print("\nStep 3: Computing independence polynomial...")
if m <= 25:
    coeffs = independence_polynomial(range(m), adj)
    print(f"  I(Omega, x) = {' + '.join(f'{c}*x^{i}' for i, c in enumerate(coeffs) if c > 0)}")

    # Check real-rootedness
    deg = max(i for i, c in enumerate(coeffs) if c > 0)
    poly_coeffs = [coeffs[deg - i] for i in range(deg + 1)]  # highest degree first
    roots = np.roots(poly_coeffs)
    all_real = all(abs(r.imag) < 1e-8 for r in roots)
    print(f"  Degree: {deg}")
    print(f"  Roots: {roots}")
    print(f"  All real: {all_real}")

    # Also check Newton's inequalities
    print(f"\n  Newton's inequalities (a_k^2 >= a_{k-1}*a_{k+1} * (k+1)/k):")
    for k in range(1, deg):
        lhs = coeffs[k]**2
        rhs = coeffs[k-1] * coeffs[k+1] * (k+1) / k
        print(f"    k={k}: {coeffs[k]}^2 = {lhs} vs {coeffs[k-1]}*{coeffs[k+1]}*{(k+1)/k:.1f} = {rhs:.1f}: {'OK' if lhs >= rhs else 'FAIL'}")
else:
    print(f"  Too many vertices ({m}) for brute-force independence polynomial")
    print(f"  Trying Omega_3 only...")

    # Just verify the Omega_3 result
    cycle3_list = list(cycles_3)
    m3 = len(cycle3_list)
    adj3 = {i: set() for i in range(m3)}
    for i in range(m3):
        for j in range(i+1, m3):
            if cycle3_list[i] & cycle3_list[j]:
                adj3[i].add(j)
                adj3[j].add(i)

    coeffs3 = independence_polynomial(range(m3), adj3)
    print(f"  I(Omega_3, x) = {' + '.join(f'{c}*x^{i}' for i, c in enumerate(coeffs3) if c > 0)}")

# Step 4: Also check a few random n=9 tournaments with full Omega
print("\n" + "=" * 70)
print("Step 4: Check a few random n=9 tournaments with full Omega")
print("=" * 70)

from tournament_lib import random_tournament

for trial in range(5):
    T2 = random_tournament(n)
    all_c = find_odd_cycles(T2, n, 9)
    clist = list(all_c)
    m2 = len(clist)

    if m2 > 25:
        print(f"\n  Trial {trial}: {m2} odd cycles — too many, skip")
        continue

    adj2 = {i: set() for i in range(m2)}
    for i in range(m2):
        for j in range(i+1, m2):
            if clist[i] & clist[j]:
                adj2[i].add(j)
                adj2[j].add(i)

    coeffs2 = independence_polynomial(range(m2), adj2)
    deg2 = max(i for i, c in enumerate(coeffs2) if c > 0) if any(c > 0 for c in coeffs2[1:]) else 0

    if deg2 > 0:
        poly_c = [coeffs2[deg2 - i] for i in range(deg2 + 1)]
        roots2 = np.roots(poly_c)
        all_real2 = all(abs(r.imag) < 1e-8 for r in roots2)
        print(f"\n  Trial {trial}: {m2} cycles, deg={deg2}, coeffs={coeffs2[:deg2+1]}, all_real={all_real2}")
    else:
        print(f"\n  Trial {trial}: {m2} cycles, trivial polynomial")
