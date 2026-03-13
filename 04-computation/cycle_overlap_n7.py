#!/usr/bin/env python3
"""
cycle_overlap_n7.py — opus-2026-03-13-S67j

INVESTIGATING: Why is H non-monotone in c5 for regular n=7 tournaments?

The OCF: H(T) = I(Ω(T), 2) where Ω(T) is the odd-cycle intersection graph.
Two directed cycles are ADJACENT in Ω iff they share at least one vertex.

For regular tournaments at n=7:
- c3 = 14 (constant), c5 = {28, 36, 42}, c7 = {2, 3}
- Total odd cycle count varies: c3+c5+c7 = {44, 52, 59}
- But H = {175, 171, 189} — not monotone in total cycle count!

This means the OVERLAP STRUCTURE matters:
- More cycles but more overlaps → fewer independent sets → lower H
- The independence number α(Ω) and chromatic number χ(Ω) matter

COMPUTING: For one representative from each class:
1. All directed 3-cycles and 5-cycles
2. The overlap graph Ω
3. I(Ω, 2) directly
4. The structure of Ω (clique number, independence number)
"""

import numpy as np
from itertools import permutations, combinations
import math

def build_paley(p):
    QR = set()
    for k in range(1, p):
        QR.add((k * k) % p)
    A = np.zeros((p, p), dtype=np.int8)
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in QR:
                A[i][j] = 1
    return A

def find_directed_cycles(A, length):
    """Find all directed cycles of given length as vertex sets."""
    n = A.shape[0]
    cycles = []
    for perm in permutations(range(n), length):
        valid = True
        for idx in range(length):
            if A[perm[idx]][perm[(idx+1) % length]] != 1:
                valid = False
                break
        if valid:
            # Canonical form: smallest starting vertex, smallest rotation
            canonical = min(perm[i:] + perm[:i] for i in range(length))
            cycles.append(tuple(canonical))

    # Remove duplicates
    unique_cycles = list(set(cycles))
    return unique_cycles

def independence_polynomial_at_2(adj_matrix):
    """Compute I(G, 2) = sum over independent sets S of 2^|S|."""
    n = adj_matrix.shape[0]
    total = 0
    for mask in range(1 << n):
        # Check if mask is independent set
        vertices = [i for i in range(n) if (mask >> i) & 1]
        independent = True
        for i in range(len(vertices)):
            for j in range(i+1, len(vertices)):
                if adj_matrix[vertices[i]][vertices[j]] == 1:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            total += 2 ** len(vertices)
    return total

n = 7

# =====================================================================
# PART 1: Build representatives of each class
# =====================================================================
print("=" * 70)
print("CYCLE OVERLAP ANALYSIS FOR REGULAR n=7 TOURNAMENTS")
print("=" * 70)

# Paley P_7 (H=189)
A_paley = build_paley(7)

# We need representatives of the other two classes.
# From our sampling, we know they have c5=28 (H=175) and c5=36 (H=171).
# Let's find them by random search.

import random
random.seed(42)

edges = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(edges)

def ham_path_count_dp(A):
    n = A.shape[0]
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0 or not (mask & (1 << v)):
                continue
            for u in range(n):
                if not (mask & (1 << u)) and A[v][u] == 1:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return int(np.sum(dp[(1 << n) - 1]))

reps = {189: A_paley.copy()}
for trial in range(200000):
    bits = random.randint(0, 2**m - 1)
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(edges):
        if (bits >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1

    scores = A.sum(axis=1)
    if not np.all(scores == 3):
        continue

    h = ham_path_count_dp(A)
    if h not in reps:
        reps[h] = A.copy()

    if len(reps) == 3:
        break

print(f"  Found representatives for H classes: {sorted(reps.keys())}")

# =====================================================================
# PART 2: Compute cycle structure and overlap graph for each class
# =====================================================================
for h_val in sorted(reps.keys()):
    A = reps[h_val]

    print(f"\n{'='*60}")
    print(f"  H = {h_val}")
    print(f"{'='*60}")

    # Find all directed 3-cycles
    cycles_3 = find_directed_cycles(A, 3)
    # Find all directed 5-cycles
    cycles_5 = find_directed_cycles(A, 5)
    # Find all directed 7-cycles (Hamiltonian)
    cycles_7 = find_directed_cycles(A, 7)

    print(f"  Directed cycles: c3={len(cycles_3)}, c5={len(cycles_5)}, c7={len(cycles_7)}")

    all_cycles = cycles_3 + cycles_5 + cycles_7
    nc = len(all_cycles)
    print(f"  Total cycles: {nc}")

    # Build overlap graph Ω
    # Two cycles overlap iff they share at least one vertex
    omega = np.zeros((nc, nc), dtype=int)
    for i in range(nc):
        for j in range(i+1, nc):
            si = set(all_cycles[i])
            sj = set(all_cycles[j])
            if si & sj:  # share a vertex
                omega[i][j] = 1
                omega[j][i] = 1

    # Statistics of Ω
    edge_count = np.sum(omega) // 2
    density = edge_count / (nc * (nc-1) / 2) if nc > 1 else 0

    print(f"\n  Overlap graph Ω: {nc} vertices, {edge_count} edges")
    print(f"  Density: {density:.4f}")

    # Degree sequence of Ω
    degrees = omega.sum(axis=1)
    print(f"  Degree range: [{min(degrees)}, {max(degrees)}]")
    print(f"  Mean degree: {np.mean(degrees):.2f}")

    # Independence polynomial I(Ω, 2)
    if nc <= 25:
        I_omega_2 = independence_polynomial_at_2(omega)
        print(f"\n  I(Ω, 2) = {I_omega_2}")
        print(f"  H(T) should equal I(Ω, 2) = {I_omega_2}")
        print(f"  Actual H(T) = {h_val}")
        print(f"  Match? {I_omega_2 == h_val}")
    else:
        print(f"\n  Too many cycles ({nc}) for exact I(Ω, 2) computation")
        # Use alpha_1 = number of cycles, alpha_2 = number of disjoint pairs
        disjoint_pairs = 0
        for i in range(nc):
            for j in range(i+1, nc):
                if omega[i][j] == 0:
                    disjoint_pairs += 1
        print(f"  alpha_1 (total cycles) = {nc}")
        print(f"  alpha_2 (disjoint pairs) = {disjoint_pairs}")
        print(f"  alpha_3 estimate from independent triples...")

        # Approximate: I(Ω, 2) ≈ 1 + 2*nc + 4*disjoint_pairs + ...
        approx = 1 + 2*nc + 4*disjoint_pairs
        print(f"  I(Ω, 2) ≈ 1 + 2*{nc} + 4*{disjoint_pairs} + ... = {approx}+higher")

    # Clique number (max clique in Ω = max overlap set)
    # For small graphs, use greedy
    max_clique = 1
    for v in range(nc):
        neighbors = [u for u in range(nc) if omega[v][u] == 1]
        for u in neighbors:
            common = [w for w in neighbors if omega[u][w] == 1]
            if len(common) + 2 > max_clique:
                max_clique = max(max_clique, len(common) + 2)

    # Independence number (max independent set)
    if nc <= 25:
        max_indep = 0
        for mask in range(1 << nc):
            vertices = [i for i in range(nc) if (mask >> i) & 1]
            independent = True
            for i in range(len(vertices)):
                for j in range(i+1, len(vertices)):
                    if omega[vertices[i]][vertices[j]] == 1:
                        independent = False
                        break
                if not independent:
                    break
            if independent:
                max_indep = max(max_indep, len(vertices))

        print(f"\n  Clique number ω(Ω) = {max_clique}")
        print(f"  Independence number α(Ω) = {max_indep}")
    else:
        print(f"\n  Clique number ω(Ω) ≥ {max_clique}")

    # Block structure: how do 3-cycles overlap 5-cycles?
    c3_count = len(cycles_3)
    c5_count = len(cycles_5)
    c7_count = len(cycles_7)

    # Cross-type overlaps
    overlap_33 = sum(1 for i in range(c3_count) for j in range(i+1, c3_count)
                     if omega[i][j] == 1)
    overlap_55 = sum(1 for i in range(c3_count, c3_count+c5_count)
                     for j in range(i+1, c3_count+c5_count) if omega[i][j] == 1)
    overlap_35 = sum(1 for i in range(c3_count)
                     for j in range(c3_count, c3_count+c5_count) if omega[i][j] == 1)
    overlap_37 = sum(1 for i in range(c3_count)
                     for j in range(c3_count+c5_count, nc) if omega[i][j] == 1)
    overlap_57 = sum(1 for i in range(c3_count, c3_count+c5_count)
                     for j in range(c3_count+c5_count, nc) if omega[i][j] == 1)

    print(f"\n  Overlap structure:")
    print(f"    3-3 overlaps: {overlap_33}/{c3_count*(c3_count-1)//2}")
    print(f"    5-5 overlaps: {overlap_55}/{c5_count*(c5_count-1)//2 if c5_count > 1 else 0}")
    print(f"    3-5 overlaps: {overlap_35}/{c3_count * c5_count}")
    print(f"    3-7 overlaps: {overlap_37}/{c3_count * c7_count}")
    print(f"    5-7 overlaps: {overlap_57}/{c5_count * c7_count}")

print("\n\nDONE — cycle_overlap_n7.py complete")
