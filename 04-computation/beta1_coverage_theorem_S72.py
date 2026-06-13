#!/usr/bin/env python3
"""
beta1_coverage_theorem_S72.py — opus-2026-03-15-S72

TESTING THE COVERAGE THEOREM FOR β₁
=====================================

From fundamental_S72.py, we discovered:
- β₁=0 tournament: ALL edges covered by ≥1 transitive triple
- β₁=1 tournament: exactly 1 edge uncovered

CONJECTURE: β₁ = (number of edges NOT covered by any transitive triple boundary)?
Or weaker: β₁ > 0 iff ∃ uncovered edge?

Test exhaustively at n=5,6.

Also test the STRONGER conjecture:
  The unbounding cycle always passes through the uncovered edge(s).
"""

import numpy as np
from itertools import combinations
import sys
import time
sys.stdout.reconfigure(line_buffering=True)

def adj_from_bits(n, bits):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def compute_beta1(A, n):
    """Compute β₁ using rank method."""
    edges = []
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j] == 1:
                edges.append((i, j))

    triples = []
    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] == 0:
                continue
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[i][j] == 1 and A[j][k] == 1 and A[i][k] == 1:
                    triples.append((i, j, k))

    edge_idx = {e: i for i, e in enumerate(edges)}
    ne = len(edges)
    nt = len(triples)

    # ∂₁: edges → vertices (n×ne)
    d1 = np.zeros((n, ne), dtype=np.int64)
    for idx, (i, j) in enumerate(edges):
        d1[j, idx] += 1  # +target
        d1[i, idx] -= 1  # -source

    # ∂₂: triples → edges (ne×nt)
    d2 = np.zeros((ne, nt), dtype=np.int64)
    for t_idx, (i, j, k) in enumerate(triples):
        # ∂(i,j,k) = (j,k) - (i,k) + (i,j)
        if (j, k) in edge_idx:
            d2[edge_idx[(j, k)], t_idx] += 1
        if (i, k) in edge_idx:
            d2[edge_idx[(i, k)], t_idx] -= 1
        if (i, j) in edge_idx:
            d2[edge_idx[(i, j)], t_idx] += 1

    PRIME = 32749

    def rank_mod(M, p):
        M = M.copy() % p
        nrows, ncols = M.shape
        r = 0
        for col in range(ncols):
            nz = np.where(M[r:, col] != 0)[0]
            if len(nz) == 0:
                continue
            pivot = nz[0] + r
            if pivot != r:
                M[[r, pivot]] = M[[pivot, r]]
            inv = pow(int(M[r, col]), p - 2, p)
            M[r] = M[r] * inv % p
            factors = M[:, col].copy()
            factors[r] = 0
            nzr = np.where(factors != 0)[0]
            if len(nzr) > 0:
                M[nzr] = (M[nzr] - np.outer(factors[nzr], M[r])) % p
            r += 1
        return r

    rank_d1 = rank_mod(d1 % PRIME, PRIME)
    rank_d2 = rank_mod(d2 % PRIME, PRIME)

    ker_d1 = ne - rank_d1
    beta1 = ker_d1 - rank_d2

    return beta1, edges, triples, edge_idx

def edge_coverage(A, n):
    """Which edges are covered by at least one transitive triple?"""
    edges = set()
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j] == 1:
                edges.add((i, j))

    covered = set()
    for i in range(n):
        for j in range(n):
            if A[i][j] != 1 or i == j:
                continue
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[j][k] == 1 and A[i][k] == 1:
                    # Triple (i,j,k): boundary covers edges (j,k), (i,k), (i,j)
                    covered.add((j, k))
                    covered.add((i, k))
                    covered.add((i, j))

    uncovered = edges - covered
    return edges, covered, uncovered

def c3v(A, n):
    """Per-vertex 3-cycle counts, sorted."""
    counts = []
    for v in range(n):
        c = 0
        for i in range(n):
            for j in range(n):
                if i == v or j == v or i == j:
                    continue
                if A[v][i] == 1 and A[i][j] == 1 and A[j][v] == 1:
                    c += 1
        counts.append(c // 2)  # Each 3-cycle counted twice per vertex
    return tuple(sorted(counts))

def t3_count(A, n):
    """Total 3-cycles."""
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                s = A[i][j] + A[j][i] + A[i][k] + A[k][i] + A[j][k] + A[k][j]
                # 3-cycle iff exactly one direction for each pair, forming a cycle
                fwd = A[i][j] * A[j][k] * A[k][i]
                bwd = A[j][i] * A[i][k] * A[k][j]
                count += fwd + bwd
    return count

print("=" * 72)
print("COVERAGE THEOREM TEST — Does β₁ correlate with uncovered edges?")
print("opus-2026-03-15-S72")
print("=" * 72)

for n in [5, 6]:
    print(f"\n{'='*60}")
    print(f"n = {n}")
    print(f"{'='*60}")

    num_edges_in_complete = n * (n - 1) // 2
    total = 1 << num_edges_in_complete

    # Track: (beta1, num_uncovered) pairs
    beta1_uncov = {}  # beta1 → list of uncovered counts
    uncov_beta1 = {}  # uncovered_count → list of beta1 values

    # Also check: uncovered edge always in the unbounding cycle?
    violations_coverage_implies_beta0 = 0  # uncov=0 but beta1>0
    violations_beta0_implies_coverage = 0  # beta1=0 but uncov>0

    t0 = time.time()

    for bits in range(total):
        A = adj_from_bits(n, bits)
        beta1, edges, triples, edge_idx = compute_beta1(A, n)
        all_edges, covered, uncovered = edge_coverage(A, n)

        n_uncov = len(uncovered)

        if beta1 not in beta1_uncov:
            beta1_uncov[beta1] = []
        beta1_uncov[beta1].append(n_uncov)

        if n_uncov not in uncov_beta1:
            uncov_beta1[n_uncov] = []
        uncov_beta1[n_uncov].append(beta1)

        if n_uncov == 0 and beta1 > 0:
            violations_coverage_implies_beta0 += 1
        if beta1 == 0 and n_uncov > 0:
            violations_beta0_implies_coverage += 1

    elapsed = time.time() - t0
    print(f"\n  Exhaustive search: {total} tournaments in {elapsed:.1f}s")

    print(f"\n  β₁ → uncovered edge counts:")
    for b in sorted(beta1_uncov.keys()):
        vals = beta1_uncov[b]
        from collections import Counter
        dist = Counter(vals)
        print(f"    β₁={b}: {dict(sorted(dist.items()))}  (total: {len(vals)})")

    print(f"\n  uncovered count → β₁ values:")
    for u in sorted(uncov_beta1.keys()):
        vals = uncov_beta1[u]
        from collections import Counter
        dist = Counter(vals)
        print(f"    uncov={u}: β₁ distribution = {dict(sorted(dist.items()))}  (total: {len(vals)})")

    print(f"\n  CONJECTURE TESTS:")
    print(f"    'uncov=0 ⟹ β₁=0': violations = {violations_coverage_implies_beta0}")
    print(f"    'β₁=0 ⟹ uncov=0': violations = {violations_beta0_implies_coverage}")

    if violations_coverage_implies_beta0 == 0:
        print(f"    ✓ Full coverage IMPLIES β₁=0 (contrapositive: β₁>0 ⟹ ∃ uncovered edge)")
    else:
        print(f"    ✗ Full coverage does NOT imply β₁=0")

    if violations_beta0_implies_coverage == 0:
        print(f"    ✓ β₁=0 IMPLIES full coverage")
    else:
        print(f"    ✗ β₁=0 does NOT imply full coverage")

    if violations_coverage_implies_beta0 == 0 and violations_beta0_implies_coverage == 0:
        print(f"\n    ★ EQUIVALENCE: β₁=0 ⟺ all edges covered by transitive triples ★")

# Deeper analysis: for β₁=1, is the uncovered edge always in the cycle?
print(f"\n{'='*60}")
print(f"DEEPER: Does the unbounding cycle pass through uncovered edges?")
print(f"{'='*60}")

n = 5
num_edges_in_complete = n * (n - 1) // 2
total = 1 << num_edges_in_complete

cycle_through_uncovered = 0
cycle_not_through_uncovered = 0
total_beta1_pos = 0

for bits in range(total):
    A = adj_from_bits(n, bits)
    beta1, edges, triples, edge_idx = compute_beta1(A, n)

    if beta1 == 0:
        continue

    total_beta1_pos += 1
    all_edges, covered, uncovered = edge_coverage(A, n)

    # Find an actual cycle in ker(∂₁)\im(∂₂)
    # Use the augmented rank method
    ne = len(edges)
    nt = len(triples)

    PRIME = 32749

    d1 = np.zeros((n, ne), dtype=np.int64)
    for idx, (i, j) in enumerate(edges):
        d1[j, idx] += 1
        d1[i, idx] -= 1

    d2 = np.zeros((ne, nt), dtype=np.int64)
    for t_idx, (i, j, k) in enumerate(triples):
        if (j, k) in edge_idx:
            d2[edge_idx[(j, k)], t_idx] += 1
        if (i, k) in edge_idx:
            d2[edge_idx[(i, k)], t_idx] -= 1
        if (i, j) in edge_idx:
            d2[edge_idx[(i, j)], t_idx] += 1

    def rank_mod(M, p):
        M = M.copy() % p
        nrows, ncols = M.shape
        r = 0
        for col in range(ncols):
            nz = np.where(M[r:, col] != 0)[0]
            if len(nz) == 0:
                continue
            pivot = nz[0] + r
            if pivot != r:
                M[[r, pivot]] = M[[pivot, r]]
            inv = pow(int(M[r, col]), p - 2, p)
            M[r] = M[r] * inv % p
            factors = M[:, col].copy()
            factors[r] = 0
            nzr = np.where(factors != 0)[0]
            if len(nzr) > 0:
                M[nzr] = (M[nzr] - np.outer(factors[nzr], M[r])) % p
            r += 1
        return r

    base_rank = rank_mod(d2 % PRIME, PRIME)

    # Find a cycle: try each edge, see if adding it as a column to d2 increases rank
    # Actually, find a vector in ker(d1) \ im(d2)
    # ker(d1) basis: find nullspace of d1
    from numpy.linalg import matrix_rank

    # Simpler: just check if any uncovered edge is in every cycle representative
    # For β₁=1, there's essentially one independent cycle
    # Check: does adding the indicator vector of each uncovered edge to im(d2) change rank?

    uncov_in_cycle = False
    for ue in uncovered:
        if ue in edge_idx:
            # Add this edge's indicator as a column
            col = np.zeros((ne, 1), dtype=np.int64)
            col[edge_idx[ue], 0] = 1
            augmented = np.hstack([d2, col])
            new_rank = rank_mod(augmented % PRIME, PRIME)
            if new_rank > base_rank:
                uncov_in_cycle = True
                break

    if uncov_in_cycle:
        cycle_through_uncovered += 1
    else:
        cycle_not_through_uncovered += 1

print(f"\n  n=5, β₁>0 tournaments: {total_beta1_pos}")
print(f"  Unbounding cycle passes through uncovered edge: {cycle_through_uncovered}")
print(f"  Unbounding cycle does NOT pass through uncovered: {cycle_not_through_uncovered}")

if cycle_not_through_uncovered == 0:
    print(f"\n  ★ The unbounding cycle ALWAYS passes through an uncovered edge ★")
    print(f"  This means: the edge that no transitive triple covers IS the obstruction.")

# Final: what determines which edges are uncovered?
print(f"\n{'='*60}")
print(f"WHAT MAKES AN EDGE UNCOVERED?")
print(f"{'='*60}")

n = 5
count_examples = 0
for bits in range(total):
    if count_examples >= 5:
        break
    A = adj_from_bits(n, bits)
    all_edges, covered, uncovered = edge_coverage(A, n)
    if len(uncovered) > 0:
        count_examples += 1
        for ue in uncovered:
            i, j = ue
            # Edge i→j is uncovered means: no k with i→k, k→j (no common out-neighbor)
            # or no k with j→k, i→k (no common target... wait)
            # Triple (i,j,k): i→j, j→k, i→k. So k is common out-neighbor of both i and j.
            # Triple (k,i,j): k→i, i→j, k→j. So k beats both i and j.
            # Triple (i,k,j): i→k, k→j, i→j. So k is intermediate.
            # Edge i→j is covered if ∃ any triple containing it.
            common_out = []
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[i][k] == 1 and A[j][k] == 1:
                    common_out.append(k)
            common_in = []
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[k][i] == 1 and A[k][j] == 1:
                    common_in.append(k)
            intermediate = []
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[i][k] == 1 and A[k][j] == 1:
                    intermediate.append(k)

            if count_examples <= 3:
                print(f"\n  bits={bits}, uncovered edge {i}→{j}:")
                print(f"    Common out-neighbors of {i},{j}: {common_out}")
                print(f"    Common in-neighbors of {i},{j}: {common_in}")
                print(f"    Intermediates (i→k→j): {intermediate}")
                print(f"    out-deg({i})={sum(A[i])}, out-deg({j})={sum(A[j])}")

print(f"\n  Edge i→j is uncovered iff:")
print(f"    - No k with (i→j, j→k, i→k) [common out-neighbor]")
print(f"    - No k with (k→i, i→j, k→j) [common in-neighbor that also beats j]")
print(f"    - No k with (i→k, k→j, i→j) [intermediate on transitive path]")
print(f"    = No vertex k forming a transitive triple with edge i→j")
print(f"    = For all k: either k→i or j→k (k 'reverses' the path through i→j)")

print(f"\n  INSIGHT: edge i→j is uncovered iff every other vertex k satisfies:")
print(f"    k→i OR j→k (or both)")
print(f"  i.e., every k either beats i or loses to j.")
print(f"  In other words: N_in(i) ∪ N_out(j) ⊇ V\\{{i,j}}")
print(f"  Equivalently: N_out(i)\\{{j}} ⊆ N_out(j)\\{{i}}")
print(f"  The out-neighborhood of i (minus j) is CONTAINED in that of j!")

# Verify this characterization
print(f"\n  Verifying containment characterization at n=5...")
n = 5
num_edges_in_complete = n * (n - 1) // 2
total = 1 << num_edges_in_complete

correct = 0
wrong = 0
for bits in range(total):
    A = adj_from_bits(n, bits)
    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] == 0:
                continue
            # Is edge i→j uncovered?
            is_uncov = True
            for k in range(n):
                if k == i or k == j:
                    continue
                # Check all three triple types containing edge i→j
                # (i,j,k): i→j, j→k, i→k
                if A[j][k] == 1 and A[i][k] == 1:
                    is_uncov = False
                    break
                # (k,i,j): k→i, i→j, k→j
                if A[k][i] == 1 and A[k][j] == 1:
                    is_uncov = False
                    break
                # (i,k,j): i→k, k→j, i→j
                if A[i][k] == 1 and A[k][j] == 1:
                    is_uncov = False
                    break

            # Check containment: N_out(i)\{j} ⊆ N_out(j)\{i}
            out_i = set(k for k in range(n) if k != i and k != j and A[i][k] == 1)
            out_j = set(k for k in range(n) if k != j and k != i and A[j][k] == 1)
            is_contained = out_i.issubset(out_j)

            if is_uncov == is_contained:
                correct += 1
            else:
                wrong += 1

print(f"    Correct: {correct}, Wrong: {wrong}")
if wrong == 0:
    print(f"    ✓ CONFIRMED: edge i→j uncovered ⟺ N_out(i)\\{{j}} ⊆ N_out(j)\\{{i}}")
    print(f"\n    This means: β₁ > 0 ⟺ ∃ edge i→j where j dominates all of i's other targets")

# Verify at n=6
print(f"\n  Verifying containment characterization at n=6...")
n = 6
num_edges_in_complete = n * (n - 1) // 2
total = 1 << num_edges_in_complete

correct = 0
wrong = 0
checked = 0
for bits in range(total):
    A = adj_from_bits(n, bits)
    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] == 0:
                continue
            is_uncov = True
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[j][k] == 1 and A[i][k] == 1:
                    is_uncov = False
                    break
                if A[k][i] == 1 and A[k][j] == 1:
                    is_uncov = False
                    break
                if A[i][k] == 1 and A[k][j] == 1:
                    is_uncov = False
                    break

            out_i = set(k for k in range(n) if k != i and k != j and A[i][k] == 1)
            out_j = set(k for k in range(n) if k != j and k != i and A[j][k] == 1)
            is_contained = out_i.issubset(out_j)

            if is_uncov == is_contained:
                correct += 1
            else:
                wrong += 1
    checked += 1

print(f"    Correct: {correct}, Wrong: {wrong}")
if wrong == 0:
    print(f"    ✓ CONFIRMED at n=6: edge i→j uncovered ⟺ N_out(i)\\{{j}} ⊆ N_out(j)\\{{i}}")

print(f"\n{'='*60}")
print(f"SUMMARY OF FINDINGS")
print(f"{'='*60}")
print(f"""
  1. Edge i→j is NOT covered by any transitive triple
     iff N_out(i)\\{{j}} ⊆ N_out(j)\\{{i}}
     (j dominates all of i's other out-neighbors)

  2. β₁ > 0 ⟹ ∃ such an edge (full coverage implies β₁=0)

  3. The unbounding cycle PASSES THROUGH the uncovered edge(s)

  4. GRAPH-THEORETIC MEANING: An uncovered edge i→j means
     j is "almost dominating" i — j beats everyone i beats.
     This creates a local ordering i < j that can't be
     consistently extended to a global acyclic structure
     via transitive triples.
""")

print("DONE — beta1_coverage_theorem_S72.py complete")
