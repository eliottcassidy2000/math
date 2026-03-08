#!/usr/bin/env python3
"""
DEEP ANALYSIS: Path homology of tournaments and circulant digraphs

Questions:
1. What exactly determines β_1 for tournaments? Is it t3 ≥ some threshold?
2. Do higher Betti numbers (β_2, β_3) ever appear for tournaments?
3. Connection between β_1 and H(T) (Hamiltonian path count)?
4. The Ω_p subspace structure for tournaments
5. Circulant digraphs: which connection sets produce which topologies?
6. Can tournaments be viewed as "deformations" of circulant digraphs?
"""
import numpy as np
from itertools import permutations, combinations
from collections import defaultdict, Counter
import random
from math import comb

# Import from v2
import sys
sys.path.insert(0, '04-computation')
from path_homology_v2 import (
    path_betti_numbers, enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix, circulant_digraph, count_3cycles, ham_path_count
)

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

random.seed(42)

# ===== SECTION 1: What determines β_1 for tournaments? =====
print("=" * 70)
print("SECTION 1: WHAT DETERMINES β_1 FOR TOURNAMENTS?")
print("=" * 70)

for n in [3, 4, 5]:
    print(f"\n--- n={n} ---")

    data = []
    for A in all_tournaments(n):
        t3 = count_3cycles(A, n)
        H = ham_path_count(A, n)
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        betti = path_betti_numbers(A, n, max_dim=min(n-1, 4))
        data.append((t3, H, scores, betti))

    # Does β_1 depend on t3?
    by_t3 = defaultdict(list)
    for t3, H, scores, betti in data:
        by_t3[t3].append((betti[1], H, scores))

    print(f"  β_1 vs t3:")
    for t3 in sorted(by_t3.keys()):
        b1_vals = [x[0] for x in by_t3[t3]]
        H_vals = [x[1] for x in by_t3[t3]]
        print(f"    t3={t3}: β_1 values = {Counter(b1_vals)}, H range = [{min(H_vals)}, {max(H_vals)}]")

    # Does β_1 depend on score sequence?
    by_scores = defaultdict(list)
    for t3, H, scores, betti in data:
        by_scores[scores].append(betti[1])

    print(f"\n  β_1 vs score sequence:")
    for scores in sorted(by_scores.keys()):
        b1_vals = Counter(by_scores[scores])
        if len(b1_vals) > 1:
            print(f"    scores={scores}: β_1 = {dict(b1_vals)} (MIXED)")
        else:
            b1 = list(b1_vals.keys())[0]
            print(f"    scores={scores}: β_1 = {b1} ({b1_vals[b1]} tournaments)")

# ===== SECTION 2: The critical condition for β_1 =====
print("\n\n" + "=" * 70)
print("SECTION 2: THE CRITICAL CONDITION FOR β_1 ≥ 1")
print("=" * 70)
print("""
β_1 ≥ 1 means there exists a directed cycle in the path complex that
doesn't bound anything. For tournaments, this should relate to the
existence of certain types of directed cycles.

In GLMY theory, a directed 3-cycle v_0→v_1→v_2→v_0 gives a 1-cycle
in the path complex IF the 2-chain (v_0,v_1,v_2) is NOT in Ω_2.

The key: (v_0,v_1,v_2) is in A_2 (allowed 2-path) iff v_0→v_1→v_2.
Its boundary is:
  ∂(v_0,v_1,v_2) = (v_1,v_2) - (v_0,v_2) + (v_0,v_1)

For (v_0,v_1,v_2) to be in Ω_2, we need ∂(v_0,v_1,v_2) ∈ A_1.
This requires (v_0,v_2) to be an allowed 1-path, i.e., v_0→v_2 ∈ E.

So the 2-path (v_0,v_1,v_2) is ∂-invariant iff v_0→v_2 (shortcut edge exists).
This is a CRUCIAL condition for tournaments:
  - The 2-path (v_0,v_1,v_2) is ∂-invariant iff v_0→v_2 (transitive)
  - If v_2→v_0 instead, the 2-path is NOT ∂-invariant

So Ω_2 consists of 2-paths where the "shortcut" edge goes the same direction.
For tournaments, (v_0→v_1→v_2) is in Ω_2 iff v_0→v_2 (transitive triple).
""")

# Verify this understanding
for n in [4]:
    print(f"\nn={n}: Verifying Ω_2 structure")
    A_list = list(all_tournaments(n))

    A = A_list[30]  # pick one
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    allowed_1 = enumerate_allowed_paths(A, n, 1)

    omega_2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
    print(f"  |A_2| = {len(allowed_2)}, dim(Ω_2) = {omega_2.shape[1]}")

    # For each allowed 2-path, check if it's ∂-invariant
    transitive_count = 0
    non_transitive_count = 0
    for path in allowed_2:
        v0, v1, v2 = path
        if A[v0][v2] == 1:  # shortcut exists
            transitive_count += 1
        else:
            non_transitive_count += 1

    print(f"  Transitive 2-paths: {transitive_count}")
    print(f"  Non-transitive 2-paths (v_2→v_0): {non_transitive_count}")
    print(f"  dim(Ω_2) should ≤ {transitive_count}")
    # Actually, Ω_2 is the span of ∂-invariant 2-paths within A_2.
    # But it might be bigger if linear combinations of non-invariant paths
    # cancel their non-allowed faces.

# ===== SECTION 3: Explicit cycle generators for β_1 =====
print("\n\n" + "=" * 70)
print("SECTION 3: EXPLICIT 1-CYCLES GENERATING β_1")
print("=" * 70)

for n in [3, 4]:
    print(f"\nn={n}:")
    for idx, A in enumerate(all_tournaments(n)):
        betti = path_betti_numbers(A, n, max_dim=2)
        if betti[1] == 0:
            continue

        t3 = count_3cycles(A, n)
        print(f"\n  Tournament #{idx} with β_1={betti[1]}, t3={t3}:")
        for i in range(n):
            row = [A[i][j] for j in range(n)]
            print(f"    {i}: → {[j for j in range(n) if A[i][j]]}")

        # Find the 1-cycles
        allowed_1 = enumerate_allowed_paths(A, n, 1)
        allowed_0 = enumerate_allowed_paths(A, n, 0)
        omega_1 = compute_omega_basis(A, n, 1, allowed_1, allowed_0)

        # Ω_1 = A_1 for tournaments (every edge deletion is a vertex, which is allowed)
        print(f"    |A_1| = {len(allowed_1)}, dim(Ω_1) = {omega_1.shape[1]}")

        # Build ∂_1 matrix and find kernel
        bd_1 = build_full_boundary_matrix(allowed_1, allowed_0)
        bd_1_omega = bd_1 @ omega_1

        U, S, Vt = np.linalg.svd(bd_1_omega, full_matrices=True)
        rank = sum(s > 1e-8 for s in S)
        kernel_vecs = Vt[rank:].T

        print(f"    ker(∂_1) dimension: {omega_1.shape[1] - rank}")

        # Show kernel vectors in terms of edges
        for k in range(kernel_vecs.shape[1]):
            vec = omega_1 @ kernel_vecs[:, k]
            nonzero = [(allowed_1[i], round(vec[i], 4)) for i in range(len(allowed_1)) if abs(vec[i]) > 1e-8]
            if len(nonzero) <= 10:
                print(f"    Cycle {k}: {nonzero}")

        if idx > 5:
            break

# ===== SECTION 4: n=6 sampling =====
print("\n\n" + "=" * 70)
print("SECTION 4: TOURNAMENTS n=6,7 (SAMPLING)")
print("=" * 70)

for n in [6, 7]:
    print(f"\nn={n}: Sampling 100 random tournaments")
    betti_dist = defaultdict(int)
    b1_by_t3 = defaultdict(list)

    for trial in range(100):
        A = random_tournament(n)
        t3 = count_3cycles(A, n)
        betti = path_betti_numbers(A, n, max_dim=min(n-1, 5))
        bt = tuple(betti)
        betti_dist[bt] += 1
        b1_by_t3[t3].append(betti[1])

        if trial < 3:
            H = ham_path_count(A, n)
            print(f"  Trial {trial}: t3={t3}, H={H}, β={betti}")

    print(f"\n  Betti distribution:")
    for bt in sorted(betti_dist.keys()):
        print(f"    β={list(bt)}: {betti_dist[bt]}")

    print(f"\n  β_1 vs t3:")
    for t3 in sorted(b1_by_t3.keys())[:10]:
        vals = Counter(b1_by_t3[t3])
        print(f"    t3={t3}: β_1 = {dict(vals)}")

# ===== SECTION 5: Circulant digraphs — higher-dimensional topology =====
print("\n\n" + "=" * 70)
print("SECTION 5: CIRCULANT DIGRAPHS — HIGHER DIMENSIONS")
print("=" * 70)

# C_n^S for various S
print("\nSystematic exploration of C_n^S:")
for n in [5, 6, 7, 8]:
    print(f"\n  n={n}:")
    # Generate all subsets of {1,...,n-1}
    for size in range(1, min(n, 5)):
        for S in combinations(range(1, n), size):
            S_list = list(S)
            A = circulant_digraph(n, S_list)
            betti = path_betti_numbers(A, n, max_dim=min(n-1, 5))
            # Only print interesting ones (non-trivial topology)
            if any(b > 0 for b in betti[1:]):
                print(f"    C_{n}^{set(S_list)}: β = {betti}")

# ===== SECTION 6: Connection between β_1 and OCF =====
print("\n\n" + "=" * 70)
print("SECTION 6: β_1 AND THE OCF")
print("=" * 70)

for n in [4, 5]:
    print(f"\nn={n}:")
    b1_by_H = defaultdict(list)
    for A in all_tournaments(n):
        H = ham_path_count(A, n)
        betti = path_betti_numbers(A, n, max_dim=2)
        b1_by_H[H].append(betti[1])

    print(f"  β_1 vs H:")
    for H in sorted(b1_by_H.keys()):
        vals = Counter(b1_by_H[H])
        print(f"    H={H}: β_1 = {dict(vals)}")

    # Correlation
    all_b1 = []
    all_H = []
    for A in all_tournaments(n):
        H = ham_path_count(A, n)
        betti = path_betti_numbers(A, n, max_dim=2)
        all_b1.append(betti[1])
        all_H.append(H)

    corr = np.corrcoef(all_b1, all_H)[0,1]
    print(f"  Corr(β_1, H) = {corr:.4f}")

# ===== SECTION 7: Dimension of Ω_p for tournaments =====
print("\n\n" + "=" * 70)
print("SECTION 7: DIMENSION OF Ω_p FOR TOURNAMENTS")
print("=" * 70)

for n in [4, 5]:
    print(f"\nn={n}:")

    dims = defaultdict(list)
    for A in all_tournaments(n):
        t3 = count_3cycles(A, n)
        for p in range(n):
            allowed_p = enumerate_allowed_paths(A, n, p)
            allowed_pm1 = enumerate_allowed_paths(A, n, p-1) if p > 0 else []
            omega_basis = compute_omega_basis(A, n, p, allowed_p, allowed_pm1)
            dim_omega = omega_basis.shape[1] if omega_basis.ndim == 2 else 0
            dims[(p, t3)].append((len(allowed_p), dim_omega))

    print(f"  (p, t3) -> (|A_p|, dim Ω_p):")
    for (p, t3) in sorted(dims.keys()):
        ap_vals = [x[0] for x in dims[(p, t3)]]
        om_vals = [x[1] for x in dims[(p, t3)]]
        if len(set(ap_vals)) == 1 and len(set(om_vals)) == 1:
            print(f"    p={p}, t3={t3}: |A_p|={ap_vals[0]}, dim(Ω_p)={om_vals[0]}")
        else:
            print(f"    p={p}, t3={t3}: |A_p| in [{min(ap_vals)},{max(ap_vals)}], dim(Ω_p) in [{min(om_vals)},{max(om_vals)}]")

print("\n\nDone.")
