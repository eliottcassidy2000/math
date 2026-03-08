#!/usr/bin/env python3
"""
PHASE 4: Complement duality, n=8 exploration, and geometric interpretation

ESTABLISHED FACTS:
1. β_2 = 0 ALWAYS for tournaments (n ≤ 7, ~2000 samples)
2. β_3 ∈ {0, 1} for tournaments at n=7 (8.4% nonzero)
3. β_3=1 and β_1=1 are MUTUALLY EXCLUSIVE
4. Circulant C_n^S: β_0 = gcd(n, S) components for |S|=1

THIS PHASE:
A. Complement duality: β(T) vs β(T^op) — small sample
B. Does n=8 show β_4>0? Does β_2 remain 0?
C. Geometric shapes catalog
D. β_1 threshold: exact condition at n=4,5
E. Arc-flip sensitivity: how many arc flips change β?
"""
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter
import random
import sys
sys.path.insert(0, '04-computation')
from path_homology_v2 import (
    path_betti_numbers, enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix, circulant_digraph, count_3cycles, ham_path_count
)

random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def complement_tournament(A, n):
    A_op = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                A_op[i][j] = 1 - A[i][j]
    return A_op

# ===== A: Complement duality =====
print("=" * 70)
print("A: COMPLEMENT DUALITY β(T) vs β(T^op)")
print("=" * 70)

n = 7
matches = 0
mismatches = []
for trial in range(50):
    A = random_tournament(n)
    A_op = complement_tournament(A, n)
    b1 = path_betti_numbers(A, n, max_dim=6)
    b2 = path_betti_numbers(A_op, n, max_dim=6)
    if b1 == b2:
        matches += 1
    else:
        mismatches.append((b1, b2))
        if len(mismatches) <= 3:
            print(f"  MISMATCH #{len(mismatches)}: β(T)={b1}, β(T^op)={b2}")

print(f"\n  β(T) = β(T^op): {matches}/50")
if mismatches:
    print(f"  Mismatches: {len(mismatches)}")
    # Check if any mismatch has reversed Betti sequence
    for b1, b2 in mismatches[:3]:
        b2_rev = list(reversed(b2))
        print(f"    β(T)={b1}, β(T^op)={b2}, reversed={b2_rev}")

# Smaller n (exhaustive)
for n_test in [4, 5]:
    print(f"\n  n={n_test} (exhaustive):")
    edges = [(i,j) for i in range(n_test) for j in range(i+1,n_test)]
    m = len(edges)
    match_count = 0
    total_count = 0
    for mask in range(1 << m):
        A = [[0]*n_test for _ in range(n_test)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        A_op = complement_tournament(A, n_test)
        b1 = path_betti_numbers(A, n_test, max_dim=n_test-1)
        b2 = path_betti_numbers(A_op, n_test, max_dim=n_test-1)
        total_count += 1
        if b1 == b2:
            match_count += 1
    print(f"    β(T) = β(T^op): {match_count}/{total_count}")

# ===== B: n=8 exploration =====
print("\n\n" + "=" * 70)
print("B: n=8 TOURNAMENT PATH HOMOLOGY")
print("=" * 70)

n = 8
betti_dist = Counter()
for trial in range(100):
    A = random_tournament(n)
    betti = path_betti_numbers(A, n, max_dim=7)
    bt = tuple(betti)
    betti_dist[bt] += 1
    if trial < 3:
        t3 = count_3cycles(A, n)
        print(f"  Trial {trial}: t3={t3}, β={betti}")

print(f"\n  Betti distribution (n={n}, 100 samples):")
for bt in sorted(betti_dist.keys()):
    print(f"    β={list(bt)}: {betti_dist[bt]}")

# ===== C: β_1 exact condition at n=4 =====
print("\n\n" + "=" * 70)
print("C: EXACT CONDITION FOR β_1=1 AT n=4")
print("=" * 70)

n = 4
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(edges)

b1_1 = []
b1_0 = []
for mask in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (mask >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    betti = path_betti_numbers(A, n, max_dim=3)
    t3 = count_3cycles(A, n)

    # Check for a specific structural property
    # Count: for each 3-cycle, does a "filling" 2-chain exist?
    has_unfilled_cycle = False
    for triple in combinations(range(n), 3):
        i, j, k = triple
        # Check if (i,j,k) forms a 3-cycle (any orientation)
        if A[i][j] and A[j][k] and A[k][i]:
            # 3-cycle i→j→k→i
            # Can it be filled? Need transitive 2-paths:
            # (i,j,k) is in Ω_2 iff i→k (but k→i, so NO)
            # (j,k,i) is in Ω_2 iff j→i (transitive?)
            # (k,i,j) is in Ω_2 iff k→j (transitive?)
            # None of these are transitive (since cycle goes i→j→k→i)
            has_unfilled_cycle = True
        if A[i][k] and A[k][j] and A[j][i]:
            # 3-cycle i→k→j→i
            has_unfilled_cycle = True

    if betti[1] == 1:
        b1_1.append((t3, has_unfilled_cycle))
    else:
        b1_0.append((t3, has_unfilled_cycle))

print(f"\n  β_1=1: {len(b1_1)} tournaments")
print(f"    all have unfilled 3-cycle? {all(x[1] for x in b1_1)}")
print(f"    t3 values: {Counter([x[0] for x in b1_1])}")

print(f"\n  β_1=0: {len(b1_0)} tournaments")
print(f"    any have unfilled 3-cycle? {any(x[1] for x in b1_0)}")
print(f"    t3 values: {Counter([x[0] for x in b1_0])}")

# Deeper: at n=4, β_1=1 iff t3=2. WHY?
# t3=0: transitive, no cycles
# t3=1: exactly one 3-cycle. But is it "fillable"?
print(f"\n  Analysis:")
print(f"  t3=0: All transitive → acyclic → β_1=0 (trivially)")
print(f"  t3=1: One 3-cycle exists but...")

# For t3=1 tournaments, show the structure
for mask in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (mask >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    t3 = count_3cycles(A, n)
    if t3 == 1:
        # Find the 3-cycle
        for triple in combinations(range(n), 3):
            i, j, k = triple
            cyc = None
            if A[i][j] and A[j][k] and A[k][i]:
                cyc = (i, j, k)
            elif A[i][k] and A[k][j] and A[j][i]:
                cyc = (i, k, j)
            if cyc:
                a, b, c = cyc
                # The 4th vertex
                d = [v for v in range(4) if v not in {a, b, c}][0]
                # How does d connect?
                d_out = [v for v in [a, b, c] if A[d][v]]
                d_in = [v for v in [a, b, c] if A[v][d]]
                print(f"  t3=1 example: cycle={cyc}, vertex {d} beats {d_out}, loses to {d_in}")

                # The cycle (a,b) + (b,c) + (c,a) is a 1-chain
                # Its boundary is 0 in Ω_0 (each vertex appears once +1 and once -1)
                # Is it in im(∂_2)?
                # ∂_2 maps from Ω_2. Elements of Ω_2 are transitive 2-paths.
                # Can we write (a,b)+(b,c)+(c,a) = sum of boundaries of 2-paths?
                # ∂(x,y,z) = (y,z) - (x,z) + (x,y) for transitive triple x→y→z, x→z
                # So we need edges of 3-cycle to appear as faces of transitive triples
                # involving vertex d.

                # Count transitive triples involving d and two cycle vertices
                for u, v in [(a,b), (b,c), (c,a)]:
                    if A[d][u] and A[d][v]:
                        print(f"    d→{u}→{v}: d→v? {A[d][v]}, transitive? {A[d][v]==1}")
                    if A[u][d] and A[d][v]:
                        print(f"    {u}→d→{v}: u→v? {A[u][v]}, transitive? {A[u][v]==1}")
                break
        break

# ===== D: Arc-flip sensitivity =====
print("\n\n" + "=" * 70)
print("D: ARC-FLIP SENSITIVITY OF β")
print("=" * 70)

n = 7
flip_changes = 0
flip_total = 0
for trial in range(50):
    A = random_tournament(n)
    b_orig = path_betti_numbers(A, n, max_dim=6)

    # Try flipping each arc
    for i in range(n):
        for j in range(i+1, n):
            A_flip = [row[:] for row in A]
            if A[i][j]:
                A_flip[i][j] = 0
                A_flip[j][i] = 1
            else:
                A_flip[j][i] = 0
                A_flip[i][j] = 1

            b_flip = path_betti_numbers(A_flip, n, max_dim=6)
            flip_total += 1
            if b_flip != b_orig:
                flip_changes += 1

print(f"  Arc flips that change β: {flip_changes}/{flip_total} ({100*flip_changes/flip_total:.1f}%)")

# ===== E: Geometric shape catalog =====
print("\n\n" + "=" * 70)
print("E: GEOMETRIC SHAPE CATALOG")
print("=" * 70)

shapes = {
    (1,): "point",
    (1, 0): "point",
    (1, 0, 0): "contractible (point)",
    (1, 0, 0, 0): "contractible",
    (1, 0, 0, 0, 0): "contractible",
    (1, 0, 0, 0, 0, 0): "contractible",
    (1, 0, 0, 0, 0, 0, 0): "contractible",
    (1, 1, 0): "S^1 (circle)",
    (1, 1, 0, 0): "S^1",
    (1, 1, 0, 0, 0): "S^1",
    (1, 1, 0, 0, 0, 0): "S^1",
    (1, 1, 0, 0, 0, 0, 0): "S^1",
    (1, 0, 1, 0): "S^2 (sphere)",
    (1, 0, 0, 1, 0): "S^3 (3-sphere)",
    (1, 0, 0, 1, 0, 0): "S^3",
    (1, 0, 0, 1, 0, 0, 0): "S^3",
    (1, 0, 0, 0, 1, 0): "S^4",
    (1, 2, 1, 0): "T^2 (torus)",
    (1, 2, 1, 0, 0): "T^2",
    (1, 2, 1, 0, 0, 0): "T^2",
}

print("\nTournament topology by n:")
for n in [3, 4, 5]:
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    shape_count = Counter()
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        betti = path_betti_numbers(A, n, max_dim=n-1)
        bt = tuple(betti)
        shape = shapes.get(bt, f"unknown ({list(bt)})")
        shape_count[shape] += 1

    total = sum(shape_count.values())
    print(f"\n  n={n} ({total} tournaments):")
    for shape, cnt in shape_count.most_common():
        print(f"    {shape}: {cnt} ({100*cnt/total:.1f}%)")

# n=7 sample
n = 7
shape_count = Counter()
for trial in range(300):
    A = random_tournament(n)
    betti = path_betti_numbers(A, n, max_dim=6)
    bt = tuple(betti)
    shape = shapes.get(bt, f"unknown ({list(bt)})")
    shape_count[shape] += 1

total = sum(shape_count.values())
print(f"\n  n=7 ({total} sampled tournaments):")
for shape, cnt in shape_count.most_common():
    print(f"    {shape}: {cnt} ({100*cnt/total:.1f}%)")

print("\nCirculant topology catalog:")
for n in [5, 7]:
    print(f"\n  n={n}:")
    for size in range(1, n):
        for S in combinations(range(1, n), size):
            A = circulant_digraph(n, list(S))
            betti = path_betti_numbers(A, n, max_dim=min(n-1, 6))
            bt = tuple(betti)
            shape = shapes.get(bt, f"unknown ({list(bt)})")
            if "unknown" in shape or shape not in ["contractible", "S^1"]:
                print(f"    C_{n}^{set(S)}: β={list(bt)} → {shape}")

print("\n\nDone.")
