#!/usr/bin/env python3
"""
β_2 = 0 STRUCTURAL PROOF: THE TOURNAMENT EXCHANGE ARGUMENT

GOAL: Prove ker(∂_2|Ω_2) ⊆ im(∂_3|Ω_3) for all tournaments.

KEY STRUCTURAL INSIGHT (from computational evidence):
- Ω_3 is NOT just the "doubly transitive" 4-paths
- Ω_3 includes LINEAR COMBINATIONS of 3-paths where non-allowed faces cancel
- The full Ω_3 is ALWAYS large enough to fill ker(∂_2|Ω_2)

APPROACH: Instead of showing individual filling, show an EXACT SEQUENCE
using the tournament exchange property.

The tournament exchange: for any {a,b,c,d} with a→b→c→d:
  Either a→c AND b→d (doubly transitive, gives Ω_3 element directly)
  Or NOT (then need linear combinations)

Let me trace through what happens with the Ω_3 chains more carefully.
"""
import numpy as np
from itertools import combinations, permutations
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix
)

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

# ===== Understanding Ω_3 structure =====
print("=" * 70)
print("Ω_3 STRUCTURE: INDIVIDUAL VS CHAIN ELEMENTS")
print("=" * 70)

n = 5
# Pick a specific tournament where ker > im for DT-paths-only
# Tournament #4 had MISMATCH with DT-only approach
A = list(all_tournaments_gen(n))[4]

print("Tournament #4 adjacency:")
for i in range(n):
    for j in range(n):
        if A[i][j]: print(f"  {i}→{j}")

a3 = enumerate_allowed_paths(A, n, 3)
a2 = enumerate_allowed_paths(A, n, 2)
a1 = enumerate_allowed_paths(A, n, 1)

print(f"\n|A_3| = {len(a3)}")
print(f"|A_2| = {len(a2)}")

# Ω_2 = transitive triples
tt = [tuple(p) for p in a2 if A[p[0]][p[2]] == 1]
print(f"|Ω_2| = {len(tt)}")
for t in tt:
    print(f"  {t[0]}→{t[1]}→{t[2]} (with {t[0]}→{t[2]})")

# DT 4-paths (individually ∂-invariant in Ω_3)
dt = [tuple(p) for p in a3 if A[p[0]][p[2]] == 1 and A[p[1]][p[3]] == 1]
print(f"\nDT 4-paths: {len(dt)}")
for p in dt:
    print(f"  {p}")

# Full Ω_3 (using compute_omega_basis)
om3 = compute_omega_basis(A, n, 3, a3, a2)
dim_om3 = om3.shape[1] if om3.ndim == 2 else 0
print(f"\ndim(Ω_3) = {dim_om3}")

# The difference: dim(Ω_3) vs number of DT paths
print(f"\nDT paths only give {len(dt)} generators")
print(f"But dim(Ω_3) = {dim_om3}")
print(f"So Ω_3 has {dim_om3 - len(dt)} extra dimensions from chain cancellation")

if dim_om3 > 0:
    print(f"\nΩ_3 basis vectors:")
    a3_arr = [tuple(p) for p in a3]
    for col in range(dim_om3):
        vec = om3[:, col]
        nonzero = [(a3_arr[j], vec[j]) for j in range(len(vec)) if abs(vec[j]) > 1e-8]
        print(f"  v{col}: ", end="")
        for path, coeff in nonzero:
            print(f" {coeff:+.3f}·{path}", end="")
        print()

        # Check which faces of this chain are NOT in A_2 → they must cancel
        boundary = {}
        for path, coeff in nonzero:
            v0, v1, v2, v3 = path
            faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
            signs = [1, -1, 1, -1]
            for face, sign in zip(faces, signs):
                if face in boundary:
                    boundary[face] += sign * coeff
                else:
                    boundary[face] = sign * coeff

        # Which faces are NOT in A_2?
        a2_set = set(tuple(p) for p in a2)
        for face, val in sorted(boundary.items()):
            if abs(val) < 1e-10: continue
            in_a2 = face in a2_set
            in_tt = face in set(tt)
            status = "A_2" if in_a2 else "NOT-A_2"
            status2 = "Ω_2" if in_tt else "NOT-Ω_2"
            print(f"    ∂ face {face}: coeff={val:+.4f}, {status}, {status2}")

# ===== Critical analysis: What are the non-DT Ω_3 elements? =====
print(f"\n\n{'='*70}")
print("NON-DT Ω_3 ELEMENTS: CANCELLATION STRUCTURE")
print("="*70)

# For each 3-path, check which faces are NOT in A_2
for path_tuple in a3:
    v0, v1, v2, v3 = path_tuple
    faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]

    a2_set_tuples = set(tuple(p) for p in a2)
    non_allowed = []
    for i, (face, sign) in enumerate(zip(faces, [1,-1,1,-1])):
        if face not in a2_set_tuples:
            non_allowed.append((i, face, sign))

    if len(non_allowed) > 0:
        # This path has non-allowed faces
        pass  # We'll summarize below

# Summary: classify 3-paths by their face pattern
print("\nFace pattern classification of 3-paths:")
from collections import Counter
patterns = Counter()
for path_tuple in a3:
    v0, v1, v2, v3 = path_tuple
    faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
    a2_set_tuples = set(tuple(p) for p in a2)

    pattern = tuple(1 if face in a2_set_tuples else 0 for face in faces)
    patterns[pattern] += 1

for pattern, count in sorted(patterns.items()):
    # pattern[i] = 1 means face i is in A_2
    # Face 0: (v1,v2,v3), always in A_2 (v1→v2→v3 are consecutive edges)
    # Face 3: (v0,v1,v2), always in A_2
    # Face 1: (v0,v2,v3), in A_2 iff v0→v2 (edge exists)
    # Face 2: (v0,v1,v3), in A_2 iff v1→v3 (edge exists)
    labels = ["(v1,v2,v3)", "(v0,v2,v3)", "(v0,v1,v3)", "(v0,v1,v2)"]
    missing = [labels[i] for i in range(4) if pattern[i] == 0]
    print(f"  pattern={pattern}: {count} paths, missing: {missing}")

# ===== Key insight =====
print(f"""

KEY INSIGHT:
For a 3-path (v0,v1,v2,v3) with edges v0→v1→v2→v3:
  Face 0 = (v1,v2,v3): ALWAYS in A_2 (edges v1→v2, v2→v3 exist)
  Face 3 = (v0,v1,v2): ALWAYS in A_2 (edges v0→v1, v1→v2 exist)
  Face 1 = (v0,v2,v3): in A_2 iff v0→v2 (edge exists)
  Face 2 = (v0,v1,v3): in A_2 iff v1→v3 (edge exists)

In a tournament, v0→v2 XOR v2→v0, and v1→v3 XOR v3→v1.

So the face pattern is:
  (1, v0→v2?, v1→v3?, 1)

There are 4 possible patterns:
  (1,1,1,1): ALL faces in A_2 (DT-path, automatically in Ω_3)
  (1,1,0,1): Face 2 missing → v3→v1 (NOT v1→v3)
  (1,0,1,1): Face 1 missing → v2→v0 (NOT v0→v2)
  (1,0,0,1): Both faces missing → v2→v0 AND v3→v1

A single path with missing faces is NOT in Ω_3.
But TWO paths can form a chain in Ω_3 if their missing faces cancel!

Cancellation requires: two paths have the SAME missing face type with
OPPOSITE signs. Since faces are step sequences in the circulant case,
this is controlled by the junk matrix.
""")

# ===== Check at n=4 which pairs cancel =====
print("=" * 70)
print("FACE CANCELLATION: WHICH 3-PATH PAIRS FORM Ω_3 CHAINS?")
print("=" * 70)

n = 4
count_examples = 0
for A_idx, A in enumerate(all_tournaments_gen(n)):
    a3 = enumerate_allowed_paths(A, n, 3)
    a2 = enumerate_allowed_paths(A, n, 2)
    a2_set = set(tuple(p) for p in a2)

    # Find 3-paths with missing faces
    paths_with_missing = []
    for path_arr in a3:
        path = tuple(path_arr)
        v0, v1, v2, v3 = path
        faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
        missing = [face for face in faces if face not in a2_set]
        if missing:
            paths_with_missing.append((path, missing))

    if len(paths_with_missing) >= 2 and count_examples < 3:
        count_examples += 1
        print(f"\nTournament #{A_idx}:")
        for path, missing in paths_with_missing:
            print(f"  path {path}: missing faces {missing}")

        # Full Ω_3
        om3 = compute_omega_basis(A, n, 3, a3, a2)
        dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

        # DT paths
        dt_count = sum(1 for p in a3 if A[p[0]][p[2]]==1 and A[p[1]][p[3]]==1)
        print(f"  DT paths: {dt_count}, dim(Ω_3)={dim_om3}")

        if dim_om3 > dt_count:
            print(f"  Extra Ω_3 elements: {dim_om3 - dt_count} (from cancellation)")

print("\nDone.")
