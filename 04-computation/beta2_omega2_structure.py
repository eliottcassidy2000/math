#!/usr/bin/env python3
"""
Ω_2 STRUCTURE: NON-TRANSITIVE-TRIPLE ELEMENTS

We discovered dim(Ω_2) > |transitive triples| for many tournaments.
This means there are Ω_2 chains involving non-transitive 2-paths.

For a non-transitive 2-path (a,b,c) with a→b, b→c, c→a:
∂(a,b,c) = (b,c) - (a,c) + (a,b)
Face (a,c): a→c? NO, c→a. So (a,c) ∉ A_1.

But (c,a) IS in A_1 (since c→a is an edge).
The non-allowed face (a,c) is a 1-path from a to c, but the edge goes c→a.

For a CHAIN of non-TT paths to be in Ω_2, the non-allowed 1-path
contributions must cancel. Two paths (a,b,c) and (a,b',c) with c→a
both produce -(a,c) face, so they ADD rather than cancel.

BUT: what if another path produces +(a,c)?
∂(d,a,c) = (a,c) - (d,c) + (d,a). Face 0 = (a,c) with sign +1.
Path (d,a,c) needs d→a and a→c. But c→a, so a→c is false!
So (d,a,c) ∉ A_2.

What about (a,c,e)? Needs a→c. Same problem.

Hmm, so the non-allowed 1-path (a,c) can only appear with sign -1 from
paths (a,*,c). It can't cancel. Unless...

Wait — (a,c) is an ORDERED pair. In A_1, the allowed 1-paths are (u,v)
where u→v. So (a,c) is allowed iff a→c, and (c,a) is allowed iff c→a.

The non-allowed face of (a,b,c) with c→a is (a,c), an ordered pair
that is NOT an allowed 1-path. For this to cancel, we need another
path whose boundary has the term (a,c) with opposite sign.

But (a,c) appears in:
- ∂(d,a,c) = +(a,c) - ... → needs d→a, a→c. Can't because c→a.
- ∂(a,c,e) = ... + (a,c) → Face 2, sign +1. Needs a→c. Can't.
- ∂(a,b,c) = -(a,c) → Face 1, sign -1. Always gives -(a,c).

So (a,c) can ONLY appear with sign -1. It CAN'T be cancelled!

But then how is dim(Ω_2) > #TT? There must be a BUG.

Let me check compute_omega_basis directly.
"""
import numpy as np
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import enumerate_allowed_paths, compute_omega_basis, boundary_coeffs

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

# ===== Check a specific tournament =====
print("=" * 70)
print("Ω_2 STRUCTURE ANALYSIS")
print("=" * 70)

n = 5
for t_idx, A in enumerate(all_tournaments_gen(n)):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    tt = [tuple(p) for p in a2 if A[p[0]][p[2]] == 1]
    ntt = [tuple(p) for p in a2 if A[p[0]][p[2]] == 0]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0

    if dim_om2 != len(tt):
        print(f"\nT#{t_idx}: dim(Ω_2)={dim_om2}, |TT|={len(tt)}, |non-TT|={len(ntt)}")
        print(f"  Edges: ", end="")
        for i in range(n):
            for j in range(n):
                if A[i][j]: print(f"{i}→{j}", end=" ")
        print()

        a2_tuples = [tuple(p) for p in a2]

        # Show the Ω_2 basis
        for col in range(dim_om2):
            vec = om2[:, col]
            terms = [(a2_tuples[j], vec[j]) for j in range(len(a2_tuples)) if abs(vec[j]) > 1e-8]
            tt_terms = [(t,c) for t,c in terms if t in set(tt)]
            ntt_terms = [(t,c) for t,c in terms if t in set(ntt)]
            print(f"  Ω_2 basis {col}: {len(tt_terms)} TT + {len(ntt_terms)} non-TT terms")
            for triple, coeff in ntt_terms:
                a,b,c = triple
                print(f"    non-TT: {coeff:+.4f}·{triple} (a→c={A[a][c]}, c→a={A[c][a]})")

            # Verify: ∂ of this basis vector should be in A_1
            boundary = {}
            for triple, coeff in terms:
                for sign, face in boundary_coeffs(triple):
                    if face in boundary:
                        boundary[face] += sign * coeff
                    else:
                        boundary[face] = sign * coeff

            # Check which face terms are non-allowed
            a1_set = set(tuple(p) for p in a1)
            non_a1 = [(face, val) for face, val in boundary.items() if abs(val) > 1e-8 and face not in a1_set]
            if non_a1:
                print(f"    *** NON-A_1 boundary terms: {non_a1} ***")
            else:
                print(f"    ∂ ∈ A_1 ✓")

        if t_idx > 20:
            break

# ===== The key: what does compute_omega_basis actually compute? =====
print(f"\n\n{'='*70}")
print("DEBUGGING compute_omega_basis")
print("="*70)

# Let me manually compute Ω_2 for a specific tournament
A = list(all_tournaments_gen(5))[4]
print("Tournament #4:")
for i in range(5):
    for j in range(5):
        if A[i][j]: print(f"  {i}→{j}")

a1 = enumerate_allowed_paths(A, 5, 1)
a2 = enumerate_allowed_paths(A, 5, 2)
a1_set = set(tuple(p) for p in a1)
a2_tuples = [tuple(p) for p in a2]

print(f"\nA_1 (edges): {[tuple(p) for p in a1]}")
print(f"A_2 (2-paths): {a2_tuples}")

# For each 2-path, compute its boundary and find non-A_1 faces
print("\nBoundary analysis:")
for path in a2:
    path_t = tuple(path)
    a,b,c = path_t
    faces = [(b,c), (a,c), (a,b)]
    signs = [1, -1, 1]
    non_a1_faces = []
    for face, sign in zip(faces, signs):
        if face not in a1_set:
            non_a1_faces.append((face, sign))
    is_tt = A[a][c] == 1
    if non_a1_faces:
        print(f"  {path_t}: {'TT' if is_tt else 'non-TT'}, non-A_1 faces: {non_a1_faces}")
    else:
        print(f"  {path_t}: TT, all faces in A_1 ✓")

# Now manually compute Ω_2 = ker(projection to non-A_1 faces)
# This is what compute_omega_basis does
print("\nManual Ω_2 computation:")
non_a1_faces_all = {}
na_count = 0
for j, path in enumerate(a2):
    for sign, face in boundary_coeffs(tuple(path)):
        if len(set(face)) == len(face) and face not in a1_set:
            if face not in non_a1_faces_all:
                non_a1_faces_all[face] = na_count
                na_count += 1

print(f"Non-A_1 face types: {na_count}")
for face, idx in sorted(non_a1_faces_all.items(), key=lambda x: x[1]):
    print(f"  {face}: index {idx}")

# Build projection matrix P (same as compute_omega_basis)
P = np.zeros((na_count, len(a2)))
for j, path in enumerate(a2):
    for sign, face in boundary_coeffs(tuple(path)):
        if face in non_a1_faces_all:
            P[non_a1_faces_all[face], j] += sign

print(f"\nProjection matrix P ({na_count} x {len(a2)}):")
for i in range(na_count):
    face = [f for f, idx in non_a1_faces_all.items() if idx == i][0]
    row = P[i, :]
    terms = [(a2_tuples[j], row[j]) for j in range(len(a2)) if abs(row[j]) > 0]
    print(f"  face {face}: {terms}")

rank_P = np.linalg.matrix_rank(P, tol=1e-8)
print(f"\nrank(P) = {rank_P}, dim(A_2) = {len(a2)}")
print(f"dim(Ω_2) = dim(A_2) - rank(P) = {len(a2) - rank_P}")
print(f"|TT| = {sum(1 for p in a2 if A[p[0]][p[2]]==1)}")

# Is the discrepancy from a face being shared by two non-TT paths?
# E.g., two non-TT paths (a,b,c) and (a,b',c) both have face (a,c),
# but with the SAME sign (-1). So they can't cancel each other.
# But maybe the rank of P is less than #non-A_1-faces because
# some faces are algebraically dependent?

print(f"\nNon-A_1 faces that appear in multiple 2-paths:")
for face, idx in non_a1_faces_all.items():
    row = P[idx, :]
    nonzero = [(a2_tuples[j], row[j]) for j in range(len(a2)) if abs(row[j]) > 0]
    if len(nonzero) > 1:
        print(f"  {face}: {nonzero}")

print("\nDone.")
