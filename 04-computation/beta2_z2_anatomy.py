#!/usr/bin/env python3
"""
beta2_z2_anatomy.py - Understand the structure of Z2 elements in tournaments

Z2 = ker(bd2|_Om2): elements of Om2 whose boundary vanishes.

For a 2-path (a,b,c):
  bd2(a,b,c) = (b,c) - (a,c) + (a,b)

So z = sum alpha_i (a_i,b_i,c_i) in Z2 means:
  For each edge (u,v): sum of alpha_i where (u,v) appears as a face = 0

The faces of (a,b,c) are: (b,c) with sign +1, (a,c) with sign -1, (a,b) with sign +1

And z in Om2 means: for each backward pair (a,c) with c->a,
  sum of alpha over all (a,b,c) = 0

Question: What do Z2 elements look like concretely?
Can we always express them as bd3(w) for some w in Om3?

Key idea: In a tournament, consider a DIRECTED SQUARE
  alpha * (a,b1,c) - alpha * (a,b2,c) where both a->b1->c and a->b2->c
  and c->a (so both paths are "junk" type with bad face (a,c))

The boundary is:
  alpha * [(b1,c) - (a,c) + (a,b1)] - alpha * [(b2,c) - (a,c) + (a,b2)]
  = alpha * [(b1,c) - (b2,c) + (a,b1) - (a,b2)]

For this to be in Z2 we need bd2 of this to vanish... but it's already
in Om2 by construction (the (a,c) terms cancel).

And its bd2 is:
  alpha * [c - b1 - c + a + b1 - a - c + b2 + c - a - b2 + a]
  = alpha * 0 = 0

Wait, that's automatically in Z2! Let me verify more carefully.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A


# ============================================================
# Construct explicit Z2 bases at n=5
# ============================================================
print("=" * 70)
print("Z2 ELEMENT ANATOMY: n=5")
print("=" * 70)

n = 5

# Pick a specific tournament
for bits in [341, 10, 100]:
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    bd2 = build_full_boundary_matrix(a2, a1)
    bd3 = build_full_boundary_matrix(a3, a2)

    # Compute Z2 basis
    bd2_om = bd2 @ om2
    U, S_vals, Vt = np.linalg.svd(bd2_om, full_matrices=True)
    rk = sum(s > 1e-8 for s in S_vals)
    z2_om_basis = Vt[rk:].T  # in Om2 coordinates
    z2_a2_basis = om2 @ z2_om_basis  # in A2 coordinates

    dim_z2 = z2_a2_basis.shape[1] if z2_a2_basis.ndim == 2 else 0

    print(f"\nbits={bits}, scores={scores}, c3={sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n) if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]))}")
    print(f"  |A2|={len(a2)}, dim(Om2)={d_om2}, dim(Z2)={dim_z2}")
    print(f"  |A3|={len(a3)}, dim(Om3)={d_om3}")

    # Show Z2 basis elements
    for col in range(min(dim_z2, 3)):
        z = z2_a2_basis[:, col]
        nonzero = [(i, z[i]) for i in range(len(z)) if abs(z[i]) > 1e-8]
        print(f"\n  Z2 basis element {col}:")
        for i, val in nonzero:
            p = a2[i]
            face_type = "TT" if A[p[0]][p[2]] else "junk"
            print(f"    {val:+.4f} * ({p[0]},{p[1]},{p[2]}) [{face_type}]")

    # For each Z2 element, find a preimage in Om3
    if d_om3 > 0 and dim_z2 > 0:
        bd3_om = bd3 @ om3  # images in A2
        # Express Z2 elements as linear combinations of bd3(om3 basis)
        for col in range(min(dim_z2, 2)):
            z = z2_a2_basis[:, col]
            # Solve bd3_om * x = z (in A2 coordinates)
            # x are coordinates in Om3 basis
            result = np.linalg.lstsq(bd3_om, z, rcond=None)
            x = result[0]
            residual = np.linalg.norm(bd3_om @ x - z)

            if residual < 1e-8:
                # Reconstruct the Om3 element
                w = om3 @ x  # in A3 coordinates
                nonzero_w = [(i, w[i]) for i in range(len(w)) if abs(w[i]) > 1e-8]
                print(f"\n  Preimage of Z2[{col}] in Om3:")
                for i, val in nonzero_w:
                    p = a3[i]
                    dd = "DD" if A[p[0]][p[2]] and A[p[1]][p[3]] else "non-DD"
                    print(f"    {val:+.4f} * ({p[0]},{p[1]},{p[2]},{p[3]}) [{dd}]")


# ============================================================
# KEY QUESTION: Are directed squares ALWAYS in Z2?
# ============================================================
print(f"\n{'='*70}")
print("ARE DIRECTED SQUARES ALWAYS IN Z2?")
print("=" * 70)

# A directed square: (a,b1,c) - (a,b2,c) where c->a (junk pair)
# and a->b1->c, a->b2->c
#
# In Om2? Yes, the junk face (a,c) cancels: coeff of (a,c) is -1+1 = 0
# In Z2 (ker bd2)?
# bd2 of this = [(b1,c) - (a,c) + (a,b1)] - [(b2,c) - (a,c) + (a,b2)]
#              = (b1,c) - (b2,c) + (a,b1) - (a,b2)
# This is automatically 0 in vertex space:
# bd1 of this = [c-b1] - [c-b2] + [b1-a] - [b2-a] = c-b1-c+b2+b1-a-b2+a = 0
# So bd2 of directed square maps to Z1, and bd2 of directed square = 0
# Wait - bd2 maps to A1, and what we computed above is bd1(bd2(square)).
# bd1 o bd2 = 0, so this is always 0. That's trivial.
# But is bd2(square) itself = 0?

# bd2(square) = (b1,c) - (b2,c) + (a,b1) - (a,b2) in A1
# This is a 1-chain. Is it 0? Only if b1=b2 (trivial).
# So directed squares are in Om2 but NOT necessarily in Z2!

# Let me verify:
n = 5
for bits in [341]:
    A = build_adj(n, bits)
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)

    a1_idx = {p: i for i, p in enumerate(a1)}
    a2_idx = {p: i for i, p in enumerate(a2)}

    bd2 = build_full_boundary_matrix(a2, a1)

    # Find a directed square
    for a in range(n):
        for c in range(n):
            if a == c or not A[c][a]:
                continue
            intermediaries = [b for b in range(n) if b != a and b != c and A[a][b] and A[b][c]]
            if len(intermediaries) >= 2:
                b1, b2 = intermediaries[0], intermediaries[1]
                # Square: (a,b1,c) - (a,b2,c)
                i1 = a2_idx.get((a,b1,c))
                i2 = a2_idx.get((a,b2,c))
                if i1 is not None and i2 is not None:
                    sq = np.zeros(len(a2))
                    sq[i1] = 1
                    sq[i2] = -1
                    bd2_sq = bd2 @ sq
                    print(f"\n  Directed square ({a},{b1},{c}) - ({a},{b2},{c}):")
                    print(f"    bd2 = {[(a1[i], bd2_sq[i]) for i in range(len(a1)) if abs(bd2_sq[i]) > 1e-8]}")
                    is_z2 = np.linalg.norm(bd2_sq) < 1e-8
                    print(f"    In Z2? {is_z2}")
                    break
        else:
            continue
        break


# So directed squares are in Om2 but have NONZERO boundary!
# They are NOT in Z2.
# Z2 requires: bd2(z) = 0, which means the boundary edges cancel.

# What IS in Z2 then?
# Z2 elements are combinations of 2-paths whose boundary
# (as 1-chains) vanishes. These are "2-dimensional cycles".

print(f"\n{'='*70}")
print("WHAT GENERATES Z2?")
print("=" * 70)

# Z2 = Om2 intersect ker(bd2)
# From our data: dim(Z2) is about half of dim(Om2)

# Let's look at the "DD part" of Z2:
# DD paths (a,b,c) with a->c have bd2(a,b,c) = (b,c) - (a,c) + (a,b)
# All faces allowed, so DD paths are in Om2.
# But bd2(DD path) = (b,c) - (a,c) + (a,b) which is nonzero.

# So individual DD paths are NOT in Z2 either!
# Z2 requires COMBINATIONS of paths whose bd2 cancels.

# Example Z2 element: take 4 vertices {a,b,c,d} forming a 4-path a->b->c->d
# with a->c (skip edge) and b->d (skip edge)
# Then: (a,b,c) + (b,c,d) - (a,c,d) + (a,b,d)
# bd2 = [(b,c)-(a,c)+(a,b)] + [(c,d)-(b,d)+(b,c)] - [(c,d)-(a,d)+(a,c)] + [(b,d)-(a,d)+(a,b)]
#      = (b,c) - (a,c) + (a,b) + (c,d) - (b,d) + (b,c) - (c,d) + (a,d) - (a,c) + (b,d) - (a,d) + (a,b)
#      = 2(a,b) + 2(b,c) - 2(a,c)
# That's NOT zero.

# So the Z2 structure is more complex. Let me look at the actual computed basis.

n = 5
bits = 341  # Regular tournament (Paley)
A = build_adj(n, bits)
a1 = enumerate_allowed_paths(A, n, 1)
a2 = enumerate_allowed_paths(A, n, 2)

om2 = compute_omega_basis(A, n, 2, a2, a1)
d_om2 = om2.shape[1] if om2.ndim == 2 else 0
bd2 = build_full_boundary_matrix(a2, a1)

bd2_om = bd2 @ om2
U, S_vals, Vt = np.linalg.svd(bd2_om, full_matrices=True)
rk = sum(s > 1e-8 for s in S_vals)
z2_om_basis = Vt[rk:].T
z2_a2_basis = om2 @ z2_om_basis

print(f"\nRegular tournament bits=341:")
print(f"  dim(Om2)={d_om2}, rk(bd2)={rk}, dim(Z2)={z2_a2_basis.shape[1]}")

# Show ALL Z2 basis elements
for col in range(z2_a2_basis.shape[1]):
    z = z2_a2_basis[:, col]
    nonzero = [(i, z[i]) for i in range(len(z)) if abs(z[i]) > 1e-8]
    print(f"\n  Z2[{col}]: ({len(nonzero)} nonzero paths)")
    for i, val in nonzero:
        p = a2[i]
        face_type = "TT" if A[p[0]][p[2]] else "junk"
        print(f"    {val:+.6f} * ({p[0]},{p[1]},{p[2]}) [{face_type}]")

    # Which vertex subsets are involved?
    verts = set()
    for i, val in nonzero:
        for v in a2[i]:
            verts.add(v)
    print(f"    Vertices: {sorted(verts)}")


print("\n\nDone.")
