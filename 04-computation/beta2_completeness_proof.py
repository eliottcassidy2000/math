#!/usr/bin/env python3
"""
beta2_completeness_proof.py — Attempt to prove KEY LEMMA using tournament completeness

KEY LEMMA: For tournament T, vertex v: H₂(T,T\v) = max(0, β₁(T\v) - β₁(T))

The quotient complex R = Ω(T)/Ω(T\v):
  R₀ = ℝ (just vertex v)
  R₁ = n-1 dimensional (arcs involving v)
  R₂ = "Ω₂ elements involving v" (mod T\v)
  R₃ = "Ω₃ elements involving v" (mod T\v)

Structure of R₁:
  Arcs involving v: (v,w) for w ∈ out(v), (w,v) for w ∈ in(v)
  dim R₁ = n-1
  ∂₁^R: (v,w) ↦ [w] - [v] = -[v]  (since [w] = 0 in R₀)
  ∂₁^R: (w,v) ↦ [v] - [w] = +[v]
  So ∂₁^R maps every arc to ±[v]. rk(∂₁^R) = 1.
  ker(∂₁^R) has dim n-2.

Structure of R₂:
  "v-triples": Ω₂ elements involving v, modulo Ω₂(T\v).
  A 2-path (a,b,c) involves v iff v ∈ {a,b,c}.
  Types:
    (v,b,c): v is first vertex, v→b, b→c
    (a,v,c): v is middle vertex, a→v, v→c
    (a,b,v): v is last vertex, a→b, b→v

  For each type, the path is allowed iff the arcs exist (they do in tournament).
  But it's in Ω₂ iff the boundary is in Ω₁.

  For (v,b,c): ∂₂ = (b,c) - (v,c) + (v,b). Faces (v,c) and (v,b) involve v.
  For (a,v,c): ∂₂ = (v,c) - (a,c) + (a,v). Faces (v,c), (a,v) involve v; (a,c) might or might not.
  For (a,b,v): ∂₂ = (b,v) - (a,v) + (a,b). Faces (b,v), (a,v) involve v.

  In A₂, ALL such triples are allowed (tournament = complete). So we're looking at Ω₂ ∩ "involves v".

IDEA: Can we compute dim R₂ and rk(∂₂^R) explicitly?

∂₂^R: R₂ → R₁. The boundary of a v-triple in R is:
  Only the faces involving v survive in R₁ (faces not involving v are in Ω₁(T\v) = 0 in R₁).

For (v,b,c): ∂₂^R = -(v,c) + (v,b)  [the (b,c) face is in T\v, goes to 0]
For (a,v,c): ∂₂^R = (v,c) + (a,v)   [the (a,c) face is in T\v, goes to 0]
For (a,b,v): ∂₂^R = (b,v) - (a,v)   [the (a,b) face is in T\v, goes to 0]

Wait, these need to be in Ω₂ first. Let me check which A₂ elements involving v are in Ω₂.

Author: opus-2026-03-08-S45
"""
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved


def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0


n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

print("=" * 70)
print("QUOTIENT BOUNDARY MAP ANALYSIS")
print("=" * 70)

# For one specific tournament, analyze R₂ → R₁ in detail
bits = 0  # Transitive tournament
A = [[0]*n for _ in range(n)]
for idx, (i,j) in enumerate(pairs):
    if (bits >> idx) & 1: A[i][j] = 1
    else: A[j][i] = 1

print(f"\nT#{bits}:")
for i in range(n):
    out = [j for j in range(n) if A[i][j]]
    print(f"  {i} → {out}")

v = 0
print(f"\nVertex v = {v}")
print(f"  out(v) = {[j for j in range(n) if A[v][j]]}")
print(f"  in(v) = {[j for j in range(n) if A[j][v]]}")

# Compute Ω₂(T) and classify by v-involvement
ap0 = enumerate_allowed_paths(A, n, 0)
ap1 = enumerate_allowed_paths(A, n, 1)
ap2 = enumerate_allowed_paths(A, n, 2)
ap3 = enumerate_allowed_paths(A, n, 3)

om1 = compute_omega_basis(A, n, 1, ap1, ap0)
om2 = compute_omega_basis(A, n, 2, ap2, ap1)
om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0, 0))

d2 = dim_om(om2)
d3 = dim_om(om3)

print(f"\n  dim Ω₂(T) = {d2}")
print(f"  dim Ω₃(T) = {d3}")

# Classify each Ω₂ basis element
ap2_list = [tuple(x) for x in ap2]
for col in range(d2):
    vec = om2[:, col]
    terms = [(ap2_list[i], round(vec[i], 4)) for i in range(len(vec)) if abs(vec[i]) > 1e-8]
    involves_v = any(v in path for path, _ in terms)
    v_pos = set()
    for path, _ in terms:
        if v in path:
            v_pos.add(path.index(v))
    print(f"  Ω₂ e{col}: involves_v={involves_v}, v_positions={v_pos if involves_v else '-'}")
    for path, coeff in terms:
        print(f"    {coeff:+.4f} * {path}")

# Now do the BOUNDARY MAP analysis in the quotient
print(f"\n{'='*70}")
print("R₂ → R₁ MAP IN QUOTIENT COMPLEX")
print("=" * 70)

# R₁ basis: arcs involving v = {(v,w) : v→w} ∪ {(w,v) : w→v}
r1_basis = []
for w in range(n):
    if w == v:
        continue
    if A[v][w]:
        r1_basis.append(('→', v, w))  # (v,w)
    else:
        r1_basis.append(('←', w, v))  # (w,v)

print(f"R₁ basis ({len(r1_basis)} elements):")
for arc in r1_basis:
    if arc[0] == '→':
        print(f"  ({arc[1]},{arc[2]})")
    else:
        print(f"  ({arc[1]},{arc[2]})")

# For each Ω₂ element involving v, compute its image in R₁
# The image is: ∂₂(element) projected to R₁ (dropping faces not involving v)
bd2 = build_full_boundary_matrix(ap2, ap1)
ap1_list = [tuple(x) for x in ap1]
ap1_idx = {p: i for i, p in enumerate(ap1_list)}

# R₁ arcs as indices in A₁
r1_arcs = []
for w in range(n):
    if w == v:
        continue
    if A[v][w]:
        r1_arcs.append((v, w))
    else:
        r1_arcs.append((w, v))

r1_arc_idx = {arc: i for i, arc in enumerate(r1_arcs)}

print(f"\n∂₂^R computation:")
# For each Ω₂ basis element, compute boundary and project to R₁
boundary_in_R1 = []
for col in range(d2):
    vec_a2 = om2[:, col]  # in A₂ coords
    bd_vec = bd2 @ vec_a2  # in A₁ coords

    # Project to R₁: keep only arcs involving v
    r1_proj = np.zeros(len(r1_arcs))
    for i, arc in enumerate(r1_arcs):
        if arc in ap1_idx:
            r1_proj[i] = bd_vec[ap1_idx[arc]]

    # Check if this Ω₂ element involves v (in the Ω basis, not necessarily directly)
    involves = any(abs(r1_proj[i]) > 1e-8 for i in range(len(r1_arcs)))
    if involves:
        terms = [(r1_arcs[i], round(r1_proj[i], 4)) for i in range(len(r1_arcs)) if abs(r1_proj[i]) > 1e-8]
        boundary_in_R1.append((col, r1_proj))
        print(f"  ∂₂^R(e{col}) = {terms}")

print(f"\n  {len(boundary_in_R1)} Ω₂ elements with nonzero R₁ projection")

# Build the ∂₂^R matrix
if boundary_in_R1:
    bd2_R = np.column_stack([proj for _, proj in boundary_in_R1])
    rk_bd2_R = np.linalg.matrix_rank(bd2_R, tol=1e-8)
    print(f"  rk(∂₂^R) = {rk_bd2_R}")
    print(f"  dim R₂ = {len(boundary_in_R1)}")
    print(f"  ker(∂₂^R) = {len(boundary_in_R1) - rk_bd2_R}")
else:
    print(f"  rk(∂₂^R) = 0")

# Now check multiple tournaments
print(f"\n{'='*70}")
print("SYSTEMATIC: ∂₂^R RANK FOR ALL n=5 (T,v)")
print("=" * 70)

results = Counter()
for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    om1 = compute_omega_basis(A, n, 1, ap1, enumerate_allowed_paths(A, n, 0))
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0, 0))
    d2 = dim_om(om2)

    for v in range(n):
        if d2 == 0:
            results[('d_out', sum(A[v]), 'dim_R2', 0, 'rk', 0)] += 1
            continue

        bd2 = build_full_boundary_matrix(ap2, ap1)
        ap1_list = [tuple(x) for x in ap1]
        ap1_idx = {p: i for i, p in enumerate(ap1_list)}

        # R₁ arcs
        r1_arcs = []
        for w in range(n):
            if w == v: continue
            if A[v][w]: r1_arcs.append((v, w))
            else: r1_arcs.append((w, v))

        # Project each Ω₂ basis to R₁
        bd2_R_cols = []
        for col in range(d2):
            vec_a2 = om2[:, col]
            bd_vec = bd2 @ vec_a2
            r1_proj = np.zeros(len(r1_arcs))
            for i, arc in enumerate(r1_arcs):
                if arc in ap1_idx:
                    r1_proj[i] = bd_vec[ap1_idx[arc]]
            if np.max(np.abs(r1_proj)) > 1e-8:
                bd2_R_cols.append(r1_proj)

        dim_R2 = len(bd2_R_cols)
        if dim_R2 > 0:
            bd2_R = np.column_stack(bd2_R_cols)
            rk_bd2_R = np.linalg.matrix_rank(bd2_R, tol=1e-8)
        else:
            rk_bd2_R = 0

        d_out = sum(A[v])
        ker_bd2_R = dim_R2 - rk_bd2_R
        results[(d_out, dim_R2, rk_bd2_R, ker_bd2_R)] += 1

print(f"\n(d_out, dim_R₂, rk_∂₂^R, ker_∂₂^R): count")
for key in sorted(results.keys()):
    print(f"  {key}: {results[key]}")

# KEY: ker(∂₂^R) is what we need to kill with im(∂₃^R)
# H₂(R) = ker(∂₂^R) - rk(∂₃^R restricted to R₃)
# For β₂=0: H₂(R) = max(0, β₁(T\v) - β₁(T))

print("\nDone.")
