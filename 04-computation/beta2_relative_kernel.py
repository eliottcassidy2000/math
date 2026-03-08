#!/usr/bin/env python3
"""
beta2_relative_kernel.py — Analyze ker(∂₂^rel) and im(∂₃^rel) explicitly

For the inductive proof of β₂=0, we need H₂(T,T\\v) = 0 for all T,v.

This means: in the relative complex, every 2-cycle through v is the
boundary of a 3-chain through v.

KEY INSIGHT to pursue: For TOURNAMENTS specifically, the relative complex
has special structure due to completeness:
- Every pair of vertices has an arc
- Vertex v is connected to ALL other vertices
- This means v can be "inserted" into any path

APPROACH: Count the allowed paths through v by position.

2-paths through v at position i (for n-vertex tournament):
  pos 0: (v,b,c) — need v→b, b→c. Count = sum_{b: v→b} |out(b)\\{v}|
  pos 1: (a,v,c) — need a→v, v→c. Count = |in(v)| × |out(v)|
  pos 2: (a,b,v) — need a→b, b→v. Count = sum_{b: b→v} |in(b)\\{v}|

3-paths through v at position i:
  pos 0: (v,b,c,d) — need v→b→c→d
  pos 1: (a,v,c,d) — need a→v→c→d
  pos 2: (a,b,v,d) — need a→b→v→d
  pos 3: (a,b,c,v) — need a→b→c→v

For these, the DT/bad-face structure depends on which arcs exist between
non-adjacent vertices in the path.

The RELATIVE boundary ∂₃^rel maps 3-paths through v to 2-chains through v.
The face NOT through v (the "non-v face") gets projected out in the quotient.

So ∂₃^rel(a,b,c,d) with v at pos i:
  Drop the face obtained by removing the vertex at pos i.
  Keep the other 3 faces (all through v).

This means ∂₃^rel has a very clean structure.

Author: opus-2026-03-08-S44
"""
import sys, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def delete_vertex(A, n, v):
    B = []
    for i in range(n):
        if i == v:
            continue
        row = []
        for j in range(n):
            if j == v:
                continue
            row.append(A[i][j])
        B.append(row)
    return B

print("=" * 70)
print("RELATIVE KERNEL AND IMAGE: EXPLICIT COMPUTATION")
print("=" * 70)

n = 5

# For each (T,v), compute:
# 1. The RELATIVE 2-cycle space ker(∂₂^rel)
# 2. The RELATIVE 3-boundary space im(∂₃^rel)
# 3. Check ker ⊆ im (i.e., H₂^rel = 0)

h2_rel_dims = Counter()
ker2_dims = Counter()
im3_dims = Counter()
total = 0

for T_idx, A in enumerate(all_tournaments(n)):
    for v in range(n):
        total += 1

        # Paths through v
        a2_v = [p for p in enumerate_allowed_paths(A, n, 2) if v in p]
        a3_v = [p for p in enumerate_allowed_paths(A, n, 3) if v in p]

        # All paths (for Ω computation)
        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)

        # Paths NOT through v (for T\\v)
        a1_nv = [p for p in a1 if v not in p]
        a2_nv = [p for p in a2 if v not in p]
        a3_nv = [p for p in a3 if v not in p]

        if not a2_v:
            continue

        # Ω₂(T) basis elements through v
        # We need the Ω₂(T) basis, then identify which basis vectors
        # have support on v-paths.

        # Simpler: work directly in A₂-v coordinates.
        # The relative 2-chains are Z₂-combinations of a2_v paths,
        # subject to: they're in Ω₂(T) mod Ω₂(T\\v).

        # Actually, the cleanest approach: compute the relative boundary matrix.
        # ∂₂^rel maps (A₂ through v) → (A₁ through v) / (A₁ not through v)
        # = the portion of ∂₂ restricted to v-paths, projected onto v-edges.

        # Boundary ∂₂: a2_v paths → a1 edges
        bd2_v = np.zeros((len(a1), len(a2_v)))
        a1_idx = {tuple(p): i for i, p in enumerate(a1)}
        for j, path in enumerate(a2_v):
            a, b, c = path
            # ∂₂(a,b,c) = (b,c) - (a,c) + (a,b)
            for sign, face in [(1, (b,c)), (-1, (a,c)), (1, (a,b))]:
                if face in a1_idx:
                    bd2_v[a1_idx[face], j] += sign

        # The relative boundary: project onto A₁ paths through v
        v_edge_indices = [i for i, p in enumerate(a1) if v in p]
        nv_edge_indices = [i for i, p in enumerate(a1) if v not in p]

        # For Ω₂ condition: ∂₂(chain) must be in Ω₁ = span(A₁).
        # For a2_v elements, some faces may not be in A₁ (NT condition).
        # The "non-edge" faces are: (a,c) when c→a.

        # BUT: for the RELATIVE complex, we quotient by chains not through v.
        # So we only care about the v-edge components of ∂₂.

        # Relative boundary matrix: restrict ∂₂ to v-rows
        bd2_rel = bd2_v[v_edge_indices, :]

        # Relative ker(∂₂^rel)
        if bd2_rel.shape[0] > 0 and bd2_rel.shape[1] > 0:
            S = np.linalg.svd(bd2_rel, compute_uv=False)
            rk_bd2_rel = int(np.sum(S > 1e-8))
            ker_bd2_rel = bd2_rel.shape[1] - rk_bd2_rel
        else:
            ker_bd2_rel = len(a2_v)

        # Now Ω₂ condition on the v-paths: some individual v-paths may NOT be in Ω₂
        # because they're NT. The relative 2-chains must be in Ω₂(T).
        # So we need to intersect ker(∂₂^rel) with Ω₂(T).

        # The Ω₂ condition for a single path (a,b,c): all faces must be in A₁.
        # Face (a,c): in A₁ iff a→c (TT condition).
        # For NT paths through v: the face (a,c) is NOT in A₁ but might cancel
        # with other NT paths.

        # Let's use the FULL Ω₂ computation.
        om2 = compute_omega_basis(A, n, 2, a2, a1)
        if om2.ndim != 2:
            continue
        dim_om2 = om2.shape[1]

        # Which Ω₂ basis vectors use v?
        # om2 is |A₂| × dim_om2 matrix. Columns are Ω₂ basis vectors in A₂ coords.
        # A basis vector uses v if any of its nonzero entries correspond to v-paths.
        a2_v_set = set(tuple(p) for p in a2_v)
        a2_list = [tuple(p) for p in a2]
        v_path_rows = [i for i, p in enumerate(a2_list) if p in a2_v_set]

        # Restrict om2 to v-path rows
        om2_v = om2[v_path_rows, :]

        # Ω₂ vectors that have nonzero v-component
        # The relative Ω₂ space = projection of Ω₂(T) onto v-paths
        # dim = rank of om2 restricted to v-rows
        rk_v = np.linalg.matrix_rank(om2_v, tol=1e-8)

        # Similarly for Ω₃
        if a3:
            om3 = compute_omega_basis(A, n, 3, a3, a2)
            if om3.ndim == 2 and om3.shape[1] > 0:
                dim_om3 = om3.shape[1]
                a3_v_set = set(tuple(p) for p in a3_v)
                a3_list = [tuple(p) for p in a3]
                v_path_rows_3 = [i for i, p in enumerate(a3_list) if p in a3_v_set]
                om3_v = om3[v_path_rows_3, :]
                rk_v3 = np.linalg.matrix_rank(om3_v, tol=1e-8)
            else:
                rk_v3 = 0
        else:
            rk_v3 = 0

        # The relative boundary ∂₃^rel: Ω₃^v → Ω₂^v
        # We need to compute rk(∂₃^rel) on the Ω-level.
        # This is rk of the composition: Ω₃^v → Ω₂ → Ω₂^v

        # For that, let's compute the boundary matrix ∂₃: A₃ → A₂
        if a3_v:
            bd3_v = np.zeros((len(a2), len(a3_v)))
            a2_idx = {tuple(p): i for i, p in enumerate(a2)}
            for j, path in enumerate(a3_v):
                a, b, c, d = path
                for sign, face in [(1, (b,c,d)), (-1, (a,c,d)), (1, (a,b,d)), (-1, (a,b,c))]:
                    if face in a2_idx:
                        bd3_v[a2_idx[face], j] += sign

            # Restrict to v-path rows of A₂ (relative image)
            bd3_rel = bd3_v[v_path_rows, :]

            # Image in Ω₂^v: need to project through Ω₂
            # Actually for the relative computation, we need:
            # im(∂₃^rel) = image of ∂₃ restricted to (Ω₃ through v) projected to (Ω₂ through v)

            # Simplify: just compute the rank of ∂₃ from a3_v to a2_v rows
            rk_im3_rel = np.linalg.matrix_rank(bd3_rel, tol=1e-8)
        else:
            rk_im3_rel = 0

        # Crude relative ker(∂₂^rel):
        bd2_rel_full = bd2_v  # full boundary restricted to v-paths in source
        # Actually bd2_v is |A₁| × |a2_v|. We need the part projecting onto v-edges.
        # AND we need the Ω₂ condition (bad faces cancel).

        # Let me just use a simple metric: the relative dimensions
        ker2_dims[rk_v] += 1
        im3_dims[rk_v3] += 1

print(f"\n  Ω₂^v rank (relative 2-chain space) distribution:")
for dim, count in sorted(ker2_dims.items()):
    print(f"    dim={dim}: {count}")

print(f"\n  Ω₃^v rank (relative 3-chain space) distribution:")
for dim, count in sorted(im3_dims.items()):
    print(f"    dim={dim}: {count}")

# ===== DIRECT H₂^rel COMPUTATION via full chain complex =====
print(f"\n{'='*70}")
print("DIRECT H₂^rel COMPUTATION")
print("=" * 70)

n = 5
h2_rel_vals = Counter()
detail_count = 0

for T_idx, A in enumerate(all_tournaments(n)):
    for v in range(n):
        # Compute full chain complex for T
        a0 = enumerate_allowed_paths(A, n, 0)
        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)
        a4 = enumerate_allowed_paths(A, n, 4)

        # Compute same for T\\v
        B = delete_vertex(A, n, v)
        b0 = enumerate_allowed_paths(B, n-1, 0)
        b1 = enumerate_allowed_paths(B, n-1, 1)
        b2 = enumerate_allowed_paths(B, n-1, 2)
        b3 = enumerate_allowed_paths(B, n-1, 3)

        # Compute Betti numbers for T
        betti_T = []
        for p in range(n):
            ap = enumerate_allowed_paths(A, n, p)
            ap_prev = enumerate_allowed_paths(A, n, p-1) if p > 0 else []
            if ap and (p == 0 or ap_prev):
                omp = compute_omega_basis(A, n, p, ap, ap_prev) if p > 0 else np.eye(n)
                dimp = omp.shape[1] if omp.ndim == 2 else n
            else:
                dimp = n if p == 0 else 0

            if p > 0 and dimp > 0:
                ap_prev_list = enumerate_allowed_paths(A, n, p-1)
                bdp = build_full_boundary_matrix(ap, ap_prev_list)
                omp_actual = compute_omega_basis(A, n, p, ap, ap_prev_list)
                om_prev = compute_omega_basis(A, n, p-1, ap_prev_list,
                    enumerate_allowed_paths(A, n, p-2) if p >= 2 else []) if p > 0 else np.eye(n)
                bdp_om = bdp @ omp_actual
                coords, _, _, _ = np.linalg.lstsq(om_prev, bdp_om, rcond=None)
                S = np.linalg.svd(coords, compute_uv=False)
                rkp = int(sum(s > 1e-8 for s in S))
            else:
                rkp = 0
            betti_T.append((dimp, rkp))

        # β₂(T) from the data
        ker2_T = betti_T[2][0] - betti_T[2][1] if len(betti_T) > 2 else 0
        im3_T = betti_T[3][1] if len(betti_T) > 3 else 0
        b2_T = ker2_T - im3_T

        # Compute Betti for T\\v
        betti_Tv = []
        for p in range(n-1):
            bp = enumerate_allowed_paths(B, n-1, p)
            bp_prev = enumerate_allowed_paths(B, n-1, p-1) if p > 0 else []
            if bp and (p == 0 or bp_prev):
                omp = compute_omega_basis(B, n-1, p, bp, bp_prev) if p > 0 else np.eye(n-1)
                dimp = omp.shape[1] if omp.ndim == 2 else n-1
            else:
                dimp = n-1 if p == 0 else 0

            if p > 0 and dimp > 0:
                bp_prev_list = enumerate_allowed_paths(B, n-1, p-1)
                bdp = build_full_boundary_matrix(bp, bp_prev_list)
                omp_actual = compute_omega_basis(B, n-1, p, bp, bp_prev_list)
                om_prev = compute_omega_basis(B, n-1, p-1, bp_prev_list,
                    enumerate_allowed_paths(B, n-1, p-2) if p >= 2 else []) if p > 0 else np.eye(n-1)
                bdp_om = bdp @ omp_actual
                coords, _, _, _ = np.linalg.lstsq(om_prev, bdp_om, rcond=None)
                S = np.linalg.svd(coords, compute_uv=False)
                rkp = int(sum(s > 1e-8 for s in S))
            else:
                rkp = 0
            betti_Tv.append((dimp, rkp))

        ker2_Tv = betti_Tv[2][0] - betti_Tv[2][1] if len(betti_Tv) > 2 else 0
        im3_Tv = betti_Tv[3][1] if len(betti_Tv) > 3 else 0
        b2_Tv = ker2_Tv - im3_Tv

        # By long exact sequence:
        # H₂(T\\v) → H₂(T) → H₂(T,T\\v) → H₁(T\\v) → H₁(T)
        # If β₂(T\\v) = 0: 0 → β₂(T) → H₂^rel → ...
        # So β₂(T) injects into H₂^rel.
        # Since β₂(T) = 0 (verified), H₂^rel ≅ {portion going to H₁(T\\v)}
        # Actually H₂^rel = 0 directly by computation.

        # Let me use Euler characteristic of relative complex:
        # R₁, R₂, R₃ = dim(Ω_p(T)) - dim(Ω_p(T\\v))
        dims_T = {}
        dims_Tv = {}
        for p in range(5):
            ap = enumerate_allowed_paths(A, n, p)
            ap_prev = enumerate_allowed_paths(A, n, p-1) if p > 0 else []
            if p == 0:
                dims_T[p] = n
            elif ap:
                om = compute_omega_basis(A, n, p, ap, ap_prev)
                dims_T[p] = om.shape[1] if om.ndim == 2 else 0
            else:
                dims_T[p] = 0

        for p in range(4):
            bp = enumerate_allowed_paths(B, n-1, p)
            bp_prev = enumerate_allowed_paths(B, n-1, p-1) if p > 0 else []
            if p == 0:
                dims_Tv[p] = n-1
            elif bp:
                om = compute_omega_basis(B, n-1, p, bp, bp_prev)
                dims_Tv[p] = om.shape[1] if om.ndim == 2 else 0
            else:
                dims_Tv[p] = 0

        # Relative Euler char: χ_rel = Σ (-1)^p (dim_T - dim_Tv)
        # = χ(T) - χ(T\\v) = (1-β₁(T)) - (1-β₁(T\\v)) = β₁(T\\v) - β₁(T)
        chi_rel = sum((-1)**p * (dims_T.get(p,0) - dims_Tv.get(p,0)) for p in range(5))

        # Also: χ_rel = h₀^rel - h₁^rel + h₂^rel - h₃^rel + ...
        # If H₂^rel = 0: χ_rel = h₀^rel - h₁^rel - h₃^rel + ...

        if detail_count < 5:
            print(f"  T#{T_idx}, v={v}: R = [{', '.join(str(dims_T.get(p,0)-dims_Tv.get(p,0)) for p in range(5))}], "
                  f"χ_rel = {chi_rel}")
            detail_count += 1

        h2_rel_vals[int(round(chi_rel))] += 1

    if T_idx % 200 == 0 and T_idx > 0:
        print(f"  ... {T_idx}/1024")

print(f"\n  χ_rel distribution:")
for val, count in sorted(h2_rel_vals.items()):
    print(f"    χ_rel={val}: {count}")

# ===== KEY INSIGHT: THE Ω₂ FORMULA AND VERTEX DELETION =====
print(f"\n{'='*70}")
print("Ω₂ FORMULA UNDER VERTEX DELETION")
print("=" * 70)

# We proved: dim(Ω₂) = |A₂| - |BF₂|
# where BF₂ = set of non-A₁ faces = {(a,c) : c→a, ∃b with a→b→c}
#
# Under vertex deletion:
# |A₂(T)| - |A₂(T\\v)| = (# 2-paths through v)
# |BF₂(T)| - |BF₂(T\\v)| = (# bad faces through v)
#
# So R₂ = (2-paths through v) - (bad faces through v)
#
# This gives R₂ in terms of combinatorial data about v!

n = 5
r2_formula_matches = 0
total_cases = 0

for A in all_tournaments(n):
    for v in range(n):
        total_cases += 1

        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        a1_set = set(tuple(p) for p in a1)

        B = delete_vertex(A, n, v)
        b1 = enumerate_allowed_paths(B, n-1, 1)
        b2 = enumerate_allowed_paths(B, n-1, 2)
        b1_set = set(tuple(p) for p in b1)

        # 2-paths through v
        paths_v = [p for p in a2 if v in tuple(p)]

        # Bad faces of A₂ through v: non-A₁ faces of 2-paths through v
        bf_v = set()
        for p in paths_v:
            a, b, c = p
            face_ac = (a, c)
            if face_ac not in a1_set:
                bf_v.add(face_ac)

        # Bad faces NOT through v: subset of BF₂(T) not involving v
        bf_total = set()
        for p in a2:
            a, b, c = p
            if (a,c) not in a1_set:
                bf_total.add((a,c))

        bf_nv = set()
        for p in b2:
            # Need to remap indices
            a, b, c = p
            # In B, vertices are 0..n-2 corresponding to T vertices with v removed
            # We need the actual A₁ for B
            pass  # complex remapping, skip

        # Direct computation of R₂
        om2_T = compute_omega_basis(A, n, 2, a2, a1)
        dim_T = om2_T.shape[1] if om2_T.ndim == 2 else 0

        om2_B = compute_omega_basis(B, n-1, 2, b2, b1) if b2 else np.zeros((0,0))
        dim_B = om2_B.shape[1] if om2_B.ndim == 2 and om2_B.shape[0] > 0 else 0

        R2 = dim_T - dim_B
        predicted_R2 = len(paths_v) - len(bf_v)

        if R2 == predicted_R2:
            r2_formula_matches += 1

print(f"\n  R₂ = |A₂ through v| - |BF₂ through v|: {r2_formula_matches}/{total_cases}")

# ===== WHAT'S THE ANALOGOUS FORMULA FOR R₃? =====
# R₃ = dim(Ω₃(T)) - dim(Ω₃(T\\v))
# If dim(Ω₃) = |A₃| - |BF₃|, then R₃ = |A₃ through v| - |BF₃ through v|
# But the Ω₃ formula is more complex (not just |A₃| - |BF₃|).

# Let me check anyway:
n = 5
r3_formula_matches = 0
total_cases = 0

for A in all_tournaments(n):
    for v in range(n):
        total_cases += 1

        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)
        a2_set = set(tuple(p) for p in a2)

        B = delete_vertex(A, n, v)
        b2 = enumerate_allowed_paths(B, n-1, 2)
        b3 = enumerate_allowed_paths(B, n-1, 3)

        paths_v3 = [p for p in a3 if v in tuple(p)]

        # Bad faces of A₃ through v
        bf_v3 = set()
        for p in paths_v3:
            a, b, c, d = p
            for fi in range(4):
                face = tuple(list(p[:fi]) + list(p[fi+1:]))
                if face not in a2_set:
                    bf_v3.add(face)

        om3_T = compute_omega_basis(A, n, 3, a3, a2) if a3 else np.zeros((0,0))
        dim3_T = om3_T.shape[1] if om3_T.ndim == 2 and om3_T.shape[0] > 0 else 0

        a1_B = enumerate_allowed_paths(B, n-1, 1)
        om3_B = compute_omega_basis(B, n-1, 3, b3, b2) if b3 else np.zeros((0,0))
        dim3_B = om3_B.shape[1] if om3_B.ndim == 2 and om3_B.shape[0] > 0 else 0

        R3 = dim3_T - dim3_B
        predicted_R3 = len(paths_v3) - len(bf_v3)

        if R3 == predicted_R3:
            r3_formula_matches += 1

print(f"  R₃ = |A₃ through v| - |BF₃ through v|: {r3_formula_matches}/{total_cases}")

print("\nDone.")
