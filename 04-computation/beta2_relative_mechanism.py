#!/usr/bin/env python3
"""
beta2_relative_mechanism.py — WHY does H₂(T,T\\v)=0 for tournaments?

KEY OBSERVATIONS:
1. For transitive T_n: R_p = C(n-1, p) (simplex boundary)
   The relative complex is the simplicial chain complex of Δ^{n-2}.
   H₂(Δ^{n-2}) = 0 for n≥4 (simplex is contractible).

2. For general tournaments: R differs from binomial but H₂^rel=0 persists.

3. Tournament completeness is KEY: every pair has an arc.
   For vertex v with out-degree d and in-degree n-1-d:
   - 2-paths starting at v: (v,b,c) for b∈out(v), c∈out(b)\\{v}
   - 2-paths through v at pos 1: (a,v,c) for a∈in(v), c∈out(v)
   - 2-paths ending at v: (a,b,v) for b∈in(v), a∈in(b)\\{v}

The position-1 paths (a,v,c) are crucial: they connect in(v) to out(v).
In a tournament, |in(v)|+|out(v)| = n-1, so the number of pos-1 paths
is |in(v)|*|out(v)| ≤ ((n-1)/2)² = (n-1)²/4.

STRATEGY: Show that for EACH relative 2-cycle z through v, there exists
a relative 3-chain w through v with ∂₃^rel(w) = z.

For this, we need enough 3-paths through v, with the RIGHT boundary structure.

INSIGHT: A 3-path through v at position i has:
  - 3 faces through v (relative part)
  - 1 face NOT through v (quotient part, ignored)

The boundary ∂₃^rel only sees the 3 faces through v.
This is a SIMPLER map than ∂₃ itself.

Let me compute ker(∂₂^rel) and im(∂₃^rel) EXACTLY and find the pattern.

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
        if i == v: continue
        row = []
        for j in range(n):
            if j == v: continue
            row.append(A[i][j])
        B.append(row)
    return B

print("=" * 70)
print("RELATIVE MECHANISM: EXACT ker/im COMPUTATION")
print("=" * 70)

n = 5

# For each (T,v), compute the EXACT relative chain complex and verify H₂^rel=0.
# Use the quotient construction: R_p = Ω_p(T) / Ω_p(T\\v).
# Represent R_p by choosing a complement to Ω_p(T\\v) in Ω_p(T).

h2_rel_check = 0
h2_rel_fail = 0
total = 0
ker_im_data = []

for T_idx, A in enumerate(all_tournaments(n)):
    if T_idx % 200 == 0 and T_idx > 0:
        print(f"  ... {T_idx}/1024")

    for v in range(n):
        total += 1

        # --- Compute Ω_p for T and T\\v ---
        B = delete_vertex(A, n, v)

        # Level p: compute Ω₂, Ω₃ for both
        a1_T = enumerate_allowed_paths(A, n, 1)
        a2_T = enumerate_allowed_paths(A, n, 2)
        a3_T = enumerate_allowed_paths(A, n, 3)

        b1 = enumerate_allowed_paths(B, n-1, 1)
        b2 = enumerate_allowed_paths(B, n-1, 2)
        b3 = enumerate_allowed_paths(B, n-1, 3)

        om2_T = compute_omega_basis(A, n, 2, a2_T, a1_T)
        om3_T = compute_omega_basis(A, n, 3, a3_T, a2_T) if a3_T else np.zeros((len(a3_T) if a3_T else 0, 0))

        om2_B = compute_omega_basis(B, n-1, 2, b2, b1) if b2 else np.zeros((0, 0))
        om3_B = compute_omega_basis(B, n-1, 3, b3, b2) if b3 else np.zeros((0, 0))

        dim_om2_T = om2_T.shape[1] if om2_T.ndim == 2 else 0
        dim_om3_T = om3_T.shape[1] if om3_T.ndim == 2 and om3_T.shape[0] > 0 else 0
        dim_om2_B = om2_B.shape[1] if om2_B.ndim == 2 and om2_B.shape[0] > 0 else 0
        dim_om3_B = om3_B.shape[1] if om3_B.ndim == 2 and om3_B.shape[0] > 0 else 0

        R2 = dim_om2_T - dim_om2_B
        R3 = dim_om3_T - dim_om3_B

        if R2 == 0:
            h2_rel_check += 1
            continue

        # --- Build relative boundary matrices ---
        # We need the boundary map ∂₂: Ω₂(T) → Ω₁(T) in Ω-coordinates,
        # and then restrict/project to the relative part.
        #
        # Simpler approach: use the full boundary matrix in A_p coordinates
        # and identify the v-components.

        # Create index maps
        a2_T_list = [tuple(p) for p in a2_T]
        a3_T_list = [tuple(p) for p in a3_T]

        # v-paths: those containing v
        v_idx_2 = [i for i, p in enumerate(a2_T_list) if v in p]
        nv_idx_2 = [i for i, p in enumerate(a2_T_list) if v not in p]
        v_idx_3 = [i for i, p in enumerate(a3_T_list) if v in p]

        # Boundary ∂₃: A₃ → A₂
        if a3_T:
            bd3 = build_full_boundary_matrix(a3_T, a2_T)
            # ∂₃ restricted to v-paths in source, projected to v-paths in target
            bd3_rel = bd3[np.ix_(v_idx_2, v_idx_3)]
        else:
            bd3_rel = np.zeros((len(v_idx_2), 0))

        # Boundary ∂₂: A₂ → A₁
        bd2 = build_full_boundary_matrix(a2_T, a1_T)
        a1_T_list = [tuple(p) for p in a1_T]
        v_idx_1 = [i for i, p in enumerate(a1_T_list) if v in p]
        # ∂₂ restricted to v-paths source, projected to v-edges target
        bd2_rel = bd2[np.ix_(v_idx_1, v_idx_2)]

        # Now: the relative complex uses Ω-coordinates, not A-coordinates.
        # But for the relative part, we need:
        # ker(∂₂^rel) in the v-path subspace of Ω₂
        # im(∂₃^rel) in the same subspace

        # The issue: Ω₂ ≠ A₂ for tournaments with NT triples.
        # The v-path restriction of Ω₂ is not just the v-rows of om2_T.

        # However, for the relative complex, we're working modulo Ω₂(T\\v).
        # The quotient R₂ = Ω₂(T)/Ω₂(T\\v) is a quotient space.

        # Let me use a direct approach: compute H₂ of the relative complex
        # by computing H₂(T) and H₂(T\\v) and using the LES.

        # Actually, the simplest approach: compute β₂ of both T and T\\v,
        # then use the LES.

        # For n=5, β₂(T) = 0 and β₂(T\\v) = 0 (4-vertex tournaments).
        # LES: 0 → 0 → H₂^rel → H₁(T\\v) → H₁(T) → ...
        # So H₂^rel injects into H₁(T\\v).
        # If H₂^rel ≠ 0, then there's a nontrivial map H₂^rel → H₁(T\\v).

        # For H₂^rel = 0: the connecting homomorphism δ: H₂^rel → H₁(T\\v) is zero.
        # This means every relative 2-cycle lifts to an absolute 2-chain in T\\v
        # whose boundary is an absolute 1-cycle... this is getting complex.

        # Let me just compute numerically using the raw relative boundary matrices.
        # ker(∂₂^rel) and im(∂₃^rel) in A-coordinates for v-paths.

        # ker(∂₂^rel): v-path combinations with zero v-edge boundary
        if bd2_rel.shape[0] > 0 and bd2_rel.shape[1] > 0:
            _, S2, Vt2 = np.linalg.svd(bd2_rel, full_matrices=True)
            rk2 = int(np.sum(S2 > 1e-8))
            # Kernel of bd2_rel: last (cols - rk) rows of Vt2
            ker2 = Vt2[rk2:].T  # columns are kernel basis
            dim_ker2 = ker2.shape[1] if ker2.ndim == 2 else 0
        else:
            dim_ker2 = len(v_idx_2)
            ker2 = np.eye(len(v_idx_2))

        # im(∂₃^rel): range of bd3_rel
        if bd3_rel.shape[1] > 0:
            rk3 = np.linalg.matrix_rank(bd3_rel, tol=1e-8)
        else:
            rk3 = 0

        # H₂^rel = dim(ker2) - rk3
        # But this is in A-coordinates, not Ω-coordinates!
        # For the correct computation, we'd need to intersect with Ω₂.
        # However, for tournaments: Ω₂ ≠ A₂ in general.

        # Correct approach: work entirely in Ω-coordinates.
        # This requires embedding T\\v paths into T's path space.
        # Complex but doable...

        # SIMPLER: Use the formula H₂^rel = β₂(T) - β₂(T\\v) + (connecting maps)
        # From the LES: ... H₂(T\\v) → H₂(T) → H₂^rel → H₁(T\\v) → H₁(T) → ...
        # With β₂(T)=β₂(T\\v)=0:
        # 0 → 0 → H₂^rel → H₁(T\\v) → H₁(T) → H₁^rel → H₀(T\\v) → ...

        # So H₂^rel is the kernel of δ: H₂^rel → H₁(T\\v).
        # Hmm, that's the definition, not helpful.

        # Actually: H₂^rel ≅ ker(H₁(T\\v) → H₁(T)) since the sequence is
        # 0 → H₂^rel → H₁(T\\v) → H₁(T)
        # So H₂^rel ≅ ker(i_*: H₁(T\\v) → H₁(T)) where i is inclusion.

        # This is GREAT! H₂^rel = 0 iff the inclusion-induced map
        # H₁(T\\v) → H₁(T) is INJECTIVE.

        # Let's verify: does inclusion T\\v ↪ T induce injection on H₁?

        # Compute β₁(T) and β₁(T\\v)
        om1_T = compute_omega_basis(A, n, 1, a1_T, enumerate_allowed_paths(A, n, 0))
        bd1_T = build_full_boundary_matrix(a1_T, enumerate_allowed_paths(A, n, 0))
        bd1_om_T = bd1_T @ om1_T
        S_bd1 = np.linalg.svd(bd1_om_T, compute_uv=False)
        rk_bd1_T = int(sum(s > 1e-8 for s in S_bd1))

        bd2_om_T = bd2 @ om2_T
        coords_bd2, _, _, _ = np.linalg.lstsq(om1_T, bd2_om_T, rcond=None)
        S_bd2 = np.linalg.svd(coords_bd2, compute_uv=False)
        im2_T = int(sum(s > 1e-8 for s in S_bd2))

        beta1_T = (om1_T.shape[1] - rk_bd1_T) - im2_T

        om1_B = compute_omega_basis(B, n-1, 1, b1, enumerate_allowed_paths(B, n-1, 0))
        bd1_B = build_full_boundary_matrix(b1, enumerate_allowed_paths(B, n-1, 0))
        bd1_om_B = bd1_B @ om1_B
        S_bd1_B = np.linalg.svd(bd1_om_B, compute_uv=False)
        rk_bd1_B = int(sum(s > 1e-8 for s in S_bd1_B))

        if b2 and om2_B.ndim == 2 and om2_B.shape[1] > 0:
            bd2_B = build_full_boundary_matrix(b2, b1)
            bd2_om_B = bd2_B @ om2_B
            coords_bd2_B, _, _, _ = np.linalg.lstsq(om1_B, bd2_om_B, rcond=None)
            S_bd2_B = np.linalg.svd(coords_bd2_B, compute_uv=False)
            im2_B = int(sum(s > 1e-8 for s in S_bd2_B))
        else:
            im2_B = 0

        beta1_B = (om1_B.shape[1] - rk_bd1_B) - im2_B

        # H₂^rel = ker(i_*: H₁(T\\v) → H₁(T))
        # dim H₂^rel = β₁(T\\v) - rk(i_*)
        # where rk(i_*) = β₁(T\\v) - dim(ker i_*)

        # For H₂^rel = 0: need i_* injective, i.e., ker(i_*) = 0.
        # This means: every 1-cycle in T\\v that becomes a boundary in T
        # was already a boundary in T\\v.

        # Equivalently: β₁(T\\v) ≤ β₁(T) + something from H₁^rel...
        # Actually, the LES gives:
        # 0 → H₂^rel → H₁(T\\v) →^{i_*} H₁(T) → H₁^rel → H₀(T\\v) → H₀(T)
        # Since H₀(T\\v) = Z and H₀(T) = Z (both connected):
        # H₁^rel → Z → Z is eventually exact.
        # rk(H₁(T) → H₁^rel) + rk(H₁^rel → H₀(T\\v)) = dim H₁^rel

        # From 0 → H₂^rel → H₁(T\\v) → H₁(T):
        # dim H₂^rel = dim ker(i_*)
        # dim im(i_*) = β₁(T\\v) - dim H₂^rel

        # We know H₂^rel = 0 computationally. So i_* is injective:
        # β₁(T\\v) ≤ β₁(T). Let's verify!

        ker_im_data.append({
            'beta1_T': beta1_T,
            'beta1_Tv': beta1_B,
            'R2': R2, 'R3': R3,
            'v': v, 'dv': sum(A[v]),
            'injective': beta1_B <= beta1_T
        })

        if beta1_B <= beta1_T:
            h2_rel_check += 1
        else:
            h2_rel_fail += 1

# KEY CHECK: Is i_*: H₁(T\\v) → H₁(T) always injective?
print(f"\n  H₂^rel=0 iff i_* injective iff β₁(T\\v) ≤ β₁(T)")
print(f"  Injective: {h2_rel_check}/{total} ({100*h2_rel_check/total:.1f}%)")
print(f"  Non-injective: {h2_rel_fail}/{total} ({100*h2_rel_fail/total:.1f}%)")

# Wait — earlier we found β₁(T\\v) > β₁(T) in 840/5120 cases!
# But H₂^rel = 0 in ALL cases. This means my LES argument is WRONG somewhere.

# Let me recheck: the map H₂(T\\v) → H₂(T) might not be zero even when both are 0.
# Actually: H₂(T\\v) = 0 and H₂(T) = 0, so the map IS zero.
# The LES gives: 0 → H₂(T) → H₂^rel → H₁(T\\v) → H₁(T)
# i.e., 0 → 0 → H₂^rel → H₁(T\\v) → H₁(T)
# So H₂^rel ↪ H₁(T\\v) is exact at H₂^rel, meaning H₂^rel ≅ ker(i_*).

# If β₁(T\\v) > β₁(T), then i_* CANNOT be injective (since it maps a
# bigger space into a smaller one, UNLESS the map is not between vector spaces
# of those dimensions).

# WAIT: i_* is a LINEAR MAP H₁(T\\v) → H₁(T). Its kernel has dimension
# ≥ β₁(T\\v) - β₁(T) when β₁(T\\v) > β₁(T).
# But H₂^rel = ker(i_*) would then be ≥ β₁(T\\v) - β₁(T) > 0.
# This contradicts H₂^rel = 0!

# Unless my LES is WRONG. Let me check the GLMY long exact sequence more carefully.

print(f"\n  Cases with β₁(T\\v) > β₁(T):")
violations = [(d['beta1_T'], d['beta1_Tv'], d['v'], d['dv'])
              for d in ker_im_data if d['beta1_Tv'] > d['beta1_T']]
print(f"  Count: {len(violations)}")

beta_diff = Counter()
for d in ker_im_data:
    beta_diff[d['beta1_Tv'] - d['beta1_T']] += 1
print(f"\n  β₁(T\\v) - β₁(T) distribution:")
for diff, count in sorted(beta_diff.items()):
    print(f"    diff={diff}: {count}")

# This means either:
# 1. The LES for GLMY path homology doesn't work the same as for simplicial
# 2. My computation of β₁ is wrong
# 3. The LES gives 0 → H₂^rel → H₁(T\\v) →^{i_*} H₁(T) but i_* may NOT be
#    the inclusion-induced map on homology
#
# Actually (3) is the issue. In the GLMY theory, the LES for a pair (X,Y)
# is defined via the SHORT EXACT SEQUENCE 0 → Ω_*(Y) → Ω_*(X) → R_* → 0.
# The connecting homomorphism δ: H_p^rel → H_{p-1}(Y) is NOT the same as
# the inclusion-induced map. It involves a zig-zag construction.
#
# So H₂^rel = 0 does NOT mean i_* is injective.
# The actual LES is:
#   H_2(Y) → H_2(X) → H_2(X,Y) →^δ H_1(Y) → H_1(X) → H_1(X,Y) → ...
#
# With H₂(X)=H₂(Y)=0:
#   0 → 0 → H₂(X,Y) →^δ H₁(Y) → H₁(X) → ...
# So δ is injective from H₂(X,Y) into H₁(Y).
# And ker(H₁(Y) → H₁(X)) = im(δ) = H₂(X,Y).
# So H₂(X,Y) ≅ ker(H₁(Y) → H₁(X)).
#
# The map H₁(Y) → H₁(X) is indeed the inclusion-induced map.
# So H₂^rel = ker(i_*) EXACTLY.
#
# If β₁(T\\v) > β₁(T), then dim ker(i_*) ≥ β₁(T\\v) - β₁(T) > 0.
# This would mean H₂^rel > 0, contradicting our computation!
#
# RESOLUTION: My β₁ computation might have errors. Let me recheck.

print(f"\n{'='*70}")
print("DOUBLE-CHECKING β₁ COMPUTATION")
print("=" * 70)

# Use a simple direct computation
n = 5
discrepancies = 0
for T_idx, A in enumerate(all_tournaments(n)):
    if T_idx >= 10:
        break

    for v in range(n):
        B = delete_vertex(A, n, v)

        # β₁(T) via direct rank computation
        a1_T = enumerate_allowed_paths(A, n, 1)
        a2_T = enumerate_allowed_paths(A, n, 2)
        a0_T = enumerate_allowed_paths(A, n, 0)

        om1 = compute_omega_basis(A, n, 1, a1_T, a0_T)
        om2 = compute_omega_basis(A, n, 2, a2_T, a1_T)

        # ∂₁ on Ω₁ → Ω₀
        bd1 = build_full_boundary_matrix(a1_T, a0_T)
        bd1_om = bd1 @ om1
        rk1 = np.linalg.matrix_rank(bd1_om, tol=1e-8)

        # ∂₂ on Ω₂ → Ω₁
        bd2 = build_full_boundary_matrix(a2_T, a1_T)
        bd2_om = bd2 @ om2
        # Image of ∂₂ in Ω₁ coordinates
        coords, _, _, _ = np.linalg.lstsq(om1, bd2_om, rcond=None)
        rk2_in_om1 = np.linalg.matrix_rank(coords, tol=1e-8)

        ker1 = om1.shape[1] - rk1
        beta1_T = ker1 - rk2_in_om1

        # β₁(T\\v)
        b1 = enumerate_allowed_paths(B, n-1, 1)
        b2 = enumerate_allowed_paths(B, n-1, 2)
        b0 = enumerate_allowed_paths(B, n-1, 0)

        om1_B = compute_omega_basis(B, n-1, 1, b1, b0)
        bd1_B = build_full_boundary_matrix(b1, b0)
        bd1_om_B = bd1_B @ om1_B
        rk1_B = np.linalg.matrix_rank(bd1_om_B, tol=1e-8)

        if b2:
            om2_B = compute_omega_basis(B, n-1, 2, b2, b1)
            bd2_B = build_full_boundary_matrix(b2, b1)
            bd2_om_B = bd2_B @ om2_B
            coords_B, _, _, _ = np.linalg.lstsq(om1_B, bd2_om_B, rcond=None)
            rk2_B = np.linalg.matrix_rank(coords_B, tol=1e-8)
        else:
            rk2_B = 0

        ker1_B = om1_B.shape[1] - rk1_B
        beta1_Tv = ker1_B - rk2_B

        if beta1_Tv > beta1_T:
            scores = sorted(sum(A[i]) for i in range(n))
            scores_B = sorted(sum(B[i]) for i in range(n-1))
            print(f"  T#{T_idx}, v={v}: β₁(T)={beta1_T}, β₁(T\\v)={beta1_Tv}, "
                  f"scores_T={scores}, scores_Tv={scores_B}")
            discrepancies += 1

if discrepancies == 0:
    print(f"  No β₁(T\\v) > β₁(T) in first 10 tournaments")
else:
    print(f"\n  Found {discrepancies} cases with β₁(T\\v) > β₁(T)")
    print(f"  This means either:")
    print(f"  1. H₂^rel ≠ 0 (contradicts computation)")
    print(f"  2. The LES argument has an error")
    print(f"  3. GLMY LES works differently than simplicial LES")

print("\nDone.")
