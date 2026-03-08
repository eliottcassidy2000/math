#!/usr/bin/env python3
"""
beta2_ses_check.py — Check if Ω_*(T\\v) → Ω_*(T) is injective

For the LES of relative homology to work, we need:
  0 → Ω_p(T\\v) → Ω_p(T) → R_p → 0

This requires Ω_p(T\\v) ↪ Ω_p(T) to be injective.

In simplicial homology, subcomplex inclusion is trivially injective.
But in GLMY path homology, Ω_p is a QUOTIENT of A_p by higher-boundary
conditions. The inclusion A_p(T\\v) ⊆ A_p(T) is obvious, but the
Ω-construction might not preserve this.

Specifically: c ∈ Ω_p(T\\v) means c is an A_p(T\\v) chain with
∂_p(c) ∈ Ω_{p-1}(T\\v). Viewed as an element of A_p(T), we need
∂_p(c) ∈ Ω_{p-1}(T). Is Ω_{p-1}(T\\v) ⊆ Ω_{p-1}(T)?

For p=1: Ω₁(T\\v) = A₁(T\\v) ⊆ A₁(T) = Ω₁(T). ✓
For p=2: Ω₂(T\\v) ⊆ Ω₂(T)?
  c ∈ Ω₂(T\\v) means ∂₂(c) ∈ Ω₁(T\\v) = A₁(T\\v).
  Viewing c in A₂(T): ∂₂(c) is the same chain (faces of T\\v paths are
  T\\v edges, which are also T edges).
  So ∂₂(c) ∈ A₁(T\\v) ⊆ A₁(T) = Ω₁(T). ✓

By induction: Ω_p(T\\v) ⊆ Ω_p(T). So the map IS injective.

Then the LES should hold. But we found β₁(T\\v) > β₁(T) with H₂^rel=0.

Wait — maybe my H₂^rel computation is wrong. Let me compute it directly
from the quotient complex.

Author: opus-2026-03-08-S44
"""
import sys, time
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

def compute_betti(A, n, max_p=None):
    """Compute all Betti numbers for digraph A."""
    if max_p is None:
        max_p = n - 1
    betti = []
    dims = []
    ranks = []

    for p in range(max_p + 1):
        ap = enumerate_allowed_paths(A, n, p)
        ap_prev = enumerate_allowed_paths(A, n, p-1) if p > 0 else []

        if p == 0:
            dim_p = n
        elif ap:
            om = compute_omega_basis(A, n, p, ap, ap_prev)
            dim_p = om.shape[1] if om.ndim == 2 else 0
        else:
            dim_p = 0
        dims.append(dim_p)

        if p == 0 or dim_p == 0:
            ranks.append(0)
            continue

        bd = build_full_boundary_matrix(ap, ap_prev)
        om = compute_omega_basis(A, n, p, ap, ap_prev)
        om_prev = compute_omega_basis(A, n, p-1, ap_prev,
                                       enumerate_allowed_paths(A, n, p-2) if p >= 2 else [])
        bd_om = bd @ om
        coords, _, _, _ = np.linalg.lstsq(om_prev, bd_om, rcond=None)
        S = np.linalg.svd(coords, compute_uv=False)
        rk = int(sum(s > 1e-8 for s in S))
        ranks.append(rk)

    for p in range(max_p + 1):
        ker = dims[p] - ranks[p]
        im_next = ranks[p+1] if p+1 <= max_p else 0
        betti.append(ker - im_next)

    return betti, dims, ranks

print("=" * 70)
print("DIRECT RELATIVE HOMOLOGY CHECK")
print("=" * 70)

# Compute H₂(T,T\\v) directly from the quotient complex.
# R_p = Ω_p(T) / Ω_p(T\\v)
# The induced boundary ∂_p^R: R_p → R_{p-1} is well-defined.
# H₂^R = ker(∂₂^R) / im(∂₃^R)

n = 5
violations = 0
les_violations = 0
total = 0

pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(pairs)

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    betti_T, dims_T, ranks_T = compute_betti(A, n, 4)

    for v in range(n):
        total += 1
        B = delete_vertex(A, n, v)
        betti_Tv, dims_Tv, ranks_Tv = compute_betti(B, n-1, 3)

        # Check LES consistency:
        # ... → H₂(T\\v) → H₂(T) → H₂^rel → H₁(T\\v) → H₁(T) → ...
        # 0 → 0 → H₂^rel → β₁(T\\v) → β₁(T)
        # dim(H₂^rel) = dim ker(i_*: H₁(T\\v) → H₁(T))
        # ≥ max(0, β₁(T\\v) - β₁(T))

        if betti_Tv[1] > betti_T[1]:
            # This forces H₂^rel ≥ β₁(T\\v) - β₁(T) > 0
            # But kind-pasteur verified H₂^rel = 0 at n≤9.
            # CONTRADICTION?
            les_violations += 1
            if les_violations <= 5:
                scores = sorted(sum(A[i]) for i in range(n))
                scores_B = sorted(sum(B[i]) for i in range(n-1))
                print(f"  LES violation: T(scores={scores}), v={v}, "
                      f"β₁(T)={betti_T[1]}, β₁(T\\v)={betti_Tv[1]}, "
                      f"scores_Tv={scores_B}")
                print(f"    dims_T={dims_T}, ranks_T={ranks_T}, betti_T={betti_T}")
                print(f"    dims_Tv={dims_Tv}, ranks_Tv={ranks_Tv}, betti_Tv={betti_Tv}")

    if bits % 200 == 0 and bits > 0:
        print(f"  ... {bits}/{1<<m}, LES violations so far: {les_violations}")

print(f"\n  Total: {total}, LES apparent violations: {les_violations}")

if les_violations > 0:
    print(f"""
  RESOLUTION: The LES of a pair in GLMY theory requires:
  0 → Ω_*(Y) →^i Ω_*(X) →^π R_* → 0 to be EXACT.

  We showed Ω_p(T\\v) ⊆ Ω_p(T) (injection), and R_p is the quotient.
  So the SES should be exact.

  But maybe the BOUNDARY MAP on the quotient is NOT the naive quotient map.
  In GLMY theory, the boundary maps ∂_p depend on the GRAPH structure.
  The quotient boundary ∂^R might differ from the standard quotient.

  Actually: the boundary ∂_p on Ω_p is independent of the graph — it's the
  standard simplicial face map. The graph only affects WHICH paths are allowed
  and WHICH elements are in Ω. So the SES of chain complexes IS exact, and
  the standard LES of homology DOES apply.

  This means our computation of H₂^rel = 0 must be WRONG, or β₁ computation
  has errors. Let me check a specific case...
""")

    # Check a specific case
    print("  --- Checking specific LES violation case ---")
    # T#4, v=4 from earlier: scores=[1,1,2,2,4], β₁(T)=0, β₁(T\\v)=1
    A = [[0]*5 for _ in range(5)]
    # bits=4 → bit pattern 00100 → pair index 2 is set
    # Pair index: (0,1)=0, (0,2)=1, (0,3)=2, (0,4)=3, (1,2)=4, (1,3)=5, (1,4)=6, (2,3)=7, (2,4)=8, (3,4)=9
    for idx, (i,j) in enumerate(pairs):
        if (4 >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    v = 4
    B = delete_vertex(A, n, v)

    print(f"  T adjacency:")
    for i in range(n):
        print(f"    {A[i]}")
    print(f"  T\\4 adjacency:")
    for i in range(n-1):
        print(f"    {B[i]}")

    betti_T, dims_T, ranks_T = compute_betti(A, n, 4)
    betti_Tv, dims_Tv, ranks_Tv = compute_betti(B, n-1, 3)

    print(f"  β(T) = {betti_T}")
    print(f"  β(T\\v) = {betti_Tv}")
    print(f"  Ω dims T = {dims_T}")
    print(f"  Ω dims T\\v = {dims_Tv}")

    # Compute H₂^rel DIRECTLY from the quotient complex
    # R_p = Ω_p(T) / Ω_p(T\\v)
    # Need explicit Ω bases

    for p in range(1, 4):
        ap_T = enumerate_allowed_paths(A, n, p)
        ap_prev_T = enumerate_allowed_paths(A, n, p-1)
        bp = enumerate_allowed_paths(B, n-1, p)
        bp_prev = enumerate_allowed_paths(B, n-1, p-1)

        om_T = compute_omega_basis(A, n, p, ap_T, ap_prev_T)
        om_B = compute_omega_basis(B, n-1, p, bp, bp_prev) if bp else np.zeros((0,0))

        dim_T = om_T.shape[1] if om_T.ndim == 2 else 0
        dim_B = om_B.shape[1] if om_B.ndim == 2 and om_B.shape[0] > 0 else 0

        print(f"\n  Ω_{p}(T): dim={dim_T}, Ω_{p}(T\\v): dim={dim_B}, R_{p}={dim_T-dim_B}")

        # Embed T\\v paths into T paths
        # T\\v vertices are 0,1,2,3 (after removing v=4)
        # These correspond to T vertices 0,1,2,3 directly.
        # So bp paths map directly to the same paths in ap_T.

        if dim_B > 0 and dim_T > 0:
            # om_B is |B_p| × dim_B, om_T is |A_p| × dim_T
            # Create embedding matrix: B_p → A_p
            ap_T_list = [tuple(x) for x in ap_T]
            bp_list = [tuple(x) for x in bp]
            ap_T_idx = {p: i for i, p in enumerate(ap_T_list)}

            embed = np.zeros((len(ap_T), len(bp)))
            for j, path in enumerate(bp_list):
                if path in ap_T_idx:
                    embed[ap_T_idx[path], j] = 1
                else:
                    print(f"    PATH {path} from T\\v NOT in T — shouldn't happen for vertex deletion!")

            # Embed Ω_p(T\\v) into A_p(T) coordinates
            embed_om = embed @ om_B  # |A_p(T)| × dim_B

            # Check if this lands in Ω_p(T): project onto om_T
            coords, res, _, _ = np.linalg.lstsq(om_T, embed_om, rcond=None)
            recon = om_T @ coords
            error = np.max(np.abs(embed_om - recon)) if embed_om.size > 0 else 0
            rk_embed = np.linalg.matrix_rank(coords, tol=1e-8)
            print(f"    Embedding error: {error:.2e}, rank of embedding in Ω(T): {rk_embed}")
            if rk_embed < dim_B:
                print(f"    *** INJECTION FAILS: rk={rk_embed} < dim(Ω_B)={dim_B} ***")

print("\nDone.")
