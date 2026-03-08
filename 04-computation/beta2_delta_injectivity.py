#!/usr/bin/env python3
"""
beta2_delta_injectivity.py — Test injectivity of δ: H₂(T,T\v) → H₁(T\v)

PROOF STRATEGY:
  Base: β₂(T) = 0 for n ≤ 3
  Step: For n-vertex tournament T, pick any vertex v.
    - By induction: β₂(T\v) = 0
    - LES: 0 = H₂(T\v) → H₂(T) → H₂(T,T\v) → H₁(T\v)
    - H₂(T\v)=0 makes H₂(T) → H₂(T,T\v) injective
    - If δ: H₂(T,T\v) → H₁(T\v) is injective, then im(H₂(T)→H₂(T,T\v)) ⊂ ker(δ) = 0
    - So H₂(T) = 0.

So β₂=0 reduces to: δ is injective for all T, v.

Equivalently: H₂(T,T\v) injects into H₁(T\v) via the connecting map.

This script verifies this at n=5,6,7.

Author: opus-2026-03-08-S49
"""
import sys, time, random
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


def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


def check_delta_injective(A, n, v):
    """
    Check if δ: H₂(T, T\\v) → H₁(T\\v) is injective.

    The relative chain complex:
    Ω_p(T, T\\v) = Ω_p(T) / Ω_p(T\\v)
    = { [α] ∈ Ω_p(T) : α uses vertex v } / (boundaries from Ω_{p+1} using v)

    The connecting map δ:
    For [z] ∈ H₂(T, T\\v), lift z to ẑ ∈ Ω₂(T), compute ∂₂(ẑ) ∈ Ω₁(T).
    Since [z] is a relative cycle, ∂₂(ẑ) ∈ Ω₁(T\\v).
    δ([z]) = [∂₂(ẑ)] ∈ H₁(T\\v).

    δ injective ⟺ if ∂₂(ẑ) = ∂₂(w) for some w ∈ Ω₂(T\\v), then z = ∂₃(σ) for σ using v.

    Direct computation: build the relative chain complex and compute.
    """
    # Build T\\v
    others = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]

    # Chain complexes for T and T\\v
    ap0_T = enumerate_allowed_paths(A, n, 0)
    ap1_T = enumerate_allowed_paths(A, n, 1)
    ap2_T = enumerate_allowed_paths(A, n, 2)
    ap3_T = enumerate_allowed_paths(A, n, 3)

    ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)
    ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
    ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)
    ap3_sub = enumerate_allowed_paths(A_sub, n1, 3)

    # Remap sub paths to T's vertex labels
    remap = {i: others[i] for i in range(n1)}

    # Ω bases
    om1_T = compute_omega_basis(A, n, 1, ap1_T, ap0_T)
    om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T) if ap2_T else np.zeros((0,0))
    om3_T = compute_omega_basis(A, n, 3, ap3_T, ap2_T) if ap3_T else np.zeros((0,0))

    om1_sub = compute_omega_basis(A_sub, n1, 1, ap1_sub, ap0_sub)
    om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))

    d2_T = dim_om(om2_T)
    d3_T = dim_om(om3_T)

    if d2_T == 0:
        return True, 0  # H₂(T,T\v)=0 trivially

    # Identify which Ω₂(T) elements use vertex v
    # An Ω₂(T) element is in Ω₂(T\v) iff it's a linear combination of paths NOT using v.
    ap2_T_list = [tuple(p) for p in ap2_T]
    uses_v = np.array([1 if v in path else 0 for path in ap2_T_list])

    # The relative Ω₂(T,T\v): we need to work in Ω₂(T) coords.
    # An element α ∈ Ω₂(T) is "relative" if its A₂-expansion uses v.
    # In Ω₂(T) coords: c ∈ R^{d2_T}, the A₂ form is om2_T @ c.
    # The "uses v" part is (uses_v^T @ om2_T) @ c ... hmm, this is complicated.

    # Better approach: directly compute relative H₂ via quotient.
    # Ω₂(T) has basis {e_1,...,e_{d2}} in Ω₂ coords.
    # Ω₂(T\v) is the subspace spanned by Ω₂ elements that don't use v.

    # Find the subspace of Ω₂(T) corresponding to T\v.
    # A path in T\v corresponds to a path in T not using v.
    # We need to express Ω₂(T\v) inside Ω₂(T).

    # Map: for each Ω₂(T\v) basis element (in T\v coords), find its expression in Ω₂(T) coords.
    ap2_sub_list = [tuple(remap[x] for x in p) for p in ap2_sub] if ap2_sub else []

    if ap2_sub and dim_om(om2_sub) > 0:
        d2_sub = dim_om(om2_sub)
        # Embed Ω₂(T\v) into A₂(T) space
        embed = np.zeros((len(ap2_T_list), d2_sub))
        for j in range(d2_sub):
            # j-th Ω₂(T\v) basis in A₂(T\v) coords
            a2_sub_vec = om2_sub[:, j]
            # Map to A₂(T) coords
            for k, path_sub in enumerate(ap2_sub):
                path_T = tuple(remap[x] for x in path_sub)
                if path_T in ap2_T_list:
                    idx_T = ap2_T_list.index(path_T)
                    embed[idx_T, :] += np.zeros(d2_sub)  # initialize
                    embed[idx_T, j] = a2_sub_vec[k]

        # Now express embed in Ω₂(T) coords
        # embed = om2_T @ phi, so phi = lstsq(om2_T, embed)
        phi = np.linalg.lstsq(om2_T, embed, rcond=None)[0]  # d2_T × d2_sub
    else:
        phi = np.zeros((d2_T, 0))
        d2_sub = 0

    # Relative Ω₂ = Ω₂(T) / im(phi)
    # dim(relative Ω₂) = d2_T - rk(phi)
    rk_phi = np.linalg.matrix_rank(phi, tol=1e-8) if phi.shape[1] > 0 else 0

    # Boundary maps in Ω coords
    bd2_T = build_full_boundary_matrix(ap2_T, ap1_T)
    coords2_T = np.linalg.lstsq(om1_T, bd2_T @ om2_T, rcond=None)[0]  # d1_T × d2_T

    if d3_T > 0:
        bd3_T = build_full_boundary_matrix(ap3_T, ap2_T)
        coords3_T = np.linalg.lstsq(om2_T, bd3_T @ om3_T, rcond=None)[0]  # d2_T × d3_T
    else:
        coords3_T = np.zeros((d2_T, 0))

    # The relative boundary ∂₂^rel: Ω₂(T)/Ω₂(T\v) → Ω₁(T)/Ω₁(T\v)
    # For the connecting map δ, we need:
    # 1. Relative Z₂ = ker(∂₂^rel)
    # 2. For each element of relative Z₂, compute δ (= ∂₂ of a lift, projected to T\v)
    # 3. Check if δ maps into B₁(T\v) only when element is in relative B₂

    # Simpler approach: compute H₂(T, T\v) and the connecting map directly.

    # H₂(T, T\v) = ker(d₂: relative → relative) / im(d₃: relative → relative)

    # Use the quotient description:
    # Working in Ω₂(T), modulo im(phi).
    # Find complement of im(phi): SVD of phi to get basis for im(phi)^perp

    if rk_phi > 0:
        U_phi, S_phi, Vt_phi = np.linalg.svd(phi, full_matrices=True)
        # Columns of U_phi[:, :rk_phi] span im(phi) in Ω₂(T)
        # The quotient Ω₂(T)/Ω₂(T\v) is represented by U_phi[:, rk_phi:]
        Q = U_phi[:, rk_phi:]  # d2_T × (d2_T - rk_phi), columns represent relative Ω₂
    else:
        Q = np.eye(d2_T)

    d_rel = Q.shape[1]
    if d_rel == 0:
        return True, 0

    # ∂₂ in relative coords: ∂₂(Qx) projected to relative Ω₁
    # But we don't need the full relative boundary. We just need ker and im.

    # ∂₂^rel: R^{d_rel} → Ω₁(T)/Ω₁(T\v)
    # coords2_rel = coords2_T @ Q  (in full Ω₁(T) coords)
    coords2_rel = coords2_T @ Q

    # The image of Ω₁(T\v) in Ω₁(T): similar embedding
    ap1_sub_list = [tuple(remap[x] for x in p) for p in ap1_sub] if ap1_sub else []
    ap1_T_list = [tuple(p) for p in ap1_T]
    d1_T = dim_om(om1_T)
    d1_sub = dim_om(om1_sub)

    if d1_sub > 0:
        embed1 = np.zeros((len(ap1_T_list), d1_sub))
        for j in range(d1_sub):
            a1_sub_vec = om1_sub[:, j]
            for k, path_sub in enumerate(ap1_sub):
                path_T = tuple(remap[x] for x in path_sub)
                if path_T in ap1_T_list:
                    idx_T = ap1_T_list.index(path_T)
                    embed1[idx_T, j] = a1_sub_vec[k]
        psi = np.linalg.lstsq(om1_T, embed1, rcond=None)[0]  # d1_T × d1_sub
    else:
        psi = np.zeros((d1_T, 0))

    rk_psi = np.linalg.matrix_rank(psi, tol=1e-8) if psi.shape[1] > 0 else 0

    # For the relative ∂₂, we quotient the target by im(psi)
    if rk_psi > 0:
        U_psi, S_psi, _ = np.linalg.svd(psi, full_matrices=True)
        R = U_psi[:, rk_psi:]  # complement of Ω₁(T\v) in Ω₁(T)
    else:
        R = np.eye(d1_T)

    # ∂₂^rel in quotient coords
    coords2_rel_q = R.T @ coords2_rel  # (d1_T - rk_psi) × d_rel
    rk_d2_rel = np.linalg.matrix_rank(coords2_rel_q, tol=1e-8)

    # ker(∂₂^rel) = relative Z₂
    z2_rel_dim = d_rel - rk_d2_rel

    if z2_rel_dim == 0:
        return True, 0  # H₂(T,T\v) = 0

    # Now find im(∂₃^rel)
    # ∂₃: Ω₃(T) → Ω₂(T), then project to quotient
    if d3_T > 0:
        # Which Ω₃ elements use v?
        ap3_T_list = [tuple(p) for p in ap3_T]
        ap3_sub_list = [tuple(remap[x] for x in p) for p in ap3_sub] if ap3_sub else []

        om3_sub = compute_omega_basis(A_sub, n1, 3, ap3_sub, ap2_sub) if ap3_sub else np.zeros((0,0))
        d3_sub = dim_om(om3_sub)

        if d3_sub > 0:
            embed3 = np.zeros((len(ap3_T_list), d3_sub))
            for j in range(d3_sub):
                a3_sub_vec = om3_sub[:, j]
                for k, path_sub in enumerate(ap3_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap3_T_list:
                        idx_T = ap3_T_list.index(path_T)
                        embed3[idx_T, j] = a3_sub_vec[k]
            chi = np.linalg.lstsq(om3_T, embed3, rcond=None)[0]
        else:
            chi = np.zeros((d3_T, 0))

        rk_chi = np.linalg.matrix_rank(chi, tol=1e-8) if chi.shape[1] > 0 else 0

        # ∂₃^rel maps relative Ω₃ → relative Ω₂
        # In Ω₂(T) coords: im(coords3_T) projected to quotient by im(phi)
        # First get the image of all of Ω₃ in the quotient
        d3_proj = Q.T @ coords3_T  # d_rel × d3_T

        # Restrict to relative Ω₃ elements (those not in T\v)
        if rk_chi > 0:
            U_chi, S_chi, _ = np.linalg.svd(chi, full_matrices=True)
            Q3 = U_chi[:, rk_chi:]
            d3_proj_rel = d3_proj @ Q3
        else:
            d3_proj_rel = d3_proj
    else:
        d3_proj_rel = np.zeros((d_rel, 0))

    rk_d3_rel = np.linalg.matrix_rank(d3_proj_rel, tol=1e-8) if d3_proj_rel.shape[1] > 0 else 0

    h2_rel = z2_rel_dim - rk_d3_rel

    if h2_rel == 0:
        return True, 0  # H₂(T,T\v) = 0, so δ vacuously injective

    # Now check δ-injectivity
    # Find basis for H₂(T,T\v)
    U_rel, S_rel, Vt_rel = np.linalg.svd(coords2_rel_q, full_matrices=True)
    z2_rel_basis = Vt_rel[rk_d2_rel:]  # rows, in relative Ω₂ coords

    # Find which Z₂^rel directions are in im(∂₃^rel)
    if d3_proj_rel.shape[1] > 0:
        proj_z2 = z2_rel_basis @ d3_proj_rel
        rk_proj = np.linalg.matrix_rank(proj_z2, tol=1e-8)
        U_p, S_p, _ = np.linalg.svd(proj_z2, full_matrices=True)
        # H₂^rel basis (uncovered directions)
        h2_rel_basis = z2_rel_basis[rk_proj:] if rk_proj < z2_rel_dim else np.zeros((0, d_rel))
    else:
        h2_rel_basis = z2_rel_basis

    if h2_rel_basis.shape[0] == 0:
        return True, 0

    # For each H₂^rel direction, compute δ:
    # Lift to Ω₂(T) via Q: lift = Q @ h2_dir
    # Apply ∂₂: boundary = coords2_T @ Q @ h2_dir
    # Project to Ω₁(T\v) via psi: δ_val = psi^T @ boundary... wait
    # Actually δ maps to H₁(T\v), so we need ∂₂(lift) projected to Ω₁(T\v)

    # δ(h) = [∂₂(Q @ h)] in H₁(T\v)
    # ∂₂(Q @ h) = coords2_T @ Q @ h (in Ω₁(T) coords)
    # This should be in Ω₁(T\v) since h is a relative cycle
    # (∂₂^rel(h) = 0, meaning ∂₂(Q@h) is in im(psi) up to quotient)
    # Wait, ∂₂^rel(h) = 0 means R^T @ coords2_T @ Q @ h = 0,
    # i.e., coords2_T @ Q @ h ∈ im(psi) (in Ω₁(T\v))

    # So δ(h) = psi^{-1}(coords2_T @ Q @ h) expressed in Ω₁(T\v) coords
    # Then check if this is in Z₁(T\v) (it should be), and if it's in B₁(T\v)

    # Express ∂₂(lift(h)) in Ω₁(T\v) coords
    delta_images = []
    for k in range(h2_rel_basis.shape[0]):
        h = h2_rel_basis[k]
        bd = coords2_T @ Q @ h  # in Ω₁(T) coords

        # Express in Ω₁(T\v) coords
        if d1_sub > 0:
            delta_sub = np.linalg.lstsq(psi, bd, rcond=None)[0]
            err = np.max(np.abs(psi @ delta_sub - bd))
            if err > 1e-6:
                # Not in im(psi)??? This shouldn't happen for a relative cycle
                pass
            delta_images.append(delta_sub)
        else:
            delta_images.append(np.zeros(0))

    if not delta_images or d1_sub == 0:
        return True, h2_rel

    # Check: are the delta images independent in H₁(T\v)?
    # First check they're in Z₁(T\v)
    bd1_sub = build_full_boundary_matrix(ap1_sub, ap0_sub)
    rk_d1_sub = np.linalg.matrix_rank(bd1_sub, tol=1e-8)

    # Z₁(T\v) = ker(∂₁^{T\v})
    # B₁(T\v) = im(∂₂^{T\v})
    if dim_om(om2_sub) > 0:
        bd2_sub = build_full_boundary_matrix(ap2_sub, ap1_sub)
        coords2_sub = np.linalg.lstsq(om1_sub, bd2_sub @ om2_sub, rcond=None)[0]
        rk_d2_sub = np.linalg.matrix_rank(coords2_sub, tol=1e-8)
    else:
        rk_d2_sub = 0

    z1_sub = d1_sub - rk_d1_sub
    beta1_sub = z1_sub - rk_d2_sub

    # Express delta images as elements of Ω₁(T\v) and check their class in H₁
    # H₁(T\v) = Z₁(T\v) / B₁(T\v)
    # We need to check: dim of span of delta images in H₁ = h2_rel

    # Get Z₁ basis
    U1, S1, Vt1 = np.linalg.svd(bd1_sub @ om1_sub, full_matrices=True)
    z1_basis = Vt1[rk_d1_sub:]  # rows, in Ω₁(T\v) coords

    # Project delta images into Z₁
    D = np.array(delta_images)  # h2_rel × d1_sub
    D_z1 = D @ z1_basis.T  # h2_rel × z1_dim

    # Get B₁ basis in Z₁ coords
    if rk_d2_sub > 0:
        B1_om = coords2_sub  # d1_sub × d2_sub
        B1_z1 = z1_basis @ B1_om  # z1_dim × d2_sub (project B₁ into Z₁ coords... wait)
        # Actually B₁ vectors are already in Ω₁ coords. Project to Z₁ coords:
        # We want: im(∂₂) in Z₁ basis.
        # ∂₂(Ω₂) = coords2_sub columns. Express in z1_basis:
        B1_z1 = z1_basis @ coords2_sub  # z1_sub × d2_sub
        rk_B1 = np.linalg.matrix_rank(B1_z1, tol=1e-8)
    else:
        B1_z1 = np.zeros((z1_sub, 0))
        rk_B1 = 0

    # delta images modulo B₁ in Z₁
    if rk_B1 > 0:
        # Combined: [B1_z1 | D_z1^T]
        combined = np.hstack([B1_z1, D_z1.T])
        rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)
        delta_rk_in_H1 = rk_combined - rk_B1
    else:
        delta_rk_in_H1 = np.linalg.matrix_rank(D_z1, tol=1e-8)

    # δ is injective iff delta_rk_in_H1 = h2_rel
    is_injective = (delta_rk_in_H1 == h2_rel)

    return is_injective, h2_rel


print("=" * 70)
print("δ-INJECTIVITY TEST FOR β₂=0 INDUCTION")
print("=" * 70)

for n in [4, 5, 6]:
    m = n*(n-1)//2
    total = 1 << m

    if n <= 6:
        samples = range(total)
    else:
        random.seed(42)
        samples = random.sample(range(total), min(2000, total))

    fail_count = 0
    h2_rel_nonzero = 0
    total_checked = 0
    exists_injective_v = 0
    all_v_injective = 0

    t0 = time.time()
    for idx, bits in enumerate(samples):
        A = build_adj(n, bits)

        # For each vertex v, check δ-injectivity
        any_injective = False
        all_injective = True
        any_h2_rel = False

        for v in range(n):
            inj, h2_rel = check_delta_injective(A, n, v)
            if h2_rel > 0:
                any_h2_rel = True
                if inj:
                    any_injective = True
                else:
                    all_injective = False
                    fail_count += 1
            else:
                # H₂(T,T\v)=0, vacuously injective
                any_injective = True

        if any_h2_rel:
            h2_rel_nonzero += 1

        if any_injective:
            exists_injective_v += 1
        if all_injective and any_h2_rel:
            all_v_injective += 1

        total_checked += 1

        if isinstance(samples, range) and bits % 5000 == 0 and bits > 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {bits}/{total} ({elapsed:.0f}s)")

    elapsed = time.time() - t0
    print(f"\nn={n}: {total_checked} tournaments in {elapsed:.0f}s")
    print(f"  H₂(T,T\\v) > 0 for some v: {h2_rel_nonzero}")
    print(f"  δ injective for ALL v: {all_v_injective}/{h2_rel_nonzero}")
    print(f"  δ injective for SOME v: {exists_injective_v}/{total_checked}")
    print(f"  δ-injectivity failures (individual (T,v) pairs): {fail_count}")

    if fail_count == 0:
        print(f"  ✓ δ is injective for ALL T, ALL v at n={n}")
    else:
        print(f"  ✗ δ fails for {fail_count} (T,v) pairs")


print("\nDone.")
