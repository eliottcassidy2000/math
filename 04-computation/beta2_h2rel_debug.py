#!/usr/bin/env python3
"""
beta2_h2rel_debug.py — Debug the H₂^rel computation failure

426 cases give H₂^rel=0 when LES demands H₂^rel=1.
The LES argument: β₂(T)=0, β₁(T)=0, β₁(T\v)=1 forces H₂(T,T\v)=1.

Let's trace through one specific case to find the bug.

Author: opus-2026-03-08-S45
"""
import sys
import numpy as np
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


n = 5
pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]

# bits=9 = 0b01001, v=1
A = [[0] * n for _ in range(n)]
for idx, (i, j) in enumerate(pairs):
    if (9 >> idx) & 1:
        A[i][j] = 1
    else:
        A[j][i] = 1

v = 1
B = delete_vertex(A, n, v)

print("T adjacency:")
for i in range(n):
    print(f"  {i} -> {[j for j in range(n) if A[i][j]]}")

print(f"\nT\\{v} adjacency (vertices {[i for i in range(n) if i != v]}):")
for i in range(n - 1):
    print(f"  {i} -> {[j for j in range(n - 1) if B[i][j]]}")

# ===== Step 1: Verify Ω_p at each level =====
print("\n" + "=" * 60)
print("STEP 1: Ω bases")
print("=" * 60)

for p in range(4):
    ap_T = enumerate_allowed_paths(A, n, p)
    ap_Tv = enumerate_allowed_paths(B, n - 1, p)
    print(f"\n  p={p}:")
    print(f"    |A_{p}(T)| = {len(ap_T)}")
    print(f"    |A_{p}(T\\v)| = {len(ap_Tv)}")

    if p == 0:
        print(f"    dim Ω_{p}(T) = {n}")
        print(f"    dim Ω_{p}(T\\v) = {n - 1}")
        continue

    ap_prev_T = enumerate_allowed_paths(A, n, p - 1)
    ap_prev_Tv = enumerate_allowed_paths(B, n - 1, p - 1)

    om_T = compute_omega_basis(A, n, p, ap_T, ap_prev_T) if ap_T else np.zeros((0, 0))
    om_Tv = compute_omega_basis(B, n - 1, p, ap_Tv, ap_prev_Tv) if ap_Tv else np.zeros((0, 0))

    d_T = om_T.shape[1] if om_T.ndim == 2 else 0
    d_Tv = om_Tv.shape[1] if om_Tv.ndim == 2 and om_Tv.shape[0] > 0 else 0

    print(f"    dim Ω_{p}(T) = {d_T}")
    print(f"    dim Ω_{p}(T\\v) = {d_Tv}")

# ===== Step 2: Check embedding Ω₂(T\v) ⊆ Ω₂(T) =====
print("\n" + "=" * 60)
print("STEP 2: Verify embedding Ω₂(T\\v) ⊆ Ω₂(T)")
print("=" * 60)

# Compute all needed bases
omega = {}
allowed = {}
for p in range(5):
    allowed[('T', p)] = enumerate_allowed_paths(A, n, p)
    if p == 0:
        omega[('T', p)] = np.eye(n)
    elif allowed[('T', p)]:
        omega[('T', p)] = compute_omega_basis(A, n, p, allowed[('T', p)],
                                                allowed[('T', p - 1)])
    else:
        omega[('T', p)] = np.zeros((0, 0))

for p in range(4):
    allowed[('Tv', p)] = enumerate_allowed_paths(B, n - 1, p)
    if p == 0:
        omega[('Tv', p)] = np.eye(n - 1)
    elif allowed[('Tv', p)]:
        omega[('Tv', p)] = compute_omega_basis(B, n - 1, p, allowed[('Tv', p)],
                                                 allowed[('Tv', p - 1)])
    else:
        omega[('Tv', p)] = np.zeros((0, 0))


def dim_om(key):
    om = omega[key]
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0


def embed_paths(p):
    """Embed A_p(T\v) into A_p(T) via inclusion map."""
    T_list = [tuple(x) for x in allowed[('T', p)]]
    Tv_list = [tuple(x) for x in allowed[('Tv', p)]]
    T_idx = {path: i for i, path in enumerate(T_list)}

    incl = np.zeros((len(T_list), len(Tv_list)))
    for j, path in enumerate(Tv_list):
        if path in T_idx:
            incl[T_idx[path], j] = 1
        else:
            print(f"    WARNING: path {path} from T\\v NOT in T!")
    return incl


# For p=2: embed Ω₂(T\v) into A₂(T), check if it's in col_span of Ω₂(T)
om_T2 = omega[('T', 2)]
om_Tv2 = omega[('Tv', 2)]
incl2 = embed_paths(2)

emb2_in_A = incl2 @ om_Tv2  # Ω₂(T\v) in A₂(T) coordinates
print(f"  Ω₂(T\\v) embedded in A₂(T): shape {emb2_in_A.shape}")
print(f"  Ω₂(T) basis in A₂(T): shape {om_T2.shape}")

# Check: is emb2_in_A in column span of om_T2?
coords, residuals, _, _ = np.linalg.lstsq(om_T2, emb2_in_A, rcond=None)
reconstruction = om_T2 @ coords
error = np.max(np.abs(emb2_in_A - reconstruction))
print(f"  Embedding error (should be ~0): {error:.2e}")
print(f"  Coordinates of Ω₂(T\\v) in Ω₂(T) basis:")
print(f"  {coords}")

# ===== Step 3: Check that embedded Ω₂(T\v) elements have ∂₂ ∈ Ω₁(T\v) =====
print("\n" + "=" * 60)
print("STEP 3: Verify ∂₂(embedded Ω₂(T\\v)) ⊆ Ω₁(T\\v)")
print("=" * 60)

bd2_A = build_full_boundary_matrix(allowed[('T', 2)], allowed[('T', 1)])
boundary_of_emb2 = bd2_A @ emb2_in_A  # in A₁(T) coordinates
print(f"  ∂₂(emb Ω₂(T\\v)): shape {boundary_of_emb2.shape}")

# Embed Ω₁(T\v) into A₁(T)
incl1 = embed_paths(1)
om_Tv1 = omega[('Tv', 1)]
emb1_in_A = incl1 @ om_Tv1  # Ω₁(T\v) in A₁(T) coordinates
print(f"  Ω₁(T\\v) embedded in A₁(T): shape {emb1_in_A.shape}")

# Check: is boundary_of_emb2 in col_span of emb1_in_A?
coords_b, _, _, _ = np.linalg.lstsq(emb1_in_A, boundary_of_emb2, rcond=None)
recon_b = emb1_in_A @ coords_b
error_b = np.max(np.abs(boundary_of_emb2 - recon_b))
print(f"  ∂₂ landing error: {error_b:.2e}")
if error_b > 1e-8:
    print("  *** ∂₂(embedded Ω₂(T\\v)) NOT in Ω₁(T\\v)! ***")
    print(f"  Boundary vectors:")
    for col in range(boundary_of_emb2.shape[1]):
        print(f"    col {col}: {boundary_of_emb2[:, col]}")
    print(f"  Residual:")
    for col in range(boundary_of_emb2.shape[1]):
        r = boundary_of_emb2[:, col] - recon_b[:, col]
        if np.max(np.abs(r)) > 1e-8:
            print(f"    col {col}: {r}")

# ===== Step 4: The preimage computation =====
print("\n" + "=" * 60)
print("STEP 4: Preimage of Ω₁(T\\v) under ∂₂|_{Ω₂(T)}")
print("=" * 60)

# bd2_om = ∂₂ restricted to Ω₂(T), in A₁(T) coords
bd2_om = bd2_A @ om_T2  # |A₁(T)| × dim_Ω₂(T)
print(f"  ∂₂|_{{Ω₂}}: shape {bd2_om.shape}")

# System: bd2_om @ c = emb1_in_A @ y
# Equivalently: [bd2_om | -emb1_in_A] [c; y] = 0
d2T = dim_om(('T', 2))
d1Tv = emb1_in_A.shape[1]

M = np.hstack([bd2_om, -emb1_in_A])
print(f"  System matrix: {M.shape}")

U, S, Vt = np.linalg.svd(M, full_matrices=True)
tol = 1e-8
rk_M = int(sum(s > tol for s in S))
print(f"  Rank of M: {rk_M}")
print(f"  Singular values: {S}")

null_dim = M.shape[1] - rk_M
print(f"  Null space dim: {null_dim}")

null_space = Vt[rk_M:].T
c_part = null_space[:d2T, :]
y_part = null_space[d2T:, :]

print(f"  c_part shape: {c_part.shape}, rank: {np.linalg.matrix_rank(c_part, tol=1e-8)}")
print(f"  y_part shape: {y_part.shape}, rank: {np.linalg.matrix_rank(y_part, tol=1e-8)}")

# Preimage in A₂(T) coordinates
preimage_A = om_T2 @ c_part
rk_pre = np.linalg.matrix_rank(preimage_A, tol=1e-8)
print(f"  Preimage in A₂: shape {preimage_A.shape}, rank: {rk_pre}")

# Check: is emb2_in_A in column span of preimage_A?
combined = np.hstack([preimage_A, emb2_in_A])
rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)
print(f"  rank(preimage ∪ emb2) = {rk_combined}")
print(f"  rank(preimage) = {rk_pre}")
print(f"  rank(emb2) = {np.linalg.matrix_rank(emb2_in_A, tol=1e-8)}")
print(f"  Is emb2 ⊆ preimage? {rk_combined == rk_pre}")

if rk_combined != rk_pre:
    print("\n  *** BUG FOUND: emb2 NOT in preimage! ***")
    print("  This means some Ω₂(T\\v) elements, embedded in Ω₂(T),")
    print("  do NOT have ∂₂ ∈ Ω₁(T\\v).")
    print("  But we checked this above and got error ~0!")
    print("  So either the check is wrong or the preimage computation is wrong.")

    # Let me verify by direct multiplication
    print("\n  Direct check: ∂₂(emb2) in emb1 column space?")
    for col in range(emb2_in_A.shape[1]):
        x = emb2_in_A[:, col]
        bx = bd2_A @ x
        # Is bx in col_span of emb1_in_A?
        coords_test, _, _, _ = np.linalg.lstsq(emb1_in_A, bx, rcond=None)
        recon_test = emb1_in_A @ coords_test
        err_test = np.max(np.abs(bx - recon_test))
        print(f"    col {col}: ||∂₂(x) - proj|| = {err_test:.2e}")

    # Also check: is emb2_in_A really in Ω₂(T)?
    print("\n  Is emb2_in_A in col_span of om_T2?")
    coords_check, _, _, _ = np.linalg.lstsq(om_T2, emb2_in_A, rcond=None)
    recon_check = om_T2 @ coords_check
    err_check = np.max(np.abs(emb2_in_A - recon_check))
    print(f"    Error: {err_check:.2e}")

    # Extract the specific emb2 columns not in preimage
    for col in range(emb2_in_A.shape[1]):
        x_A = emb2_in_A[:, col]
        # Express x in Ω₂(T) coordinates
        c_x, _, _, _ = np.linalg.lstsq(om_T2, x_A, rcond=None)
        # Check if c_x is in column span of c_part
        if c_part.shape[1] > 0:
            coords_c, _, _, _ = np.linalg.lstsq(c_part, c_x, rcond=None)
            recon_c = c_part @ coords_c
            err_c = np.max(np.abs(c_x - recon_c))
            print(f"    emb2 col {col} in c_part? Error: {err_c:.2e}")
            if err_c > 1e-6:
                print(f"      c_x = {c_x}")
                # What does ∂₂(x) look like?
                bx = bd2_om @ c_x  # in A₁(T) coords
                print(f"      ∂₂(x) in A₁(T) = {bx}")
                # Is this in emb1_in_A col span?
                coords_b, _, _, _ = np.linalg.lstsq(emb1_in_A, bx, rcond=None)
                recon_b = emb1_in_A @ coords_b
                err_b = np.max(np.abs(bx - recon_b))
                print(f"      ||∂₂(x) - proj onto Ω₁(T\\v)|| = {err_b:.2e}")

# ===== Step 5: Alternative - work purely in Ω-coordinates =====
print("\n" + "=" * 60)
print("STEP 5: Alternative computation in Ω-coordinates")
print("=" * 60)

# Express embeddings in Ω(T) coordinates
emb2_coords, _, _, _ = np.linalg.lstsq(om_T2, emb2_in_A, rcond=None)
print(f"  emb2 in Ω₂(T) coords: shape {emb2_coords.shape}")
print(f"  emb2 rank: {np.linalg.matrix_rank(emb2_coords, tol=1e-8)}")

om_T1 = omega[('T', 1)]
emb1_coords, _, _, _ = np.linalg.lstsq(om_T1, emb1_in_A, rcond=None)
print(f"  emb1 in Ω₁(T) coords: shape {emb1_coords.shape}")
print(f"  emb1 rank: {np.linalg.matrix_rank(emb1_coords, tol=1e-8)}")

# ∂₂ in Ω-coordinates: bd_omega = om_T1^{-1} bd2_A om_T2
bd2_omega = np.linalg.lstsq(om_T1, bd2_om, rcond=None)[0]
print(f"  ∂₂ in Ω coords: shape {bd2_omega.shape}")
print(f"  ∂₂ rank: {np.linalg.matrix_rank(bd2_omega, tol=1e-8)}")

# Preimage in Ω-coordinates:
# bd2_omega @ c = emb1_coords @ y
M_omega = np.hstack([bd2_omega, -emb1_coords])
U_o, S_o, Vt_o = np.linalg.svd(M_omega, full_matrices=True)
rk_Mo = int(sum(s > 1e-8 for s in S_o))
null_omega = Vt_o[rk_Mo:].T
c_omega = null_omega[:d2T, :]
print(f"  M_omega rank: {rk_Mo}, null dim: {null_omega.shape[1]}")
print(f"  c_omega rank: {np.linalg.matrix_rank(c_omega, tol=1e-8)}")

# Check: is emb2_coords in c_omega?
comb_omega = np.hstack([c_omega, emb2_coords])
print(f"  rank(preimage ∪ emb2) in Ω-coords: {np.linalg.matrix_rank(comb_omega, tol=1e-8)}")
print(f"  rank(preimage) in Ω-coords: {np.linalg.matrix_rank(c_omega, tol=1e-8)}")

# ===== Step 6: Direct check using rational/exact arithmetic =====
print("\n" + "=" * 60)
print("STEP 6: Direct rational check")
print("=" * 60)

# Let me just directly verify the LES components.
# H₁(T) and H₁(T\v) computation
def compute_homology_groups(A_mat, n_vert, max_p=3):
    """Compute homology groups with explicit cycle/boundary spaces."""
    results = {}
    for p in range(max_p + 1):
        ap = enumerate_allowed_paths(A_mat, n_vert, p)
        ap_prev = enumerate_allowed_paths(A_mat, n_vert, p - 1) if p > 0 else []
        ap_next = enumerate_allowed_paths(A_mat, n_vert, p + 1) if p < max_p else []

        if p == 0:
            om = np.eye(n_vert)
        elif ap:
            om = compute_omega_basis(A_mat, n_vert, p, ap, ap_prev)
        else:
            om = np.zeros((0, 0))

        d = om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

        # ker(∂_p)
        if p == 0 or d == 0:
            ker_dim = d
            bd_rk = 0
        else:
            bd = build_full_boundary_matrix(ap, ap_prev)
            om_prev = np.eye(n_vert) if p == 1 else (
                compute_omega_basis(A_mat, n_vert, p - 1, ap_prev,
                                     enumerate_allowed_paths(A_mat, n_vert, p - 2)) if ap_prev else np.zeros((0, 0)))
            bd_om = bd @ om
            coords_bd = np.linalg.lstsq(om_prev, bd_om, rcond=None)[0]
            S_bd = np.linalg.svd(coords_bd, compute_uv=False)
            bd_rk = int(sum(s > 1e-8 for s in S_bd))
            ker_dim = d - bd_rk

        # im(∂_{p+1})
        if ap_next:
            om_next = compute_omega_basis(A_mat, n_vert, p + 1, ap_next, ap) if ap_next else np.zeros((0, 0))
            d_next = om_next.shape[1] if om_next.ndim == 2 else 0
            if d_next > 0 and d > 0:
                bd_next = build_full_boundary_matrix(ap_next, ap)
                bd_next_om = bd_next @ om_next
                coords_next = np.linalg.lstsq(om, bd_next_om, rcond=None)[0]
                S_next = np.linalg.svd(coords_next, compute_uv=False)
                im_rk = int(sum(s > 1e-8 for s in S_next))
            else:
                im_rk = 0
        else:
            im_rk = 0

        results[p] = {
            'dim': d, 'ker': ker_dim, 'im_next': im_rk,
            'beta': ker_dim - im_rk,
            'om': om, 'ap': ap
        }
    return results

print("Computing H_*(T):")
H_T = compute_homology_groups(A, n, 3)
for p in range(4):
    r = H_T[p]
    print(f"  H_{p}(T): dim Ω={r['dim']}, ker ∂={r['ker']}, im ∂_{{p+1}}={r['im_next']}, β={r['beta']}")

print("\nComputing H_*(T\\v):")
H_Tv = compute_homology_groups(B, n - 1, 3)
for p in range(4):
    r = H_Tv[p]
    print(f"  H_{p}(T\\v): dim Ω={r['dim']}, ker ∂={r['ker']}, im ∂_{{p+1}}={r['im_next']}, β={r['beta']}")

# The LES connecting map ∂_*: H₂(T,T\v) → H₁(T\v) should map surjectively
# onto ker(i_*: H₁(T\v) → H₁(T))
b1_T = H_T[1]['beta']
b1_Tv = H_Tv[1]['beta']
b2_T = H_T[2]['beta']
print(f"\nβ₁(T) = {b1_T}, β₁(T\\v) = {b1_Tv}, β₂(T) = {b2_T}")
print(f"LES: ... → H₂(T) → H₂(T,T\\v) → H₁(T\\v) → H₁(T) → ...")
print(f"         ... → {b2_T} → H₂^rel → {b1_Tv} → {b1_T} → ...")
if b2_T == 0:
    print(f"Since β₂(T)=0, H₂^rel = ker(i_*: H₁(T\\v) → H₁(T))")
    print(f"With β₁(T)=0, i_* is zero map, so H₂^rel = β₁(T\\v) = {b1_Tv}")

# ===== Step 7: Explicitly compute the connecting map =====
print("\n" + "=" * 60)
print("STEP 7: Explicit connecting map ∂_*")
print("=" * 60)

# The connecting map ∂_*: H₂(T,T\v) → H₁(T\v) is the snake lemma map.
# For a cycle z ∈ R₂ (= Ω₂(T)/Ω₂(T\v)):
# 1. Lift z to x ∈ Ω₂(T)
# 2. Compute ∂₂(x) ∈ Ω₁(T)
# 3. Since π(∂₂(x)) = ∂₂^R(z) = 0, we have ∂₂(x) ∈ Ω₁(T\v)
# 4. ∂_*(z) = [∂₂(x)] ∈ H₁(T\v) = Ω₁(T\v)/im(∂₂|_{T\v})

# For this to make sense, we need R₂ cycles.
# But we computed dim(ker ∂₂^R) = 3 and dim(im ∂₃^R) = 3, so H₂^rel = 0.
# This means EVERY R₂ cycle is a R₂ boundary. So there are no non-trivial H₂^rel classes.

# BUT the LES says H₂^rel = 1. So EITHER:
# a) ker ∂₂^R = 3 is correct but im ∂₃^R = 3 is wrong (should be 2), or
# b) ker ∂₂^R = 3 is wrong (should be 4), or
# c) The quotient complex setup is wrong

# Let me compute everything from scratch with INTEGERS

# List all Ω₂(T) basis elements explicitly
ap2_T = [tuple(x) for x in allowed[('T', 2)]]
ap1_T = [tuple(x) for x in allowed[('T', 1)]]
ap3_T = [tuple(x) for x in allowed[('T', 3)]]

print(f"\nAllowed 2-paths of T ({len(ap2_T)}):")
for p in ap2_T:
    print(f"  {p}")

om2_T = omega[('T', 2)]
print(f"\nΩ₂(T) basis ({om2_T.shape[1]} vectors in A₂ coords):")
for col in range(om2_T.shape[1]):
    v = om2_T[:, col]
    nz = [(ap2_T[i], v[i]) for i in range(len(v)) if abs(v[i]) > 1e-10]
    print(f"  e{col}: {nz}")

# Ω₂(T\v)
ap2_Tv = [tuple(x) for x in allowed[('Tv', 2)]]
ap1_Tv = [tuple(x) for x in allowed[('Tv', 1)]]
om2_Tv = omega[('Tv', 2)]
print(f"\nΩ₂(T\\v) basis ({om2_Tv.shape[1] if om2_Tv.ndim == 2 else 0} vectors):")
if om2_Tv.ndim == 2:
    for col in range(om2_Tv.shape[1]):
        v_vec = om2_Tv[:, col]
        nz = [(ap2_Tv[i], v_vec[i]) for i in range(len(v_vec)) if abs(v_vec[i]) > 1e-10]
        print(f"  f{col}: {nz}")

# Boundary map ∂₂ in explicit form
bd2 = build_full_boundary_matrix(ap2_T, ap1_T)
print(f"\n∂₂ matrix (A₁ × A₂): {bd2.shape}")

# ∂₂ restricted to Ω₂(T)
bd2_om = bd2 @ om2_T
print(f"∂₂|_Ω₂(T) (A₁ × dim_Ω₂): {bd2_om.shape}")

# Image of each Ω₂ basis element under ∂₂
print(f"\n∂₂ of each Ω₂(T) basis element:")
for col in range(bd2_om.shape[1]):
    bv = bd2_om[:, col]
    nz = [(ap1_T[i], bv[i]) for i in range(len(bv)) if abs(bv[i]) > 1e-10]
    print(f"  ∂₂(e{col}) = {nz}")

# Which of these land in the embedded Ω₁(T\v)?
print(f"\nEdges of T\\v (= A₁(T\\v)):")
for p in ap1_Tv:
    print(f"  {p}")

print(f"\nEdges through v={v}:")
for i, p in enumerate(ap1_T):
    if v in p:
        print(f"  {p} (index {i})")

# Ω₁ = A₁, so Ω₁(T\v) embedded in A₁(T) = edges NOT through v
print(f"\n∂₂(Ω₂(T)) components through v={v}:")
for col in range(bd2_om.shape[1]):
    bv = bd2_om[:, col]
    through_v = [(ap1_T[i], bv[i]) for i in range(len(bv)) if abs(bv[i]) > 1e-10 and v in ap1_T[i]]
    if through_v:
        print(f"  ∂₂(e{col}) through v: {through_v}")
    else:
        print(f"  ∂₂(e{col}) through v: NONE (fully in T\\v)")

# Elements with ∂₂ ∈ Ω₁(T\v) are those whose boundary has NO component on v-edges
# i.e., the "through v" part is zero
print(f"\n∂₂ in Ω₁(T\\v) iff all edge-v components are zero.")
print(f"Elements of Ω₂(T) with ∂₂ ∈ Ω₁(T\\v):")
preimage_count = 0
for col in range(bd2_om.shape[1]):
    bv = bd2_om[:, col]
    through_v_vals = [bv[i] for i in range(len(bv)) if v in ap1_T[i]]
    if all(abs(x) < 1e-10 for x in through_v_vals):
        print(f"  e{col}: YES")
        preimage_count += 1
    else:
        print(f"  e{col}: NO (v-components: {through_v_vals})")

print(f"\nCount of individual basis elements in preimage: {preimage_count}")
print(f"But LINEAR COMBINATIONS might also work!")

# Build the constraint: v-edge components of ∂₂(x) = 0
v_edge_indices = [i for i in range(len(ap1_T)) if v in ap1_T[i]]
print(f"\nEdge indices through v: {v_edge_indices}")
constraint_rows = bd2_om[v_edge_indices, :]
print(f"Constraint matrix (v-edge rows of ∂₂|_Ω₂): shape {constraint_rows.shape}")
print(f"Constraint matrix:")
print(constraint_rows)

# Null space of constraint = preimage
S_c = np.linalg.svd(constraint_rows, compute_uv=False)
rk_c = int(sum(s > 1e-8 for s in S_c))
preimage_dim = bd2_om.shape[1] - rk_c
print(f"\nRank of constraint: {rk_c}")
print(f"Preimage dimension: {preimage_dim}")
print(f"Ω₂(T\\v) dimension: {dim_om(('Tv', 2))}")
print(f"ker(∂₂^R) = preimage - Ω₂(T\\v) = {preimage_dim} - {dim_om(('Tv', 2))} = {preimage_dim - dim_om(('Tv', 2))}")

# Now check im(∂₃^R)
if ap3_T:
    bd3 = build_full_boundary_matrix(ap3_T, ap2_T)
    om3_T = omega[('T', 3)]
    d3T = om3_T.shape[1] if om3_T.ndim == 2 else 0
    if d3T > 0:
        bd3_om = bd3 @ om3_T  # Ω₃(T) boundary in A₂(T) coords
        # Express in Ω₂(T) coords
        bd3_in_Om2, _, _, _ = np.linalg.lstsq(om2_T, bd3_om, rcond=None)
        print(f"\n∂₃|_Ω₃(T) in Ω₂(T) coords: {bd3_in_Om2.shape}")

        # im(∂₃^R) = (im(∂₃) + Ω₂(T\v)) / Ω₂(T\v) inside ker(∂₂^R)
        # But we computed ker(∂₂^R) as the null space of constraint_rows.
        # Let me compute everything in Ω₂(T) coordinates.

        # Null space of constraint_rows = preimage
        U_c, S_c_full, Vt_c = np.linalg.svd(constraint_rows, full_matrices=True)
        rk_c = int(sum(s > 1e-8 for s in S_c_full))
        preimage_basis = Vt_c[rk_c:].T  # d2T × preimage_dim

        print(f"Preimage basis: {preimage_basis.shape}")

        # emb2 in Ω₂(T) coords
        emb2_om = emb2_coords  # computed earlier
        print(f"emb2 in Ω₂: {emb2_om.shape}, rank {np.linalg.matrix_rank(emb2_om, tol=1e-8)}")

        # Check emb2 ⊆ preimage
        comb_check = np.hstack([preimage_basis, emb2_om])
        print(f"rank(preimage ∪ emb2) = {np.linalg.matrix_rank(comb_check, tol=1e-8)}")
        print(f"rank(preimage) = {np.linalg.matrix_rank(preimage_basis, tol=1e-8)}")

        # ker(∂₂^R) = preimage / emb2
        dim_ker_R = preimage_dim - dim_om(('Tv', 2))

        # im(∂₃) in Ω₂ coords = bd3_in_Om2 columns
        # im(∂₃^R) = (span(bd3_in_Om2) + span(emb2_om)) / span(emb2_om) restricted to preimage
        im3_and_emb2 = np.hstack([emb2_om, bd3_in_Om2])
        rk_im3_emb2 = np.linalg.matrix_rank(im3_and_emb2, tol=1e-8)
        rk_emb2_only = np.linalg.matrix_rank(emb2_om, tol=1e-8)
        dim_im_R = rk_im3_emb2 - rk_emb2_only
        print(f"\nim(∂₃^R) = rank(im∂₃ ∪ emb2) - rank(emb2) = {rk_im3_emb2} - {rk_emb2_only} = {dim_im_R}")

        # But we need im(∂₃^R) WITHIN ker(∂₂^R)
        # Since ∂₂∘∂₃=0, im(∂₃) ⊆ ker(∂₂) (absolute).
        # Do all ∂₃ images land in the preimage (i.e., constraint_rows @ bd3_in_Om2 = 0)?
        test = constraint_rows @ bd3_in_Om2
        print(f"∂₃ images satisfy v-constraint? max entry: {np.max(np.abs(test)):.2e}")
        # They should, since ∂₂∘∂₃ = 0 implies ALL components of ∂₂(∂₃(x)) are zero,
        # including the v-edge components.

        h2_rel = dim_ker_R - dim_im_R
        print(f"\nH₂^rel = {dim_ker_R} - {dim_im_R} = {h2_rel}")

print("\nDone.")
