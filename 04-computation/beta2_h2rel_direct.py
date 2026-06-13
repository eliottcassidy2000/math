#!/usr/bin/env python3
"""
beta2_h2rel_direct.py — DIRECT computation of H₂(T,T\\v) for specific cases

We have an apparent contradiction:
- LES says H₂(T,T\\v) ≅ ker(i_*: H₁(T\\v) → H₁(T))
- β₁(T\\v) = 1, β₁(T) = 0 → ker should be 1-dim
- Earlier scripts claimed H₂(T,T\\v) = 0

Let me compute H₂^rel DIRECTLY from the quotient chain complex.

Author: opus-2026-03-08-S44
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

print("=" * 70)
print("DIRECT H₂^rel COMPUTATION")
print("=" * 70)

# Specific case: T with scores [1,1,2,2,4], v=4
# T: vertex 4 beats all others
n = 5
A = [[0, 0, 0, 1, 0],
     [1, 0, 0, 0, 0],
     [1, 1, 0, 0, 0],
     [0, 1, 1, 0, 0],
     [1, 1, 1, 1, 0]]
v = 4
B = delete_vertex(A, n, v)

print(f"\nT adjacency (v=4 beats everyone):")
for i in range(n):
    arcs = [j for j in range(n) if A[i][j]]
    print(f"  {i} → {arcs}")

print(f"\nT\\{v} adjacency:")
for i in range(n-1):
    arcs = [j for j in range(n-1) if B[i][j]]
    print(f"  {i} → {arcs}")

# Compute Ω bases for T and T\\v at all levels
omega_T = {}
omega_Tv = {}
allowed_T = {}
allowed_Tv = {}

for p in range(5):
    allowed_T[p] = enumerate_allowed_paths(A, n, p)
    if p == 0:
        omega_T[p] = np.eye(n)
    elif allowed_T[p]:
        omega_T[p] = compute_omega_basis(A, n, p, allowed_T[p],
                                          allowed_T[p-1])
    else:
        omega_T[p] = np.zeros((0, 0))

    if p <= 3:
        allowed_Tv[p] = enumerate_allowed_paths(B, n-1, p)
        if p == 0:
            omega_Tv[p] = np.eye(n-1)
        elif allowed_Tv[p]:
            omega_Tv[p] = compute_omega_basis(B, n-1, p, allowed_Tv[p],
                                               allowed_Tv[p-1])
        else:
            omega_Tv[p] = np.zeros((0, 0))

for p in range(5):
    dim_T = omega_T[p].shape[1] if omega_T[p].ndim == 2 else 0
    dim_Tv = omega_Tv.get(p, np.zeros((0,0)))
    dim_Tv = dim_Tv.shape[1] if dim_Tv.ndim == 2 and dim_Tv.shape[0] > 0 else 0
    print(f"  Ω_{p}(T) = {dim_T}, Ω_{p}(T\\v) = {dim_Tv}, R_{p} = {dim_T - dim_Tv}")

# Build the quotient complex R_*
# Strategy: embed Ω(T\\v) into Ω(T) using the A-coordinate embedding.
# Then take a complement in Ω(T) to represent R.
# The boundary map on R is induced from the boundary map on Ω(T).

# Step 1: For each p, embed Ω_p(T\\v) into Ω_p(T)
print(f"\n--- Building quotient complex ---")

embeddings = {}
R_bases = {}

for p in range(4):
    if p not in omega_Tv or omega_Tv[p].ndim != 2 or omega_Tv[p].shape[0] == 0:
        R_bases[p] = omega_T[p].copy() if omega_T[p].ndim == 2 else np.zeros((0,0))
        continue

    om_Tv = omega_Tv[p]
    om_T = omega_T[p]
    dim_Tv = om_Tv.shape[1]
    dim_T = om_T.shape[1] if om_T.ndim == 2 else 0

    if dim_T == 0:
        R_bases[p] = np.zeros((0, 0))
        continue

    # Embed: T\\v path indices → T path indices
    allowed_T_list = [tuple(x) for x in allowed_T[p]]
    allowed_Tv_list = [tuple(x) for x in allowed_Tv[p]]
    T_idx = {path: i for i, path in enumerate(allowed_T_list)}

    embed_mat = np.zeros((len(allowed_T_list), len(allowed_Tv_list)))
    for j, path in enumerate(allowed_Tv_list):
        if path in T_idx:
            embed_mat[T_idx[path], j] = 1

    # Ω(T\\v) in T's A-coordinates
    om_Tv_in_T = embed_mat @ om_Tv  # |A_p(T)| × dim_Tv

    # Express in Ω(T) coordinates
    coords_embed, _, _, _ = np.linalg.lstsq(om_T, om_Tv_in_T, rcond=None)
    # coords_embed is dim_T × dim_Tv

    rk = np.linalg.matrix_rank(coords_embed, tol=1e-8)
    print(f"  p={p}: embed rank = {rk} (should be {dim_Tv})")

    # Find a complement: columns of om_T not in span of embedded Ω(T\\v)
    # Use SVD to find the complement
    U, S, Vt = np.linalg.svd(coords_embed, full_matrices=True)
    # The last (dim_T - rk) columns of Vt^T form the complement in Ω(T)
    complement = Vt[rk:].T  # dim_T × (dim_T - rk)

    # R_p basis in Ω(T) coordinates
    R_bases[p] = complement  # dim_T × R_p_dim

    print(f"  p={p}: R_{p} dim = {R_bases[p].shape[1]}")

# Step 2: Build boundary maps on the quotient
# ∂_p on Ω(T): compute in Ω-coordinates
print(f"\n--- Building boundary maps ---")

bd_omega = {}
for p in range(1, 5):
    if omega_T[p].ndim != 2 or omega_T[p].shape[1] == 0:
        bd_omega[p] = np.zeros((0, 0))
        continue
    if omega_T[p-1].ndim != 2 or omega_T[p-1].shape[1] == 0:
        bd_omega[p] = np.zeros((0, omega_T[p].shape[1]))
        continue

    bd = build_full_boundary_matrix(allowed_T[p], allowed_T[p-1])
    bd_om = bd @ omega_T[p]
    coords, _, _, _ = np.linalg.lstsq(omega_T[p-1], bd_om, rcond=None)
    bd_omega[p] = coords  # dim_{p-1} × dim_p

    print(f"  ∂_{p} in Ω-coords: {bd_omega[p].shape}, rank = {np.linalg.matrix_rank(bd_omega[p], tol=1e-8)}")

# Step 3: Induced boundary on quotient
# ∂_p^R: R_p → R_{p-1}
# In Ω-coords: ∂_p maps om_T column j to sum over om_T columns.
# Restricted to R: ∂_p (R_p basis) projected onto R_{p-1} complement.

bd_R = {}
for p in range(1, 4):
    if p not in R_bases or p-1 not in R_bases:
        bd_R[p] = np.zeros((0, 0))
        continue
    if R_bases[p].shape[1] == 0 or R_bases[p-1].shape[1] == 0:
        bd_R[p] = np.zeros((R_bases[p-1].shape[1] if R_bases[p-1].ndim == 2 else 0,
                            R_bases[p].shape[1] if R_bases[p].ndim == 2 else 0))
        continue

    # ∂_p maps Ω_p(T) → Ω_{p-1}(T) (in Ω-coords)
    # R_p columns are in Ω_p(T)-coords
    # Apply ∂_p: bd_omega[p] @ R_bases[p] gives image in Ω_{p-1}(T)-coords
    image = bd_omega[p] @ R_bases[p]  # dim_{p-1} × R_p_dim

    # Project onto R_{p-1} complement
    # R_{p-1} basis: R_bases[p-1] (dim_{p-1} × R_{p-1}_dim)
    # Project image onto R_{p-1}: coords in R_{p-1} basis
    coords_R, _, _, _ = np.linalg.lstsq(R_bases[p-1], image, rcond=None)
    bd_R[p] = coords_R  # R_{p-1}_dim × R_p_dim

    rk_R = np.linalg.matrix_rank(bd_R[p], tol=1e-8)
    print(f"  ∂_{p}^R: {bd_R[p].shape}, rank = {rk_R}")

# Step 4: Compute H₂^R
# H₂^R = ker(∂₂^R) / im(∂₃^R)
print(f"\n--- Computing H₂^R ---")

if 2 in bd_R and bd_R[2].shape[0] > 0 and bd_R[2].shape[1] > 0:
    S2 = np.linalg.svd(bd_R[2], compute_uv=False)
    rk_d2R = int(sum(s > 1e-8 for s in S2))
    ker_d2R = bd_R[2].shape[1] - rk_d2R
else:
    ker_d2R = R_bases[2].shape[1] if 2 in R_bases and R_bases[2].ndim == 2 else 0

if 3 in bd_R and bd_R[3].shape[0] > 0 and bd_R[3].shape[1] > 0:
    S3 = np.linalg.svd(bd_R[3], compute_uv=False)
    rk_d3R = int(sum(s > 1e-8 for s in S3))
else:
    rk_d3R = 0

h2_rel = ker_d2R - rk_d3R
print(f"  ker(∂₂^R) = {ker_d2R}")
print(f"  rk(∂₃^R) = {rk_d3R}")
print(f"  H₂^R = {h2_rel}")

if h2_rel > 0:
    print(f"\n  *** H₂^rel = {h2_rel} > 0! ***")
    print(f"  This is CONSISTENT with the LES: H₂^rel = ker(i_*) = 1")
    print(f"  The EARLIER claim that H₂^rel=0 was WRONG!")
    print(f"  This means H₂(T,T\\v) ≠ 0 in general.")
    print(f"  The inductive approach via vertex deletion does NOT work directly!")
elif h2_rel == 0:
    print(f"\n  H₂^rel = 0 ✓")
    print(f"  But then the LES gives 0 → H₂^rel → H₁(T\\v) → H₁(T)")
    print(f"  0 → 0 → 1 → 0 which is NOT exact!")
    print(f"  Unless i_*: H₁(T\\v) → H₁(T) is injective despite dim mismatch.")
    print(f"  Wait: i_* could be injective even with β₁(T)=0 if im(i_*)=0.")
    print(f"  But β₁(T\\v)=1 and ker(i_*)=0 means i_* maps a 1-dim space injectively")
    print(f"  into a 0-dim space — IMPOSSIBLE!")

# Let me also check ∂₂∘∂₃ = 0 on the quotient
if 2 in bd_R and 3 in bd_R and bd_R[2].shape[1] > 0 and bd_R[3].shape[1] > 0:
    comp = bd_R[2] @ bd_R[3]
    print(f"\n  ∂₂^R ∘ ∂₃^R max entry: {np.max(np.abs(comp)):.2e}")
else:
    print(f"\n  Cannot check ∂² = 0 (empty matrices)")

# ===== ALTERNATIVE: Direct from boundary matrices in A-coordinates =====
print(f"\n{'='*70}")
print("ALTERNATIVE: DIRECT A-COORDINATE RELATIVE COMPLEX")
print("=" * 70)

# The relative complex more carefully:
# Ω_p(T) has basis in R^{|A_p(T)|} coordinates.
# Ω_p(T\\v) embeds into R^{|A_p(T)|} coordinates.
# R_p = Ω_p(T) / image(Ω_p(T\\v))

# Boundary ∂_p on Ω_p(T) is the restriction of the face map.
# The induced map on R is: ∂_p^R([x]) = [∂_p(x)] where [·] = quotient class.

# To compute H₂^R correctly:
# 1. Build ∂₂: Ω₂(T) → Ω₁(T) and ∂₃: Ω₃(T) → Ω₂(T) in Ω-coords
# 2. Express Ω(T\\v) embedding in Ω(T)-coords at each level
# 3. ker(∂₂^R) = {x ∈ R₂ : ∂₂(x) ∈ Ω₁(T\\v)} / Ω₂(T\\v)
#    = {x ∈ Ω₂(T) : ∂₂(x) ∈ Ω₁(T\\v)} / Ω₂(T\\v)
# 4. im(∂₃^R) = (∂₃(Ω₃(T)) + Ω₂(T\\v)) / Ω₂(T\\v)

# This is the CORRECT quotient computation.

# Step 1: ∂₂ and ∂₃ in Ω-coords (already done above)
# ∂₂: dim_Ω₁ × dim_Ω₂ matrix
# ∂₃: dim_Ω₂ × dim_Ω₃ matrix

# Step 2: Embedding in Ω-coords
# embed_1: Ω₁(T\\v) in Ω₁(T)-coords, shape dim_Ω₁(T) × dim_Ω₁(T\\v)
# embed_2: Ω₂(T\\v) in Ω₂(T)-coords, shape dim_Ω₂(T) × dim_Ω₂(T\\v)

# For the specific case: v=4, n=5
# T\\v has vertices {0,1,2,3}

# Compute embedding matrices
embeds_omega = {}
for p in [1, 2, 3]:
    if omega_Tv.get(p) is None or omega_Tv[p].ndim != 2 or omega_Tv[p].shape[0] == 0:
        embeds_omega[p] = np.zeros((omega_T[p].shape[1], 0))
        continue

    allowed_T_list = [tuple(x) for x in allowed_T[p]]
    allowed_Tv_list = [tuple(x) for x in allowed_Tv[p]]
    T_idx_map = {path: i for i, path in enumerate(allowed_T_list)}

    # Map T\\v paths to T path coordinates
    embed_A = np.zeros((len(allowed_T_list), len(allowed_Tv_list)))
    for j, path in enumerate(allowed_Tv_list):
        if path in T_idx_map:
            embed_A[T_idx_map[path], j] = 1

    # T\\v Ω-basis in T's A-coords
    om_Tv_in_T_A = embed_A @ omega_Tv[p]

    # Express in T's Ω-coords
    om_T_p = omega_T[p]
    if om_T_p.ndim == 2 and om_T_p.shape[1] > 0:
        coords, _, _, _ = np.linalg.lstsq(om_T_p, om_Tv_in_T_A, rcond=None)
        embeds_omega[p] = coords  # dim_Ω_p(T) × dim_Ω_p(T\\v)
    else:
        embeds_omega[p] = np.zeros((0, omega_Tv[p].shape[1]))

# Step 3: ker(∂₂^R)
# x ∈ Ω₂(T) (dim_Ω₂(T) vector in Ω-coords)
# ∂₂(x) ∈ Ω₁(T\\v) iff ∂₂(x) = embed_1 @ y for some y
# i.e., bd_omega[2] @ x is in column space of embeds_omega[1]

print(f"\n  Embedding dims:")
for p in [1,2,3]:
    if p in embeds_omega:
        print(f"    embed_{p}: {embeds_omega[p].shape}")

# Compute preimage of Ω₁(T\\v) under ∂₂
# bd_omega[2] x = embeds_omega[1] y
# System: [bd_omega[2] | -embeds_omega[1]] [x; y] = 0

bd2 = bd_omega[2]  # dim_Ω₁ × dim_Ω₂
embed1 = embeds_omega[1]  # dim_Ω₁ × dim_Ω₁(T\\v)

M = np.hstack([bd2, -embed1])
U, S, Vt = np.linalg.svd(M, full_matrices=True)
rk_M = int(sum(s > 1e-8 for s in S))
null_M = Vt[rk_M:].T  # (dim_Ω₂ + dim_Ω₁(T\\v)) × null_dim

# Extract the x part (first dim_Ω₂ rows)
dim_Om2_T = omega_T[2].shape[1]
dim_Om1_Tv = omega_Tv[1].shape[1] if omega_Tv[1].ndim == 2 else 0

x_part = null_M[:dim_Om2_T, :]  # dim_Ω₂ × null_dim
rk_x = np.linalg.matrix_rank(x_part, tol=1e-8)

print(f"\n  Preimage of Ω₁(T\\v) under ∂₂:")
print(f"    System matrix M: {M.shape}, rank = {rk_M}")
print(f"    Null space dim: {null_M.shape[1]}")
print(f"    x-part rank: {rk_x}")

# The elements x ∈ Ω₂(T) with ∂₂(x) ∈ Ω₁(T\\v) form a subspace of dim = rk_x.
# ker(∂₂^R) = this subspace / Ω₂(T\\v)
# dim ker(∂₂^R) = rk_x - dim_Ω₂(T\\v) (assuming embedding is full rank)

embed2 = embeds_omega[2]
dim_Om2_Tv = embed2.shape[1]

# But we need to intersect: some of these x might be in Ω₂(T\\v) already.
# Project x_part onto complement of embed2 in Ω₂(T)
combined = np.hstack([embed2, x_part])
rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)
dim_ker_d2R = rk_combined - np.linalg.matrix_rank(embed2, tol=1e-8)

print(f"    dim Ω₂(T\\v) = {dim_Om2_Tv}")
print(f"    dim ker(∂₂^R) = rk(embed2 ∪ preimage) - rk(embed2) = {rk_combined} - {np.linalg.matrix_rank(embed2, tol=1e-8)} = {dim_ker_d2R}")

# Step 4: im(∂₃^R)
# im(∂₃^R) = (∂₃(Ω₃(T)) + Ω₂(T\\v)) / Ω₂(T\\v)
bd3 = bd_omega[3]  # dim_Ω₂ × dim_Ω₃
im3 = bd3  # columns are images of Ω₃ basis elements

combined_3 = np.hstack([embed2, im3])
rk_combined_3 = np.linalg.matrix_rank(combined_3, tol=1e-8)
dim_im_d3R = rk_combined_3 - np.linalg.matrix_rank(embed2, tol=1e-8)

print(f"    dim im(∂₃^R) = rk(embed2 ∪ im∂₃) - rk(embed2) = {rk_combined_3} - {np.linalg.matrix_rank(embed2, tol=1e-8)} = {dim_im_d3R}")

# H₂^R
# But we need im(∂₃^R) ⊆ ker(∂₂^R), so the actual H₂ is within ker.
# Let me check if im ⊆ ker:

# im(∂₃) ⊆ ker(∂₂) always (since ∂₂∘∂₃=0)
comp = bd2 @ bd3
print(f"    ∂₂∘∂₃ max entry: {np.max(np.abs(comp)):.2e}")

# Now compute H₂^R = ker(∂₂^R) / im(∂₃^R)
# = dim(ker(∂₂^R)) - dim(im(∂₃^R) ∩ relevant space)
# Since im ⊆ ker (as quotient classes), h2 = dim_ker - dim_im (within quotient)

# Actually we need to check: the image of ∂₃^R lands in ker(∂₂^R).
# Since ∂₂∘∂₃=0 in absolute complex, ∂₃(Ω₃) ⊆ ker(∂₂).
# Modding out by Ω₂(T\\v): im(∂₃^R) ⊆ ker(∂₂^R). ✓

# h₂^R = dim_ker_d2R - dim_im_d3R
h2_rel_direct = dim_ker_d2R - dim_im_d3R
print(f"\n  H₂^rel = ker - im = {dim_ker_d2R} - {dim_im_d3R} = {h2_rel_direct}")

if h2_rel_direct > 0:
    print(f"\n  *** H₂(T,T\\v) = {h2_rel_direct} ≠ 0 ***")
    print(f"  The earlier claim of H₂^rel=0 was WRONG!")
    print(f"  Consistent with LES: H₂^rel ≅ ker(i_*) where beta1(Tv) > beta1(T)")
elif h2_rel_direct == 0:
    print(f"\n  H₂(T,T\\v) = 0")
    print(f"  But LES says H₂^rel ≅ ker(i_*: H₁(T\\v) → H₁(T))")
    print(f"  With β₁(T\\v)=1, β₁(T)=0, this should be ≥ 1.")
    print(f"  CONTRADICTION — need to understand where the argument breaks down.")

# ===== CHECK ALL n=5 (T,v) pairs =====
print(f"\n{'='*70}")
print("SYSTEMATIC CHECK: ALL n=5 (T,v) pairs")
print("=" * 70)

def compute_h2_rel(A, n, v):
    """Compute H₂(T,T\\v) directly from quotient complex."""
    B = delete_vertex(A, n, v)

    # Compute Ω bases
    omega_T_loc = {}
    omega_Tv_loc = {}
    allowed_T_loc = {}
    allowed_Tv_loc = {}

    for p in range(5):
        allowed_T_loc[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            omega_T_loc[p] = np.eye(n)
        elif allowed_T_loc[p]:
            omega_T_loc[p] = compute_omega_basis(A, n, p, allowed_T_loc[p], allowed_T_loc[p-1])
        else:
            omega_T_loc[p] = np.zeros((0, 0))

    for p in range(4):
        allowed_Tv_loc[p] = enumerate_allowed_paths(B, n-1, p)
        if p == 0:
            omega_Tv_loc[p] = np.eye(n-1)
        elif allowed_Tv_loc[p]:
            omega_Tv_loc[p] = compute_omega_basis(B, n-1, p, allowed_Tv_loc[p], allowed_Tv_loc[p-1])
        else:
            omega_Tv_loc[p] = np.zeros((0, 0))

    # Get Ω-dimensions
    def dim_om(om):
        return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

    d2T = dim_om(omega_T_loc[2])
    d3T = dim_om(omega_T_loc[3])
    d1Tv = dim_om(omega_Tv_loc[1])
    d2Tv = dim_om(omega_Tv_loc[2])

    if d2T == 0:
        return 0

    # Build embeddings
    embeds = {}
    for p in [1, 2]:
        if dim_om(omega_Tv_loc[p]) == 0:
            embeds[p] = np.zeros((dim_om(omega_T_loc[p]), 0))
            continue
        T_list = [tuple(x) for x in allowed_T_loc[p]]
        Tv_list = [tuple(x) for x in allowed_Tv_loc[p]]
        T_idx_map = {path: i for i, path in enumerate(T_list)}
        embed_A = np.zeros((len(T_list), len(Tv_list)))
        for j, path in enumerate(Tv_list):
            if path in T_idx_map:
                embed_A[T_idx_map[path], j] = 1
        om_in_A = embed_A @ omega_Tv_loc[p]
        coords, _, _, _ = np.linalg.lstsq(omega_T_loc[p], om_in_A, rcond=None)
        embeds[p] = coords

    # Boundary maps in Ω-coords
    bd = {}
    for p in [2, 3]:
        if dim_om(omega_T_loc[p]) == 0:
            bd[p] = np.zeros((dim_om(omega_T_loc[p-1]), 0))
            continue
        bd_A = build_full_boundary_matrix(allowed_T_loc[p], allowed_T_loc[p-1])
        bd_om = bd_A @ omega_T_loc[p]
        coords, _, _, _ = np.linalg.lstsq(omega_T_loc[p-1], bd_om, rcond=None)
        bd[p] = coords

    # ker(∂₂^R): preimage of embed1 under bd2, modulo embed2
    M = np.hstack([bd[2], -embeds[1]]) if embeds[1].shape[1] > 0 else bd[2]
    U, S, Vt = np.linalg.svd(M, full_matrices=True)
    rk_M = int(sum(s > 1e-8 for s in S))
    null_M = Vt[rk_M:].T
    x_part = null_M[:d2T, :]
    rk_x = np.linalg.matrix_rank(x_part, tol=1e-8)

    combined_2 = np.hstack([embeds[2], x_part]) if embeds[2].shape[1] > 0 else x_part
    rk_comb2 = np.linalg.matrix_rank(combined_2, tol=1e-8)
    rk_emb2 = np.linalg.matrix_rank(embeds[2], tol=1e-8) if embeds[2].shape[1] > 0 else 0
    dim_ker = rk_comb2 - rk_emb2

    # im(∂₃^R): ∂₃(Ω₃) modulo embed2
    if d3T > 0:
        combined_3 = np.hstack([embeds[2], bd[3]]) if embeds[2].shape[1] > 0 else bd[3]
        rk_comb3 = np.linalg.matrix_rank(combined_3, tol=1e-8)
        dim_im = rk_comb3 - rk_emb2
    else:
        dim_im = 0

    return dim_ker - dim_im

# Test ALL
from collections import Counter as Ctr
pairs_n5 = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(pairs_n5)
h2_rel_dist = Ctr()
total = 0

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs_n5):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    for v in range(n):
        h2r = compute_h2_rel(A, n, v)
        h2_rel_dist[h2r] += 1
        total += 1

    if bits % 200 == 0 and bits > 0:
        print(f"  ... {bits}/{1<<m}")

print(f"\n  H₂(T,T\\v) distribution at n={n} ({total} pairs):")
for val, count in sorted(h2_rel_dist.items()):
    print(f"    H₂^rel = {val}: {count}")

print("\nDone.")
