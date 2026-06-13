#!/usr/bin/env python3
"""
beta2_laplacian.py — Laplacian approach to β₂=0

Key idea: β₂=0 ⟺ ∂₃∂₃ᵀ is positive definite on Z₂.

The "Hodge Laplacian" at degree 2 is:
  Δ₂ = ∂₂ᵀ∂₂ + ∂₃∂₃ᵀ

By Hodge theory, ker(Δ₂) ≅ H₂. So β₂=0 ⟺ Δ₂ > 0 on Ω₂.

But more useful: restricted to Z₂ = ker(∂₂), the Laplacian reduces to
  Δ₂|_{Z₂} = ∂₃∂₃ᵀ|_{Z₂}

So β₂=0 ⟺ ∂₃∂₃ᵀ|_{Z₂} > 0 ⟺ min eigenvalue of ∂₃∂₃ᵀ on Z₂ > 0.

We found computationally that min singular value of ∂₃ on Z₂ ≥ √2.
This means min eigenvalue of ∂₃∂₃ᵀ on Z₂ ≥ 2.

Question: can we prove this bound? What is the exact distribution?

Author: opus-2026-03-08-S49
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


print("=" * 70)
print("LAPLACIAN ANALYSIS FOR β₂=0")
print("=" * 70)

# =============================================================
# PART 1: Eigenvalue distribution of ∂₃∂₃ᵀ on Z₂
# =============================================================
print("\nPART 1: Eigenvalues of ∂₃∂₃ᵀ restricted to Z₂")
print("-" * 50)

for n in [4, 5, 6]:
    m = n*(n-1)//2
    total = 1 << m

    all_eigenvalues = []
    min_evals = []
    eval_multisets = Counter()  # rounded eigenvalue multiset → count

    t0 = time.time()
    for bits in range(total):
        A = build_adj(n, bits)

        ap0 = enumerate_allowed_paths(A, n, 0)
        ap1 = enumerate_allowed_paths(A, n, 1)
        ap2 = enumerate_allowed_paths(A, n, 2)
        ap3 = enumerate_allowed_paths(A, n, 3)

        om1 = compute_omega_basis(A, n, 1, ap1, ap0)
        om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
        om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

        d2 = dim_om(om2)
        d3 = dim_om(om3)
        if d2 == 0:
            continue

        # ∂₂ in Ω coords
        bd2 = build_full_boundary_matrix(ap2, ap1)
        bd2_om = bd2 @ om2
        coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
        rk_d2 = np.linalg.matrix_rank(coords2, tol=1e-8)
        z2_dim = d2 - rk_d2
        if z2_dim == 0:
            continue

        # Find Z₂ basis (null space of coords2)
        U, S, Vt = np.linalg.svd(coords2, full_matrices=True)
        z2_basis = Vt[rk_d2:]  # rows are Z₂ basis in Ω₂ coords

        # ∂₃ in Ω coords
        if d3 > 0:
            bd3 = build_full_boundary_matrix(ap3, ap2)
            bd3_om = bd3 @ om3
            coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]  # d2 × d3
        else:
            coords3 = np.zeros((d2, 0))

        # ∂₃∂₃ᵀ in Ω₂ coords
        L = coords3 @ coords3.T  # d2 × d2

        # Restrict to Z₂
        L_z2 = z2_basis @ L @ z2_basis.T  # z2_dim × z2_dim

        # Eigenvalues
        evals = np.sort(np.linalg.eigvalsh(L_z2))
        all_eigenvalues.extend(evals.tolist())
        min_evals.append(evals[0])

        # Round for counting
        rounded = tuple(sorted([round(e, 2) for e in evals]))
        eval_multisets[rounded] += 1

    elapsed = time.time() - t0
    print(f"\nn={n}: {len(min_evals)} tournaments with Z₂>0 ({elapsed:.0f}s)")

    if not min_evals:
        continue

    print(f"  Min eigenvalue: min={min(min_evals):.6f}, max={max(min_evals):.6f}")
    print(f"  All eigenvalues: min={min(all_eigenvalues):.6f}, max={max(all_eigenvalues):.6f}")
    print(f"  Mean eigenvalue: {np.mean(all_eigenvalues):.6f}")

    # Distribution of min eigenvalues
    min_eval_rounded = Counter([round(e, 4) for e in min_evals])
    print(f"\n  Min eigenvalue distribution:")
    for val in sorted(min_eval_rounded.keys()):
        cnt = min_eval_rounded[val]
        pct = 100*cnt/len(min_evals)
        print(f"    λ_min = {val:.4f}: {cnt} ({pct:.1f}%)")

    # Full eigenvalue distribution
    all_eval_rounded = Counter([round(e, 4) for e in all_eigenvalues])
    print(f"\n  Full eigenvalue distribution:")
    for val in sorted(all_eval_rounded.keys()):
        cnt = all_eval_rounded[val]
        pct = 100*cnt/len(all_eigenvalues)
        print(f"    λ = {val:.4f}: {cnt} ({pct:.1f}%)")

    # Eigenvalue multisets (for small n)
    if n <= 5:
        print(f"\n  Eigenvalue multisets:")
        for ms in sorted(eval_multisets.keys()):
            cnt = eval_multisets[ms]
            print(f"    {ms}: {cnt}")


# =============================================================
# PART 2: Also compute the FULL Hodge Laplacian Δ₂ = ∂₂ᵀ∂₂ + ∂₃∂₃ᵀ
# =============================================================
print(f"\n{'='*70}")
print("PART 2: Full Hodge Laplacian Δ₂")
print("-" * 50)

for n in [4, 5]:
    m = n*(n-1)//2
    total = 1 << m

    hodge_eigenvalues = []

    for bits in range(total):
        A = build_adj(n, bits)

        ap0 = enumerate_allowed_paths(A, n, 0)
        ap1 = enumerate_allowed_paths(A, n, 1)
        ap2 = enumerate_allowed_paths(A, n, 2)
        ap3 = enumerate_allowed_paths(A, n, 3)

        om1 = compute_omega_basis(A, n, 1, ap1, ap0)
        om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
        om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

        d1 = dim_om(om1)
        d2 = dim_om(om2)
        d3 = dim_om(om3)
        if d2 == 0:
            continue

        # ∂₂ in Ω coords
        bd2 = build_full_boundary_matrix(ap2, ap1)
        bd2_om = bd2 @ om2
        coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]  # d1 × d2

        # ∂₃ in Ω coords
        if d3 > 0:
            bd3 = build_full_boundary_matrix(ap3, ap2)
            bd3_om = bd3 @ om3
            coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]  # d2 × d3
        else:
            coords3 = np.zeros((d2, 0))

        # Hodge Laplacian
        Delta = coords2.T @ coords2 + coords3 @ coords3.T

        evals = np.sort(np.linalg.eigvalsh(Delta))
        hodge_eigenvalues.append(evals)

    print(f"\nn={n}: Hodge Laplacian Δ₂ eigenvalues")
    all_evals = np.concatenate(hodge_eigenvalues)
    min_evals = [e[0] for e in hodge_eigenvalues]
    print(f"  Min eigenvalue: min={min(min_evals):.6f}, max={max(min_evals):.6f}")
    print(f"  All eigenvalues: min={min(all_evals):.6f}, max={max(all_evals):.6f}")

    # β₂=0 means no zero eigenvalue
    near_zero = sum(1 for e in all_evals if abs(e) < 0.01)
    print(f"  Near-zero eigenvalues (|λ|<0.01): {near_zero}")

    # Distribution
    eval_rounded = Counter([round(e, 2) for e in all_evals])
    print(f"\n  Hodge eigenvalue distribution:")
    for val in sorted(eval_rounded.keys())[:30]:
        cnt = eval_rounded[val]
        pct = 100*cnt/len(all_evals)
        print(f"    λ = {val:.2f}: {cnt} ({pct:.1f}%)")


# =============================================================
# PART 3: Can we compute ∂₃∂₃ᵀ combinatorially?
# (∂₃∂₃ᵀ)_{ij} = Σ_k ∂₃_{ik} ∂₃_{jk}
# where i,j index Ω₂ and k indexes Ω₃
# In A-path language:
# (∂₃∂₃ᵀ)_{(a,b,c),(a',b',c')} = #{3-paths σ s.t. (a,b,c) and (a',b',c') are both faces of σ}
# with appropriate signs
# =============================================================
print(f"\n{'='*70}")
print("PART 3: Combinatorial structure of ∂₃∂₃ᵀ")
print("-" * 50)

# At the A₂ level (before Ω projection):
# (∂₃∂₃ᵀ)_{(a,b,c),(a',b',c')} counts how many A₃ paths σ = (v₀,v₁,v₂,v₃)
# have both (a,b,c) and (a',b',c') as faces, with sign product.
# Faces of σ = (v₀,v₁,v₂,v₃): (v₁,v₂,v₃), -(v₀,v₂,v₃), (v₀,v₁,v₃), -(v₀,v₁,v₂)

# For the DIAGONAL entry (a,b,c) ↔ (a,b,c):
# Count A₃ paths having (a,b,c) as a face.
# (a,b,c) appears as:
#   face₀ = (v₁,v₂,v₃) of (v₀,a,b,c): v₀→a, with sign +1
#   face₁ = (v₀,v₂,v₃) of (a,v₁,b,c) but that means v₁=... no.
# Wait, face₁ of (v₀,v₁,v₂,v₃) is -(v₀,v₂,v₃). So (a,b,c)=(v₀,v₂,v₃) means v₀=a, v₂=b, v₃=c.
# Then v₁ is anything with a→v₁, v₁→b. Sign = -1.
# Similarly for other faces.

# Let me just compute the diagonal of ∂₃∂₃ᵀ in A₂ coords.
n = 5
bits = 0  # transitive tournament
A = build_adj(n, bits)

ap2 = enumerate_allowed_paths(A, n, 2)
ap3 = enumerate_allowed_paths(A, n, 3)
bd3 = build_full_boundary_matrix(ap3, ap2)

# ∂₃∂₃ᵀ in A₂ coords
LLT = bd3 @ bd3.T

print(f"\nTransitive tournament (n=5):")
print(f"  |A₂|={len(ap2)}, |A₃|={len(ap3)}")

# Diagonal entries
diag = np.diag(LLT)
diag_dist = Counter([int(round(d)) for d in diag])
print(f"  Diagonal of ∂₃∂₃ᵀ: {dict(sorted(diag_dist.items()))}")

# What determines the diagonal value?
# (∂₃∂₃ᵀ)_{p,p} = #{A₃ paths having p as a face} (with signs²=1)
# For p = (a,b,c), this is:
#   #{x: x→a, (x,a,b,c) ∈ A₃}  (left extension)
# + #{x: a→x→b, (a,x,b,c) ∈ A₃}  (type 1 insertion, but only if (a,x,b,c) is a path)
# Wait, actually the faces of (v₀,v₁,v₂,v₃) are:
#   (v₁,v₂,v₃) with sign +1
#   (v₀,v₂,v₃) with sign -1
#   (v₀,v₁,v₃) with sign +1
#   (v₀,v₁,v₂) with sign -1

# So (a,b,c) appears as:
# face₀: (v₁,v₂,v₃) = (a,b,c) → v₀→a needed, sign +1
# face₁: (v₀,v₂,v₃) = (a,b,c) → v₁ with v₀=a, v₂=b, v₃=c → a→v₁, v₁→b, sign -1
# face₂: (v₀,v₁,v₃) = (a,b,c) → v₂ with v₀=a, v₁=b, v₃=c → b→v₂, v₂→c, sign +1
# face₃: (v₀,v₁,v₂) = (a,b,c) → v₃ with v₀=a, v₁=b, v₂=c → c→v₃, sign -1

# In sign²=1, the diagonal is:
# #{x→a: (x,a,b,c) ∈ A₃} + #{x: a→x→b, (a,x,b,c) ∈ A₃}
# + #{x: b→x→c, (a,b,x,c) ∈ A₃} + #{c→x: (a,b,c,x) ∈ A₃}

# For (x,a,b,c) ∈ A₃: need x→a (and a→b→c already). Number = d⁻(a) minus {b,c if they→a}.
# More precisely: #{x ∉ {a,b,c}: x→a and x→a,a→b,b→c means (x,a,b,c) is a 3-path}
# But also need it in A₃, which means just x→a, a→b, b→c.

# For (a,x,b,c) ∈ A₃: need a→x, x→b, b→c, x ∉ {a,b,c}. Count = #{x: a→x, x→b, x≠c}
# For (a,b,x,c) ∈ A₃: need a→b, b→x, x→c, x ∉ {a,b,c}. Count = #{x: b→x, x→c, x≠a}
# For (a,b,c,x) ∈ A₃: need a→b, b→c, c→x, x ∉ {a,b,c}. Count = d⁺_remaining(c)

# So diagonal = d⁻_S(a) + #{x: a→x, x→b, x≠c} + #{x: b→x, x→c, x≠a} + d⁺_S(c)
# where S = V\{a,b,c}

# Let's verify this
for idx, p2 in enumerate(ap2[:5]):
    a, b, c = p2
    S = [x for x in range(n) if x not in (a,b,c)]

    # Count each type
    left = sum(1 for x in S if A[x][a])  # d⁻_S(a)
    mid1 = sum(1 for x in range(n) if x not in (a,b,c) and A[a][x] and A[x][b])
    mid2 = sum(1 for x in range(n) if x not in (a,b,c) and A[b][x] and A[x][c])
    right = sum(1 for x in S if A[c][x])  # d⁺_S(c)

    total_count = left + mid1 + mid2 + right
    actual = int(round(diag[idx]))

    match = "✓" if total_count == actual else "✗"
    print(f"  ({a},{b},{c}): left={left} mid1={mid1} mid2={mid2} right={right} = {total_count} {match} (actual={actual})")


# =============================================================
# PART 4: Off-diagonal structure of ∂₃∂₃ᵀ
# When do (a,b,c) and (a',b',c') share an A₃ path?
# =============================================================
print(f"\n{'='*70}")
print("PART 4: Off-diagonal of ∂₃∂₃ᵀ")
print("-" * 50)

# Two 2-paths share a 3-path as faces iff they share 3 vertices.
# (a,b,c) and (a',b',c') are both faces of (v₀,v₁,v₂,v₃) iff:
# - They share a common set of 3 from {v₀,v₁,v₂,v₃}
# - But they occupy different face positions

# For a 4-set {v₀,v₁,v₂,v₃}, the 4 faces are:
# f₀=(v₁,v₂,v₃), f₁=(v₀,v₂,v₃), f₂=(v₀,v₁,v₃), f₃=(v₀,v₁,v₂)
# Signs: +1, -1, +1, -1

# So (∂₃∂₃ᵀ)_{fi,fj} = sign_i * sign_j for each A₃ path containing both.

# The off-diagonal entries are ±1 (since each A₃ path contributes at most once per pair).
# But actually two A₂ paths can share multiple A₃ paths!

# Let me check the actual values
offdiag_vals = Counter()
for i in range(len(ap2)):
    for j in range(i+1, len(ap2)):
        val = int(round(LLT[i,j]))
        if val != 0:
            offdiag_vals[val] += 1

print(f"Transitive tournament (n=5): Off-diagonal values of ∂₃∂₃ᵀ:")
print(f"  {dict(sorted(offdiag_vals.items()))}")
print(f"  Nonzero entries: {sum(offdiag_vals.values())}/{len(ap2)*(len(ap2)-1)//2}")

# Now check: what is ∂₃∂₃ᵀ restricted to Z₂?
# Is it always an integer matrix?
print(f"\n  ∂₃∂₃ᵀ restricted to Z₂ for a few tournaments:")

for bits in [0, 5, 21, 42]:
    A = build_adj(n, bits)

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2)
    d3 = dim_om(om3)
    if d2 == 0:
        continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk_d2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk_d2
    if z2_dim == 0:
        continue

    U, S, Vt = np.linalg.svd(coords2, full_matrices=True)
    z2_basis = Vt[rk_d2:]

    if d3 > 0:
        bd3 = build_full_boundary_matrix(ap3, ap2)
        bd3_om = bd3 @ om3
        coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]
    else:
        coords3 = np.zeros((d2, 0))

    L = coords3 @ coords3.T
    L_z2 = z2_basis @ L @ z2_basis.T
    evals = np.sort(np.linalg.eigvalsh(L_z2))

    scores = tuple(sorted(sum(A[i][j] for j in range(n) if j!=i) for i in range(n)))
    print(f"\n  T#{bits} scores={scores}, z₂={z2_dim}")
    print(f"    Eigenvalues: {[round(e, 4) for e in evals]}")
    # Is it integer-valued?
    int_check = np.max(np.abs(L_z2 - np.round(L_z2)))
    print(f"    Integer matrix? max deviation: {int_check:.6f}")


print("\nDone.")
