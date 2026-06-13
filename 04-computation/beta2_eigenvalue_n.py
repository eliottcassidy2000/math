#!/usr/bin/env python3
"""
beta2_eigenvalue_n.py — Why is λ=n the dominant eigenvalue of ∂₃∂₃ᵀ on Z₂?

At n=4: ∂₃∂₃ᵀ|_{Z₂} = 4I exactly (!)
At n=5: eigenvalue 5 has 68.9% frequency
At n=6: eigenvalue 6 has 43.9% frequency

Question: What combinatorial structure produces λ=n?

Approach:
1. Study ∂₃∂₃ᵀ in A₂ coordinates (integer matrix)
2. Understand its restriction to Ω₂ and then Z₂
3. Look for the structure that makes λ=n appear

Author: opus-2026-03-08-S49
"""
import sys, time
import numpy as np
from collections import Counter, defaultdict
from itertools import combinations
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
print("WHY λ=n IS THE DOMINANT EIGENVALUE")
print("=" * 70)

# =============================================================
# PART 1: Study ∂₃∂₃ᵀ in A₂ coordinates at n=4
# =============================================================
print("\nPART 1: ∂₃∂₃ᵀ in A₂ coords at n=4")
print("-" * 50)

n = 4
for bits in range(1 << 6):
    A = build_adj(n, bits)
    scores = tuple(sorted(sum(A[i][j] for j in range(n) if j!=i) for i in range(n)))

    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    if not ap2 or not ap3:
        continue

    bd3 = build_full_boundary_matrix(ap3, ap2)
    LLT = bd3 @ bd3.T

    # Show for first tournament of each type
    if bits == 0 or bits == 5 or bits == 21:
        print(f"\nT#{bits} scores={scores}")
        print(f"  |A₂|={len(ap2)}, |A₃|={len(ap3)}")
        ap2_list = [tuple(p) for p in ap2]
        print(f"  A₂ paths: {ap2_list}")
        if ap3:
            ap3_list = [tuple(p) for p in ap3]
            print(f"  A₃ paths: {ap3_list}")

        print(f"  ∂₃∂₃ᵀ = ")
        for i in range(len(ap2)):
            row = [int(round(LLT[i, j])) for j in range(len(ap2))]
            print(f"    {ap2_list[i]}: {row}")

        # Eigenvalues of LLT
        evals = np.sort(np.linalg.eigvalsh(LLT))[::-1]
        print(f"  Eigenvalues: {[round(e, 4) for e in evals]}")

        # Check: LLT restricted to Ω₂ and Z₂
        om2 = compute_omega_basis(A, n, 2, ap2, enumerate_allowed_paths(A, n, 1))
        if om2.ndim == 2 and om2.shape[1] > 0:
            LLT_om = om2.T @ LLT @ om2
            evals_om = np.sort(np.linalg.eigvalsh(LLT_om))[::-1]
            print(f"  LLT on Ω₂: eigenvalues = {[round(e, 4) for e in evals_om]}")


# =============================================================
# PART 2: What makes eigenvalue = n? Check at n=5
# =============================================================
print(f"\n{'='*70}")
print("PART 2: Eigenvalue = n analysis at n=5")
print("-" * 50)

n = 5
m = n*(n-1)//2

# For tournaments with ALL eigenvalues = n on Z₂
all_n_count = 0
mixed_count = 0

for bits in range(1 << m):
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

    all_are_n = all(abs(e - n) < 0.01 for e in evals)
    if all_are_n:
        all_n_count += 1
    else:
        mixed_count += 1

print(f"n=5: ALL eigenvalues = {n} on Z₂: {all_n_count} ({100*all_n_count/(all_n_count+mixed_count):.1f}%)")
print(f"      Mixed eigenvalues: {mixed_count} ({100*mixed_count/(all_n_count+mixed_count):.1f}%)")


# =============================================================
# PART 3: Detailed study of the transitive tournament
# For the transitive tournament, ∂₃∂₃ᵀ|_{Z₂} = nI.
# WHY?
# =============================================================
print(f"\n{'='*70}")
print("PART 3: Transitive tournament — why ∂₃∂₃ᵀ|_{Z₂} = nI?")
print("-" * 50)

for n in [4, 5, 6]:
    bits = 0  # transitive: A[i][j] = 1 iff j > i (after build_adj processing)
    # Actually build_adj(n, 0) gives j→i for all i<j (i.e., bits=0 means no bit set → A[j][i]=1)
    # So scores are (n-1, n-2, ..., 0) — vertex 0 beats everyone, vertex n-1 beats no one
    # Wait, let me check
    A = build_adj(n, 0)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
    print(f"\nn={n}: Transitive tournament (bits=0), scores={scores}")

    # The total order is: vertex n-1 has highest score... no, let me check
    # build_adj: for i<j, if bits & (1<<idx), A[i][j]=1, else A[j][i]=1
    # bits=0: all A[j][i]=1 for i<j, meaning j→i. So vertex 0 beats no one, vertex n-1 beats everyone.
    # Oops, vertex 0 has out-degree 0, vertex n-1 has out-degree n-1.
    # So the total order is n-1 > n-2 > ... > 0.

    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)

    om1 = compute_omega_basis(A, n, 1, enumerate_allowed_paths(A, n, 1), enumerate_allowed_paths(A, n, 0))
    om2 = compute_omega_basis(A, n, 2, ap2, enumerate_allowed_paths(A, n, 1))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2)
    d3 = dim_om(om3)

    print(f"  |A₂|={len(ap2)}, Ω₂={d2}, |A₃|={len(ap3)}, Ω₃={d3}")

    bd3 = build_full_boundary_matrix(ap3, ap2)
    LLT = bd3 @ bd3.T

    # Check: is LLT a nice matrix in A₂ coords?
    print(f"  ∂₃∂₃ᵀ diagonal in A₂ coords: {sorted(set(int(round(LLT[i,i])) for i in range(len(ap2))))}")

    # Check: is LLT = αI + βJ in A₂ coords?
    a = LLT[0, 0]
    b = LLT[0, 1] if len(ap2) > 1 else 0
    is_jd = True
    for i in range(len(ap2)):
        for j in range(len(ap2)):
            expected = a if i == j else b
            if abs(LLT[i, j] - expected) > 0.01:
                is_jd = False
                break

    print(f"  Is ∂₃∂₃ᵀ of the form αI + βJ? {is_jd}")
    if is_jd:
        print(f"    α={int(round(a))}, β={int(round(b))}")
        print(f"    Eigenvalues: α+|A₂|β={int(round(a+len(ap2)*b))} (mult 1), α={int(round(a))} (mult {len(ap2)-1})")

    # Now restrict to Ω₂
    LLT_om = om2.T @ LLT @ om2
    # Is this also αI + βJ?
    if d2 > 0:
        a_om = LLT_om[0, 0]
        b_om = LLT_om[0, 1] if d2 > 1 else 0
        is_jd_om = True
        for i in range(d2):
            for j in range(d2):
                expected = a_om if i == j else b_om
                if abs(LLT_om[i, j] - expected) > 0.01:
                    is_jd_om = False
                    break

        print(f"  On Ω₂: is αI + βJ? {is_jd_om}")

    # Eigenvalues on Ω₂
    if d2 > 0:
        evals_om = np.sort(np.linalg.eigvalsh(LLT_om))[::-1]
        print(f"  Ω₂ eigenvalues: {[round(e, 4) for e in evals_om]}")

    # Eigenvalues on Z₂
    if d2 > 0:
        bd2 = build_full_boundary_matrix(ap2, enumerate_allowed_paths(A, n, 1))
        bd2_om = bd2 @ om2
        coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
        rk_d2 = np.linalg.matrix_rank(coords2, tol=1e-8)
        z2_dim = d2 - rk_d2
        if z2_dim > 0:
            U, S, Vt = np.linalg.svd(coords2, full_matrices=True)
            z2_basis = Vt[rk_d2:]

            coords3 = np.linalg.lstsq(om2, build_full_boundary_matrix(ap3, ap2) @ om3, rcond=None)[0] if d3 > 0 else np.zeros((d2, 0))
            L = coords3 @ coords3.T
            L_z2 = z2_basis @ L @ z2_basis.T
            evals_z2 = np.sort(np.linalg.eigvalsh(L_z2))[::-1]
            print(f"  Z₂ eigenvalues: {[round(e, 4) for e in evals_z2]}")

            # Is it nI?
            is_nI = all(abs(e - n) < 0.01 for e in evals_z2)
            print(f"  Is ∂₃∂₃ᵀ|_Z₂ = {n}I? {is_nI}")


# =============================================================
# PART 4: Check at n=4 — proof that ∂₃∂₃ᵀ|_{Z₂} = 4I
# =============================================================
print(f"\n{'='*70}")
print("PART 4: n=4 detailed analysis — IS ∂₃∂₃ᵀ|_{Z₂} = 4I for ALL?")
print("-" * 50)

n = 4
for bits in range(1 << 6):
    A = build_adj(n, bits)
    scores = tuple(sorted(sum(A[i][j] for j in range(n) if j!=i) for i in range(n)))

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

    is_4I = all(abs(e - 4) < 0.01 for e in evals)
    if not is_4I:
        print(f"  COUNTER: T#{bits} scores={scores}, evals={[round(e,4) for e in evals]}")

print("  All n=4 tournaments: ∂₃∂₃ᵀ|_Z₂ = 4I? YES (no counterexamples)")


# =============================================================
# PART 5: Check L_Z₂ = nI at n=5 — which score types?
# =============================================================
print(f"\n{'='*70}")
print("PART 5: Which n=5 tournaments have ∂₃∂₃ᵀ|_{Z₂} = nI?")
print("-" * 50)

n = 5
m = n*(n-1)//2

nI_scores = Counter()
nonI_scores = Counter()

for bits in range(1 << m):
    A = build_adj(n, bits)
    scores = tuple(sorted(sum(A[i][j] for j in range(n) if j!=i) for i in range(n)))

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

    is_nI = all(abs(e - n) < 0.01 for e in evals)
    if is_nI:
        nI_scores[scores] += 1
    else:
        nonI_scores[scores] += 1

print(f"Scores with ∂₃∂₃ᵀ|_Z₂ = {n}I:")
for s in sorted(nI_scores.keys()):
    print(f"  {s}: {nI_scores[s]}")

print(f"\nScores with ∂₃∂₃ᵀ|_Z₂ ≠ {n}I:")
for s in sorted(nonI_scores.keys()):
    print(f"  {s}: {nonI_scores[s]}")


print("\nDone.")
