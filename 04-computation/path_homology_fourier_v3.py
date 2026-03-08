#!/usr/bin/env python3
"""
CORRECT Fourier decomposition for path homology of circulant digraphs (v3).

KEY FIX from v1/v2:
Ω_p is NOT span{individually ∂-invariant step sequences}.
Ω_p = {c ∈ A_p : ∂c ∈ A_{p-1}} — includes COMBINATIONS where
non-allowed faces cancel between different paths.

Algorithm for each eigenspace λ:
1. Enumerate A_p step sequences (allowed p-paths up to shift)
2. Compute boundary, separating "valid" faces (in A_{p-1}) from "junk"
3. Ω_p^(λ) = ker(junk_matrix(λ))
4. ∂_p: Ω_p^(λ) → Ω_{p-1}^(λ) via valid face contributions
5. β_p = dim(ker ∂_p / im ∂_{p+1})

This correctly handles the C_4^{1,2} case where v2 failed.
"""
import numpy as np
from itertools import combinations
from collections import Counter
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import path_betti_numbers, circulant_digraph


def enumerate_step_sequences(S_set, n, p):
    """
    Enumerate all valid p-step sequences from connection set S.
    A step sequence (s_1,...,s_p) with s_i ∈ S is valid if
    all partial sums 0, s_1, s_1+s_2, ... are distinct mod n and nonzero
    (except the initial 0).
    """
    S_list = sorted(S_set)
    if p == 0:
        return [()]

    results = []
    def backtrack(steps, partial_sums_set, current_sum):
        if len(steps) == p:
            results.append(tuple(steps))
            return
        for s in S_list:
            new_sum = (current_sum + s) % n
            if new_sum != 0 and new_sum not in partial_sums_set:
                backtrack(steps + [s],
                         partial_sums_set | {new_sum},
                         new_sum)

    backtrack([], {0}, 0)
    return results


def compute_faces(steps, n):
    """
    Compute all faces of a p-path given by step sequence.
    Returns list of (face_steps, sign, shift) where:
    - face_steps: tuple of steps for the face
    - sign: (-1)^i
    - shift: the starting vertex offset (λ^shift factor in eigenspace)
    """
    p = len(steps)
    # Compute partial sums
    partial_sums = [0]
    for s in steps:
        partial_sums.append((partial_sums[-1] + s) % n)

    faces = []
    for i in range(p + 1):
        sign = (-1) ** i
        if i == 0:
            # Delete first vertex: remaining starts at σ_1
            face_steps = tuple(steps[1:])
            shift = partial_sums[1]
        elif i == p:
            # Delete last vertex
            face_steps = tuple(steps[:-1])
            shift = 0
        else:
            # Delete middle vertex i: merge steps i-1 and i
            merged = (steps[i-1] + steps[i]) % n
            face_steps = tuple(steps[:i-1]) + (merged,) + tuple(steps[i+1:])
            shift = 0
        faces.append((face_steps, sign, shift))

    return faces


def fourier_betti_v3(S_set, n, max_dim=4):
    """
    Compute Betti numbers using correct Fourier decomposition.
    """
    roots = [np.exp(2j * np.pi * k / n) for k in range(n)]

    # Pre-compute step sequences for each dimension
    step_seqs = {}
    for p in range(max_dim + 2):
        step_seqs[p] = enumerate_step_sequences(S_set, n, p)

    # Make index sets
    seq_sets = {}
    seq_indices = {}
    for p in range(max_dim + 2):
        seq_sets[p] = set(step_seqs[p])
        seq_indices[p] = {seq: i for i, seq in enumerate(step_seqs[p])}

    # Pre-compute faces and classify them
    # For each p and each sequence, compute faces and note which are valid/junk
    face_data = {}
    junk_types = {}  # p -> set of junk face step sequences

    for p in range(1, max_dim + 2):
        face_data[p] = []
        junk_set = set()
        for seq in step_seqs[p]:
            faces = compute_faces(seq, n)
            classified = []
            for face_steps, sign, shift in faces:
                if p - 1 == 0:
                    # 0-step sequences: always valid (just a vertex)
                    is_valid = (face_steps == ())
                else:
                    is_valid = (face_steps in seq_sets[p - 1])
                if not is_valid:
                    junk_set.add(face_steps)
                classified.append((face_steps, sign, shift, is_valid))
            face_data[p].append(classified)
        junk_types[p] = sorted(junk_set)

    betti = [0] * (max_dim + 1)

    for k in range(n):
        lam = roots[k]

        # For each dimension, compute Ω_p^(λ) and boundary maps
        omega_bases = {}  # p -> matrix whose columns span Ω_p^(λ)
        omega_dims = {}

        for p in range(max_dim + 2):
            dim_Ap = len(step_seqs[p])

            if p == 0:
                # Ω_0 = A_0 always (no boundary to check)
                omega_bases[0] = np.eye(1, dtype=complex) if dim_Ap > 0 else np.zeros((0, 0), dtype=complex)
                omega_dims[0] = 1 if dim_Ap > 0 else 0
                continue

            if dim_Ap == 0:
                omega_bases[p] = np.zeros((0, 0), dtype=complex)
                omega_dims[p] = 0
                continue

            # Build junk matrix J_p(λ): constraints on A_p coefficients
            # For each junk face type, sum contributions from all A_p generators
            junk_list = junk_types[p]
            n_junk = len(junk_list)
            junk_idx = {j: i for i, j in enumerate(junk_list)}

            if n_junk == 0:
                # No junk constraints: Ω_p = A_p
                omega_bases[p] = np.eye(dim_Ap, dtype=complex)
                omega_dims[p] = dim_Ap
                continue

            # Build junk matrix: rows = junk types, cols = A_p generators
            J = np.zeros((n_junk, dim_Ap), dtype=complex)
            for j_col, classified_faces in enumerate(face_data[p]):
                for face_steps, sign, shift, is_valid in classified_faces:
                    if not is_valid and face_steps in junk_idx:
                        row = junk_idx[face_steps]
                        J[row, j_col] += sign * (lam ** shift)

            # Ω_p^(λ) = ker(J)
            # Use SVD to find null space
            if J.shape[0] > 0 and J.shape[1] > 0:
                U, S_vals, Vh = np.linalg.svd(J, full_matrices=True)
                tol = 1e-8
                null_dim = np.sum(S_vals < tol) if len(S_vals) > 0 else J.shape[1]
                # Also account for more columns than rows
                null_dim = J.shape[1] - np.sum(S_vals > tol)
                if null_dim > 0:
                    omega_basis = Vh[-null_dim:, :].T.conj()  # columns = null vectors
                else:
                    omega_basis = np.zeros((dim_Ap, 0), dtype=complex)
            else:
                omega_basis = np.eye(dim_Ap, dtype=complex)
                null_dim = dim_Ap

            omega_bases[p] = omega_basis
            omega_dims[p] = omega_basis.shape[1]

        # Now build boundary matrices ∂_p: Ω_p^(λ) → Ω_{p-1}^(λ)
        boundary_ranks = [0] * (max_dim + 2)

        for p in range(1, max_dim + 2):
            if omega_dims.get(p, 0) == 0 or omega_dims.get(p-1, 0) == 0:
                boundary_ranks[p] = 0
                continue

            dim_Ap = len(step_seqs[p])
            dim_Apm1 = len(step_seqs[p - 1])

            # Build full boundary in A-space: ∂_p^full: A_p → A_{p-1}
            # Only include valid faces (junk is zero by Ω constraint)
            B_full = np.zeros((dim_Apm1, dim_Ap), dtype=complex)
            for j_col, classified_faces in enumerate(face_data[p]):
                for face_steps, sign, shift, is_valid in classified_faces:
                    if is_valid and face_steps in seq_indices[p - 1]:
                        row = seq_indices[p - 1][face_steps]
                        B_full[row, j_col] += sign * (lam ** shift)

            # Restrict to Ω_p → Ω_{p-1}:
            # Input: Ω_p basis (columns in A_p space)
            # Output: project onto Ω_{p-1} basis

            omega_p = omega_bases[p]      # dim_Ap × dim_Ωp
            omega_pm1 = omega_bases[p-1]  # dim_Apm1 × dim_Ωpm1

            # ∂ in A-space: B_full @ omega_p  (dim_Apm1 × dim_Ωp)
            boundary_in_A = B_full @ omega_p

            # Project onto Ω_{p-1}: solve omega_pm1 @ x = boundary_in_A
            # Using least squares (should be exact up to numerics)
            if omega_pm1.shape[1] > 0:
                bd_matrix, residuals, _, _ = np.linalg.lstsq(omega_pm1, boundary_in_A, rcond=None)
                rank = np.linalg.matrix_rank(bd_matrix, tol=1e-8)
            else:
                rank = 0

            boundary_ranks[p] = rank

        # Compute β_p^(λ)
        for p in range(max_dim + 1):
            dim_omega_p = omega_dims.get(p, 0)
            rank_p = boundary_ranks[p]     # rank of ∂_p
            rank_pp1 = boundary_ranks[p + 1] if p + 1 < len(boundary_ranks) else 0

            ker_p = dim_omega_p - rank_p
            beta_p_lam = ker_p - rank_pp1
            if beta_p_lam > 0.5:
                betti[p] += round(beta_p_lam)

    return betti


# ===== TEST: Compare v3 Fourier vs full computation =====
print("=" * 70)
print("FOURIER v3 (CORRECT Ω) vs FULL COMPUTATION")
print("=" * 70)

mismatches = 0
total = 0

for n in [3, 4, 5, 6, 7]:
    print(f"\nn={n}:", flush=True)
    for size in range(1, min(n, 4)):
        for S in combinations(range(1, n), size):
            S_set = set(S)

            # Full computation
            A = circulant_digraph(n, list(S))
            max_d = min(n - 1, 5)
            betti_full = path_betti_numbers(A, n, max_dim=max_d)

            # Fourier computation
            betti_fourier = fourier_betti_v3(S_set, n, max_dim=max_d)

            total += 1
            match = all(betti_full[p] == betti_fourier[p]
                       for p in range(min(len(betti_full), len(betti_fourier))))
            if not match:
                mismatches += 1
                print(f"  MISMATCH C_{n}^{S_set}: full={betti_full}, fourier={betti_fourier}")
            elif any(b > 0 for b in betti_full[1:]):
                print(f"  ✓ C_{n}^{S_set}: β={betti_full}")

    print(f"  Checked {total} so far, {mismatches} mismatches", flush=True)

print(f"\n\nTotal: {total} tested, {mismatches} mismatches")

# ===== If correct, test large n =====
if mismatches == 0:
    print("\n\n" + "=" * 70)
    print("FOURIER v3 VALIDATED! Testing large circulants...")
    print("=" * 70)

    for n in [11, 13, 17, 19, 23]:
        print(f"\nn={n}:", flush=True)
        for size in range(1, 4):
            for S in combinations(range(1, n), size):
                S_set = set(S)
                betti = fourier_betti_v3(S_set, n, max_dim=min(n-1, 5))
                if any(b > 0 for b in betti[2:]):
                    print(f"  C_{n}^{S_set}: β={betti}")

    # Test very large n for specific connection sets
    print("\n\n" + "=" * 70)
    print("LARGE N STABILITY TESTS")
    print("=" * 70)

    for S_list, label in [({1,2}, "S={1,2}"), ({1,3}, "S={1,3}"),
                          ({1,2,3}, "S={1,2,3}"), ({1,2,4}, "S={1,2,4}")]:
        print(f"\n{label}:")
        for n in [29, 37, 41, 47, 53, 59, 67, 71]:
            if max(S_list) < n:
                betti = fourier_betti_v3(S_list, n, max_dim=min(n-1, 6))
                print(f"  n={n}: β={betti}")

print("\nDone.")
