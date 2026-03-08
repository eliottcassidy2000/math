#!/usr/bin/env python3
"""
CORRECTED Fourier decomposition for path homology of circulant digraphs.

The previous version (path_homology_fourier.py) had a bug: it computed
dim(Ω_1^(λ)) = |S| (the number of edge types), but the Ω constraint
requires ∂-invariance WITHIN the eigenspace.

KEY INSIGHT from Tang-Yau:
For circulant digraphs, the Ω_p subspace has its OWN Fourier decomposition.
The dimension of Ω_p^(λ) is NOT just the dimension of A_p^(λ).

The ∂-invariance condition ∂u ∈ A_{p-1} imposes constraints that vary
by eigenvalue λ.

APPROACH: Build the full chain complex within each eigenspace correctly.
"""
import numpy as np
from itertools import combinations, product
from collections import Counter
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import path_betti_numbers, circulant_digraph

def enumerate_step_sequences(S, p):
    """
    Enumerate all allowed p-paths as step sequences (s_1,...,s_p) where
    s_i ∈ S and all partial sums are distinct (mod n) and nonzero.

    For a circulant C_n^S, a p-path from vertex v is determined by
    steps s_1,...,s_p where:
    - Each s_i ∈ S
    - All partial sums σ_0=0, σ_1=s_1, σ_2=s_1+s_2, ... are distinct mod n
    """
    S_list = list(S)
    if p == 0:
        return [()]  # single vertex, no steps

    results = []
    def backtrack(steps, partial_sums):
        if len(steps) == p:
            results.append(tuple(steps))
            return
        for s in S_list:
            new_sum = (partial_sums[-1] + s) % n_global
            if new_sum not in partial_sums and new_sum != 0:
                # Ensure distinctness: new vertex ≠ any previous
                # Actually, partial_sums tracks cumulative offsets
                # We need new_sum ∉ partial_sums (distinct vertices)
                backtrack(steps + [s], partial_sums + [new_sum])

    backtrack([], [0])
    return results

def boundary_step_faces(steps, n):
    """
    Compute the faces of a p-path given by step sequence (s_1,...,s_p).

    The path visits vertices 0, σ_1, σ_2, ..., σ_p (cumulative sums).
    Face i (deleting vertex i) gives a (p-1)-path.

    Returns list of (face_steps, sign, offset) where:
    - face_steps is the step sequence of the face
    - sign is (-1)^i
    - offset is the starting vertex offset (face may start from σ_i)
    """
    p = len(steps)
    partial_sums = [0]
    for s in steps:
        partial_sums.append((partial_sums[-1] + s) % n)

    faces = []
    for i in range(p + 1):
        sign = (-1) ** i
        # Delete vertex at position i (which is at offset partial_sums[i])
        if i == 0:
            # Delete first vertex: new path starts at σ_1
            new_steps = tuple(steps[1:])
            offset = partial_sums[1]
        elif i == p:
            # Delete last vertex
            new_steps = tuple(steps[:-1])
            offset = 0
        else:
            # Delete middle vertex: merge steps i and i+1
            merged = (steps[i-1] + steps[i]) % n
            new_steps = tuple(steps[:i-1]) + (merged,) + tuple(steps[i+1:])
            offset = 0
        faces.append((new_steps, sign, offset))

    return faces

def fourier_betti_correct(S_set, n, max_dim=4):
    """
    Compute Betti numbers using the correct Fourier decomposition.

    For each eigenvalue λ (nth root of unity):
    1. Enumerate step sequences for A_p (allowed p-paths mod shift)
    2. Compute Ω_p^(λ) (∂-invariant subspace in eigenspace λ)
    3. Compute ∂_p on Ω_p^(λ)
    4. Sum β_p^(λ) over all λ
    """
    global n_global
    n_global = n

    roots = [np.exp(2j * np.pi * k / n) for k in range(n)]

    # Pre-compute step sequences for each dimension
    step_seqs = {}
    for p in range(max_dim + 2):  # need p+1 for im(∂_{p+1})
        step_seqs[p] = enumerate_step_sequences(S_set, p)

    betti = [0] * (max_dim + 1)

    for k in range(n):
        lam = roots[k]

        # For each dimension, compute Ω_p^(λ) and ∂_p^(λ)
        omega_dims = []
        boundary_ranks = []

        for p in range(max_dim + 2):
            seqs_p = step_seqs[p]
            if p == 0:
                seqs_pm1 = []
            else:
                seqs_pm1 = step_seqs[p - 1]

            if not seqs_p:
                omega_dims.append(0)
                boundary_ranks.append(0)
                continue

            if p == 0:
                # Ω_0 = A_0 always (0-paths are trivially ∂-invariant)
                omega_dims.append(1)  # one generator per eigenspace
                boundary_ranks.append(0)
                continue

            # For each step sequence, check if ALL its boundary faces
            # correspond to allowed (p-1)-paths (step sequences in seqs_pm1)
            # When checking faces: a face with offset σ contributes a factor λ^σ

            # Build index for (p-1)-step sequences
            seq_pm1_set = set(seqs_pm1)

            # Find ∂-invariant p-step sequences
            omega_seqs = []
            for seq in seqs_p:
                faces = boundary_step_faces(seq, n)
                all_faces_allowed = True
                for face_steps, sign, offset in faces:
                    if face_steps not in seq_pm1_set:
                        all_faces_allowed = False
                        break
                if all_faces_allowed:
                    omega_seqs.append(seq)

            dim_omega = len(omega_seqs)
            omega_dims.append(dim_omega)

            if dim_omega == 0 or not seqs_pm1:
                boundary_ranks.append(0)
                continue

            # Build boundary matrix ∂_p: Ω_p^(λ) → Ω_{p-1}^(λ)
            # But we need Ω_{p-1} not just A_{p-1}
            # For now, map into A_{p-1} and we'll figure out Ω later

            # Actually, we need to map into Ω_{p-1}. Let me find Ω_{p-1} first.
            if p - 1 == 0:
                omega_pm1_seqs = [()]  # Ω_0 has one generator
            else:
                seqs_pm2 = step_seqs[p - 2]
                seq_pm2_set = set(seqs_pm2)
                omega_pm1_seqs = []
                for seq in seqs_pm1:
                    faces = boundary_step_faces(seq, n)
                    all_ok = True
                    for face_steps, sign, offset in faces:
                        if face_steps not in seq_pm2_set:
                            all_ok = False
                            break
                    if all_ok:
                        omega_pm1_seqs.append(seq)

            omega_pm1_idx = {seq: i for i, seq in enumerate(omega_pm1_seqs)}
            dim_omega_pm1 = len(omega_pm1_seqs)

            if dim_omega_pm1 == 0:
                boundary_ranks.append(0)
                continue

            # Build ∂_p matrix in eigenspace λ
            bd_matrix = np.zeros((dim_omega_pm1, dim_omega), dtype=complex)
            for j, seq in enumerate(omega_seqs):
                faces = boundary_step_faces(seq, n)
                for face_steps, sign, offset in faces:
                    if face_steps in omega_pm1_idx:
                        i = omega_pm1_idx[face_steps]
                        # The face starting at offset contributes λ^offset
                        bd_matrix[i, j] += sign * (lam ** offset)

            rank = np.linalg.matrix_rank(bd_matrix, tol=1e-8)
            boundary_ranks.append(rank)

        # Compute β_p^(λ) = dim(Ω_p^(λ)) - rank(∂_p^(λ)) - rank(∂_{p+1}^(λ))
        for p in range(max_dim + 1):
            ker_p = omega_dims[p] - boundary_ranks[p]  # kernel of ∂_p
            im_pp1 = boundary_ranks[p + 1] if p + 1 < len(boundary_ranks) else 0
            beta_p_lam = ker_p - im_pp1
            if beta_p_lam > 0.5:  # should be integer
                betti[p] += round(beta_p_lam)

    return betti

# ===== TEST: Compare corrected Fourier vs full computation =====
print("=" * 70)
print("CORRECTED FOURIER vs FULL COMPUTATION")
print("=" * 70)

mismatches = 0
total = 0
for n_global in [4, 5, 6, 7]:
    n = n_global
    print(f"\nn={n}:", flush=True)
    for size in range(1, min(n, 4)):
        for S in combinations(range(1, n), size):
            S_set = set(S)

            # Full computation
            A = circulant_digraph(n, list(S))
            betti_full = path_betti_numbers(A, n, max_dim=min(n-1, 4))

            # Fourier computation
            betti_fourier = fourier_betti_correct(S_set, n, max_dim=min(n-1, 4))

            total += 1
            match = all(betti_full[p] == betti_fourier[p] for p in range(min(len(betti_full), len(betti_fourier))))
            if not match:
                mismatches += 1
                print(f"  MISMATCH C_{n}^{S_set}: full={betti_full}, fourier={betti_fourier}")
            elif any(b > 0 for b in betti_full[2:]):
                print(f"  MATCH C_{n}^{S_set}: β={betti_full} ✓")

    print(f"  Checked {total} so far, {mismatches} mismatches", flush=True)

print(f"\n\nTotal: {total} tested, {mismatches} mismatches")

# ===== If Fourier works, test on larger n =====
if mismatches == 0:
    print("\n\nFourier decomposition CORRECT! Testing larger n...")

    for n in [9, 11, 13]:
        n_global = n
        print(f"\nn={n}:")
        for size in range(1, 4):
            for S in combinations(range(1, n), size):
                S_set = set(S)
                betti = fourier_betti_correct(S_set, n, max_dim=min(n-1, 5))
                if any(b > 0 for b in betti[2:]):
                    print(f"  C_{n}^{S_set}: β={betti}")

print("\nDone.")
