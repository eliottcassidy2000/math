#!/usr/bin/env python3
"""
PALEY TOURNAMENT PER-EIGENSPACE ANALYSIS

For P_p (p ≡ 3 mod 4 prime), decompose path homology into eigenspaces
of the cyclic shift τ: v → v+1.

Key findings from P_7:
- All 6 non-trivial eigenspaces contribute β_4=1 each
- Ω dims are palindromic: [1,3,6,9,9,6,3]
- Only Paley has χ = p among all circulant tournaments at p=7

Questions:
1. Does this pattern extend to P_11? (predict β_8 = 10?)
2. Are Ω dims always palindromic?
3. Is χ = p always for Paley?
4. Does the Gauss sum control which dimension gets the β?
"""
import numpy as np
import sys, time
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

import os
old_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_fourier_v3 import (
    enumerate_step_sequences, compute_faces
)
sys.stdout = old_stdout

def qr(p):
    return set((a*a) % p for a in range(1, p))

def per_eigenspace_betti(S_set, n, max_dim=6):
    """Compute per-eigenspace Betti numbers for C_n^S."""

    # Pre-compute step sequences
    step_seqs = {}
    for pp in range(max_dim + 2):
        step_seqs[pp] = enumerate_step_sequences(S_set, n, pp)

    seq_sets = {}
    seq_indices = {}
    for pp in range(max_dim + 2):
        seq_sets[pp] = set(step_seqs[pp])
        seq_indices[pp] = {seq: i for i, seq in enumerate(step_seqs[pp])}

    # Pre-compute face data
    face_data = {}
    junk_types = {}
    for pp in range(1, max_dim + 2):
        face_data[pp] = []
        junk_set = set()
        for seq in step_seqs[pp]:
            faces = compute_faces(seq, n)
            classified = []
            for face_steps, sign, shift in faces:
                if pp - 1 == 0:
                    is_valid = (face_steps == ())
                else:
                    is_valid = (face_steps in seq_sets[pp - 1])
                if not is_valid:
                    junk_set.add(face_steps)
                classified.append((face_steps, sign, shift, is_valid))
            face_data[pp].append(classified)
        junk_types[pp] = sorted(junk_set)

    print(f"\nStep sequence counts for n={n}, S={sorted(S_set)}:")
    for pp in range(max_dim + 2):
        print(f"  dim {pp}: |A_{pp}|={len(step_seqs[pp])}, |junk|={len(junk_types.get(pp, []))}")

    roots = [np.exp(2j * np.pi * k / n) for k in range(n)]
    all_betti = {}
    all_omega_dims = {}

    for k in range(n):
        lam = roots[k]

        omega_dims = {}
        boundary_ranks = [0] * (max_dim + 2)

        for pp in range(max_dim + 2):
            dim_Ap = len(step_seqs[pp])
            if pp == 0:
                omega_dims[0] = 1 if dim_Ap > 0 else 0
                continue
            if dim_Ap == 0:
                omega_dims[pp] = 0
                continue

            junk_list = junk_types[pp]
            n_junk = len(junk_list)
            if n_junk == 0:
                omega_dims[pp] = dim_Ap
                continue

            junk_idx = {j: i for i, j in enumerate(junk_list)}
            J = np.zeros((n_junk, dim_Ap), dtype=complex)
            for j_col, classified_faces in enumerate(face_data[pp]):
                for face_steps, sign, shift, is_valid in classified_faces:
                    if not is_valid and face_steps in junk_idx:
                        row = junk_idx[face_steps]
                        J[row, j_col] += sign * (lam ** shift)

            U, S_vals, Vh = np.linalg.svd(J, full_matrices=True)
            null_dim = J.shape[1] - np.sum(S_vals > 1e-8)
            omega_dims[pp] = null_dim

        # Compute boundary maps between Ω spaces
        omega_bases = {}
        for pp in range(max_dim + 2):
            dim_Ap = len(step_seqs[pp])
            if pp == 0:
                omega_bases[0] = np.eye(1, dtype=complex) if dim_Ap > 0 else np.zeros((0,0), dtype=complex)
                continue
            if omega_dims.get(pp, 0) == 0:
                omega_bases[pp] = np.zeros((dim_Ap, 0), dtype=complex)
                continue
            junk_list = junk_types[pp]
            n_junk = len(junk_list)
            if n_junk == 0:
                omega_bases[pp] = np.eye(dim_Ap, dtype=complex)
                continue
            junk_idx = {j: i for i, j in enumerate(junk_list)}
            J = np.zeros((n_junk, dim_Ap), dtype=complex)
            for j_col, classified_faces in enumerate(face_data[pp]):
                for face_steps, sign, shift, is_valid in classified_faces:
                    if not is_valid and face_steps in junk_idx:
                        row = junk_idx[face_steps]
                        J[row, j_col] += sign * (lam ** shift)
            U, S_vals, Vh = np.linalg.svd(J, full_matrices=True)
            null_dim = J.shape[1] - np.sum(S_vals > 1e-8)
            if null_dim > 0:
                omega_bases[pp] = Vh[-null_dim:, :].T.conj()
            else:
                omega_bases[pp] = np.zeros((dim_Ap, 0), dtype=complex)

        for pp in range(1, max_dim + 2):
            if omega_dims.get(pp, 0) == 0 or omega_dims.get(pp-1, 0) == 0:
                boundary_ranks[pp] = 0
                continue
            dim_Ap = len(step_seqs[pp])
            dim_Apm1 = len(step_seqs[pp-1])
            B_full = np.zeros((dim_Apm1, dim_Ap), dtype=complex)
            for j_col, classified_faces in enumerate(face_data[pp]):
                for face_steps, sign, shift, is_valid in classified_faces:
                    if is_valid and face_steps in seq_indices[pp-1]:
                        row = seq_indices[pp-1][face_steps]
                        B_full[row, j_col] += sign * (lam ** shift)
            omega_p = omega_bases[pp]
            omega_pm1 = omega_bases[pp-1]
            boundary_in_A = B_full @ omega_p
            if omega_pm1.shape[1] > 0:
                bd_matrix, _, _, _ = np.linalg.lstsq(omega_pm1, boundary_in_A, rcond=None)
                rank = np.linalg.matrix_rank(bd_matrix, tol=1e-8)
            else:
                rank = 0
            boundary_ranks[pp] = rank

        betti_k = []
        for pp in range(max_dim + 1):
            dim_omega_p = omega_dims.get(pp, 0)
            rank_p = boundary_ranks[pp]
            rank_pp1 = boundary_ranks[pp+1] if pp+1 < len(boundary_ranks) else 0
            ker_p = dim_omega_p - rank_p
            beta_p = ker_p - rank_pp1
            betti_k.append(max(0, round(beta_p)))

        omega_dim_list = [omega_dims.get(i, 0) for i in range(max_dim + 1)]
        all_betti[k] = betti_k
        all_omega_dims[k] = omega_dim_list

    return all_betti, all_omega_dims

# ===== P_7 reference =====
print("=" * 70)
print("PER-EIGENSPACE BETTI NUMBERS FOR PALEY TOURNAMENTS")
print("=" * 70)

p = 7
S = qr(p)
print(f"\n{'='*50}")
print(f"P_{p}: QR = {sorted(S)}")
print(f"{'='*50}")
t0 = time.time()
betti_per_k, omega_per_k = per_eigenspace_betti(S, p, max_dim=p-1)
t1 = time.time()

# Report
total_betti = [0] * p
for k in range(p):
    nonzero = any(betti_per_k[k][i] > 0 for i in range(len(betti_per_k[k])))
    if nonzero:
        print(f"  k={k}: Ω_dims={omega_per_k[k]}, β={betti_per_k[k]}")
    for i in range(len(betti_per_k[k])):
        total_betti[i] += betti_per_k[k][i]
print(f"  TOTAL: β = {total_betti}")
chi = sum((-1)**i * total_betti[i] for i in range(len(total_betti)))
print(f"  χ = {chi}")
print(f"  Time: {t1-t0:.1f}s")

# Check palindromicity of Ω dims
print(f"\n  Palindromic Ω dims check:")
for k in range(p):
    dims = omega_per_k[k]
    is_palindrome = all(dims[i] == dims[-(i+1)] for i in range(len(dims)//2))
    if k <= 2 or not is_palindrome:
        print(f"    k={k}: {dims} palindromic={is_palindrome}")
all_same_omega = all(omega_per_k[k] == omega_per_k[1] for k in range(1, p))
print(f"  All non-trivial eigenspaces have same Ω dims? {all_same_omega}")

# ===== P_11 =====
p = 11
S = qr(p)
print(f"\n{'='*50}")
print(f"P_{p}: QR = {sorted(S)}")
print(f"{'='*50}")
t0 = time.time()
betti_per_k, omega_per_k = per_eigenspace_betti(S, p, max_dim=p-1)
t1 = time.time()

total_betti = [0] * p
for k in range(p):
    nonzero = any(betti_per_k[k][i] > 0 for i in range(len(betti_per_k[k])))
    if nonzero:
        print(f"  k={k}: Ω_dims={omega_per_k[k]}, β={betti_per_k[k]}")
    for i in range(len(betti_per_k[k])):
        total_betti[i] += betti_per_k[k][i]
print(f"  TOTAL: β = {total_betti}")
chi = sum((-1)**i * total_betti[i] for i in range(len(total_betti)))
print(f"  χ = {chi}")
print(f"  Time: {t1-t0:.1f}s")

# Check palindromicity
print(f"\n  Palindromic Ω dims check:")
for k in range(min(p, 5)):
    dims = omega_per_k[k]
    is_palindrome = all(dims[i] == dims[-(i+1)] for i in range(len(dims)//2))
    print(f"    k={k}: {dims} palindromic={is_palindrome}")
all_same_omega = all(omega_per_k[k] == omega_per_k[1] for k in range(1, p))
print(f"  All non-trivial eigenspaces have same Ω dims? {all_same_omega}")

# ===== P_3 reference =====
p = 3
S = qr(p)
print(f"\n{'='*50}")
print(f"P_{p}: QR = {sorted(S)}")
print(f"{'='*50}")
betti_per_k, omega_per_k = per_eigenspace_betti(S, p, max_dim=p-1)
total_betti = [0] * p
for k in range(p):
    print(f"  k={k}: Ω_dims={omega_per_k[k]}, β={betti_per_k[k]}")
    for i in range(len(betti_per_k[k])):
        total_betti[i] += betti_per_k[k][i]
print(f"  TOTAL: β = {total_betti}")

# ===== Summary & Pattern =====
print(f"\n\n{'='*70}")
print("PATTERN SUMMARY")
print("=" * 70)
print("""
For Paley tournament P_p (p ≡ 3 mod 4):
- P_3: β = (1,1,0), χ = 0
- P_7: β = (1,0,0,0,6,0), χ = 7
- P_11: β = (computed above)

Prediction: β_{p-3} = p-1 for p ≡ 3 mod 4, p > 3?
Or more precisely: β_{(p-3)/2 * 2} = p-1?

The non-trivial eigenspace contributing β at dimension d=(p-3)/2*2:
- p=7: d=4, β_4^(k)=1 for k=1,...,6
- p=11: d=?, expect β_d^(k)=1 for k=1,...,10
""")

print("\nDone.")
