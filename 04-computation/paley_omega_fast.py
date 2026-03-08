#!/usr/bin/env python3
"""
FAST Ω dim computation using QR factorization instead of SVD.
QR is O(mn²) vs SVD O(mn min(m,n)), much faster for rectangular matrices.

Also: since all eigenspaces have the same Ω dims for Paley,
we only need k=1.
"""
import numpy as np
from scipy import sparse
from scipy.linalg import qr as scipy_qr
import sys, time
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

import os
old_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_fourier_v3 import enumerate_step_sequences, compute_faces
sys.stdout = old_stdout

def qr_set(p):
    return set((a*a) % p for a in range(1, p))

def fast_nullity(M, tol=1e-8):
    """Compute nullity of M using QR with pivoting."""
    m, n = M.shape
    if m == 0 or n == 0:
        return n

    # QR with column pivoting
    Q, R, P = scipy_qr(M, pivoting=True)
    # Rank = number of diagonal elements of R above tolerance
    diag = np.abs(np.diag(R[:min(m,n), :min(m,n)]))
    rank = np.sum(diag > tol)
    return n - rank

def compute_omega_dims_fast(S_set, n, k, max_dim):
    """Compute Ω dims for eigenspace k using QR factorization."""
    lam = np.exp(2j * np.pi * k / n)

    step_seqs = {}
    for pp in range(max_dim + 2):
        t0 = time.time()
        step_seqs[pp] = enumerate_step_sequences(S_set, n, pp)
        t1 = time.time()
        print(f"  dim {pp}: {len(step_seqs[pp])} seqs ({t1-t0:.1f}s)", flush=True)

    seq_sets = {pp: set(step_seqs[pp]) for pp in range(max_dim + 2)}

    omega_dims = {}
    for pp in range(max_dim + 2):
        dim_Ap = len(step_seqs[pp])
        if pp == 0:
            omega_dims[0] = 1 if dim_Ap > 0 else 0
            continue
        if dim_Ap == 0:
            omega_dims[pp] = 0
            continue

        # Compute junk matrix
        junk_set = set()
        face_data = []
        t0 = time.time()
        for seq in step_seqs[pp]:
            faces = compute_faces(seq, n)
            classified = []
            for face_steps, sign, shift in faces:
                is_valid = (face_steps == ()) if pp - 1 == 0 else (face_steps in seq_sets[pp - 1])
                if not is_valid:
                    junk_set.add(face_steps)
                classified.append((face_steps, sign, shift, is_valid))
            face_data.append(classified)

        junk_list = sorted(junk_set)
        n_junk = len(junk_list)
        t1 = time.time()

        if n_junk == 0:
            omega_dims[pp] = dim_Ap
            print(f"  Ω_{pp} = {dim_Ap} (no junk, {t1-t0:.1f}s)", flush=True)
            continue

        # Build dense junk matrix
        junk_idx = {j: i for i, j in enumerate(junk_list)}
        J = np.zeros((n_junk, dim_Ap), dtype=complex)
        for j_col, classified_faces in enumerate(face_data):
            for face_steps, sign, shift, is_valid in classified_faces:
                if not is_valid and face_steps in junk_idx:
                    row = junk_idx[face_steps]
                    J[row, j_col] += sign * (lam ** shift)

        t0 = time.time()
        null_dim = fast_nullity(J)
        omega_dims[pp] = null_dim
        t1 = time.time()
        print(f"  Ω_{pp} = {null_dim} (J: {n_junk}×{dim_Ap}, {t1-t0:.1f}s)", flush=True)

    return omega_dims

# ===== P_11 =====
print("=" * 70)
print("P_11: FAST Ω DIMS (QR factorization, k=1)")
print("=" * 70)

p = 11
S = qr_set(p)
print(f"QR = {sorted(S)}\n")

dims = compute_omega_dims_fast(S, p, k=1, max_dim=p-1)
omega_list = [dims.get(i, 0) for i in range(p+1)]
print(f"\nΩ dims (k=1): {omega_list[:p+1]}")

# Alternating sum
alt_sum = sum((-1)**i * omega_list[i] for i in range(p))
print(f"Alt sum = {alt_sum} (should be 1)")

# Check palindromicity of Ω_1...Ω_{p-1}
inner = omega_list[1:p]
is_pal = all(inner[i] == inner[-(i+1)] for i in range(len(inner)//2))
print(f"Inner palindromic? {is_pal}")
print(f"Inner sequence: {inner}")

# ===== Predict Betti numbers =====
print(f"\n{'='*70}")
print("BETTI NUMBER PREDICTION")
print("="*70)

# For P_7, we know the boundary ranks
# [1, 2, 4, 5, 3, 3] for eigenspace k≠0
# Can we predict them from Ω dims?

# The boundary ranks satisfy:
# β_0 = Ω_0 - r_1 = 0 → r_1 = 1
# β_m = Ω_m - r_m - r_{m+1} = 0 for most m, = 1 for exactly one m

# If β_d = 1 and all other β = 0 (for d > 0):
# r_1 = 1
# r_m = Ω_m - r_{m-1} for m < d (from β_m = 0)
# At m=d: β_d = Ω_d - r_d - r_{d+1} = 1
# For m > d: β_m = 0 continues

# This means r_m is determined by the Ω dims:
# r_1 = 1
# r_2 = Ω_1 - r_1 = Ω_1 - 1
# r_3 = Ω_2 - r_2 = Ω_2 - Ω_1 + 1
# r_m = Ω_{m-1} - r_{m-1} = Σ_{i=0}^{m-1} (-1)^i Ω_{m-1-i} - (-1)^m

# More simply: r_m = Σ_{i=0}^{m-2} (-1)^i Ω_{m-1-i}
# (cumulative alternating sum from the right)

# For β_d = 1: Ω_d - r_d - r_{d+1} = 1
# This determines d as the unique dimension where the "excess" appears.

print("\nBoundary rank prediction for per-eigenspace Betti:")
print("If pattern is β_d=1 with all other β=0:")

r = [0]  # r_0 = 0 (no boundary into dim 0 from below)
# Actually r_1 = rank(∂_1) = Ω_0 - β_0 = 1 - 0 = 1 (if β_0=0 for k≠0)
r.append(1)  # r_1 = 1

for m in range(2, p):
    r_m = omega_list[m-1] - r[-1]
    r.append(r_m)

print(f"Predicted r = {r[1:]}")

# Check where β would be nonzero
for m in range(p):
    if m < len(r) - 1:
        ker_m = omega_list[m] - r[m]
        im_m1 = r[m+1]
        beta_m = ker_m - im_m1
    else:
        ker_m = omega_list[m] - r[m]
        beta_m = ker_m
    if beta_m != 0:
        print(f"  β_{m} = {beta_m} (ker={ker_m}, im={r[m+1] if m+1<len(r) else 0})")

# Verify for P_7
print(f"\nVerification for P_7 (Ω = [1,3,6,9,9,6,3]):")
o7 = [1, 3, 6, 9, 9, 6, 3]
r7 = [0, 1]
for m in range(2, 7):
    r7.append(o7[m-1] - r7[-1])
print(f"  r = {r7[1:]}")
for m in range(7):
    ker = o7[m] - r7[m]
    im_next = r7[m+1] if m+1 < len(r7) else 0
    beta = ker - im_next
    if beta != 0:
        print(f"  β_{m} = {beta}")

print("\nDone.")
