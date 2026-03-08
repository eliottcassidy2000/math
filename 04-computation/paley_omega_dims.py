#!/usr/bin/env python3
"""
FAST COMPUTATION OF Ω DIMS FOR PALEY TOURNAMENTS

Since all eigenspaces of P_p have the same Ω dims, we only need to
compute one eigenspace (say k=1). The per-eigenspace Betti numbers
then give the total.

For P_11, the junk matrices at high dim are huge but we can compute
the Ω dims one dimension at a time (just need SVD ranks).

OPTIMIZATION: Use sparse matrices and only compute nullity (rank),
not the full SVD.
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import svds
import sys, time
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

import os
old_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_fourier_v3 import enumerate_step_sequences, compute_faces
sys.stdout = old_stdout

def qr(p):
    return set((a*a) % p for a in range(1, p))

def compute_omega_dims_single_eigenspace(S_set, n, k, max_dim):
    """Compute Ω dims for eigenspace k of C_n^S."""
    lam = np.exp(2j * np.pi * k / n)

    step_seqs = {}
    for pp in range(max_dim + 2):
        t0 = time.time()
        step_seqs[pp] = enumerate_step_sequences(S_set, n, pp)
        t1 = time.time()
        if len(step_seqs[pp]) > 0:
            print(f"  enum dim {pp}: {len(step_seqs[pp])} seqs ({t1-t0:.1f}s)", flush=True)

    seq_sets = {}
    for pp in range(max_dim + 2):
        seq_sets[pp] = set(step_seqs[pp])

    omega_dims = {}
    for pp in range(max_dim + 2):
        dim_Ap = len(step_seqs[pp])
        if pp == 0:
            omega_dims[0] = 1 if dim_Ap > 0 else 0
            continue
        if dim_Ap == 0:
            omega_dims[pp] = 0
            continue

        # Compute junk types
        junk_set = set()
        face_data = []
        t0 = time.time()
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
            face_data.append(classified)

        junk_list = sorted(junk_set)
        n_junk = len(junk_list)
        t1 = time.time()
        print(f"  junk dim {pp}: {n_junk} types ({t1-t0:.1f}s)", flush=True)

        if n_junk == 0:
            omega_dims[pp] = dim_Ap
            continue

        # Build junk matrix (sparse)
        junk_idx = {j: i for i, j in enumerate(junk_list)}
        t0 = time.time()

        rows, cols, vals = [], [], []
        for j_col, classified_faces in enumerate(face_data):
            for face_steps, sign, shift, is_valid in classified_faces:
                if not is_valid and face_steps in junk_idx:
                    row = junk_idx[face_steps]
                    rows.append(row)
                    cols.append(j_col)
                    vals.append(sign * (lam ** shift))

        J = sparse.csr_matrix((vals, (rows, cols)), shape=(n_junk, dim_Ap), dtype=complex)
        t1 = time.time()
        print(f"  sparse J dim {pp}: {J.nnz} nonzeros, shape {J.shape} ({t1-t0:.1f}s)", flush=True)

        # Compute rank via dense SVD (sparse SVD unreliable for nullity)
        t0 = time.time()
        if n_junk < 5000 and dim_Ap < 5000:
            J_dense = J.toarray()
            S_vals = np.linalg.svd(J_dense, compute_uv=False)
            rank = np.sum(S_vals > 1e-8)
        else:
            # For very large matrices, try dense anyway (may be slow)
            print(f"    Large matrix {J.shape}, using dense SVD...", flush=True)
            J_dense = J.toarray()
            S_vals = np.linalg.svd(J_dense, compute_uv=False)
            rank = np.sum(S_vals > 1e-8)

        null_dim = dim_Ap - rank
        omega_dims[pp] = null_dim
        t1 = time.time()
        print(f"  Ω_{pp} = {null_dim} (rank {rank}/{min(n_junk,dim_Ap)}, {t1-t0:.1f}s)", flush=True)

    return omega_dims

# ===== P_11 =====
print("=" * 70)
print("Ω DIMS FOR P_11 (single eigenspace k=1)")
print("=" * 70)

p = 11
S = qr(p)
print(f"P_{p}: QR = {sorted(S)}")
dims = compute_omega_dims_single_eigenspace(S, p, k=1, max_dim=p-1)
print(f"\nΩ dims (k=1): {[dims.get(i, 0) for i in range(p+1)]}")

# Alternating sum = χ per eigenspace
omega_list = [dims.get(i, 0) for i in range(p)]
alt_sum = sum((-1)**i * omega_list[i] for i in range(len(omega_list)))
print(f"Alt sum = {alt_sum} (should be 1)")

# Check palindromic
is_pal = all(omega_list[i] == omega_list[-(i+1)] for i in range(len(omega_list)//2))
print(f"Palindromic? {is_pal}")

# Predict Betti numbers
print(f"\nIf all eigenspaces identical:")
print(f"  β_0 = 1 (from k=0)")
for d in range(1, p):
    # Need to determine β_d for single eigenspace
    # β_d = ker(∂_d) - im(∂_{d+1}) on Ω
    # For now just show Ω dims
    pass
print(f"  Total χ = 1 + (p-1)·1 = p = {p}")

print("\nDone.")
