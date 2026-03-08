#!/usr/bin/env python3
"""
P_11 Ω_7 computation — low memory approach.

The dim 7 junk matrix is 15230×8735. Full complex matrix = 2.1 GB.
QR decomposition needs ~3× that. Too much for 8GB RAM.

APPROACH: Use randomized rank estimation.
Multiply J by random vectors and check if rank < full.
Or use iterative SVD to find smallest singular values.

Actually, better approach: compute rank via LU decomposition on rows,
processing in batches to avoid full matrix allocation.

BEST APPROACH: Since we only need the RANK (not the full decomposition),
use the randomized SVD to estimate it.
"""
import numpy as np
from scipy import sparse
import sys, time, gc
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

import os
old_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_fourier_v3 import enumerate_step_sequences, compute_faces
sys.stdout = old_stdout

def qr_set(p):
    return set((a*a) % p for a in range(1, p))

p = 11
S = qr_set(p)
lam = np.exp(2j * np.pi * 1 / p)  # k=1 eigenspace

# Known: Ω_0..6 = [1, 5, 20, 70, 205, 460, 700]
# Need: Ω_7, 8, 9, 10
# Constraint: alt sum = 1, so Ω_7 - Ω_8 + Ω_9 - Ω_10 = 1-5+20-70+205-460+700 = 391

print("Enumerating step sequences...", flush=True)
step_seqs = {}
for pp in range(12):
    step_seqs[pp] = enumerate_step_sequences(S, p, pp)
    print(f"  dim {pp}: {len(step_seqs[pp])}", flush=True)

seq_sets = {pp: set(step_seqs[pp]) for pp in range(12)}

# For dims 7-10, use sparse matrix + iterative approach
for pp in [7, 8, 9, 10]:
    dim_Ap = len(step_seqs[pp])
    if dim_Ap == 0:
        print(f"Ω_{pp} = 0", flush=True)
        continue

    junk_set = set()
    face_data = []
    t0 = time.time()
    for seq in step_seqs[pp]:
        faces = compute_faces(seq, p)
        classified = []
        for face_steps, sign, shift in faces:
            is_valid = (face_steps == ()) if pp == 1 else (face_steps in seq_sets[pp - 1])
            if not is_valid:
                junk_set.add(face_steps)
            classified.append((face_steps, sign, shift, is_valid))
        face_data.append(classified)

    junk_list = sorted(junk_set)
    n_junk = len(junk_list)
    t1 = time.time()
    print(f"dim {pp}: {n_junk} junk types, {dim_Ap} seqs ({t1-t0:.1f}s)", flush=True)

    if n_junk == 0:
        print(f"Ω_{pp} = {dim_Ap}", flush=True)
        del face_data
        gc.collect()
        continue

    # Build sparse junk matrix
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

    J_sparse = sparse.csr_matrix((vals, (rows, cols)), shape=(n_junk, dim_Ap), dtype=complex)
    t1 = time.time()
    print(f"  Sparse J: {J_sparse.nnz} nonzeros ({t1-t0:.1f}s)", flush=True)

    del face_data, junk_list, junk_idx, rows, cols, vals
    gc.collect()

    # APPROACH 1: Randomized rank estimation
    # Multiply J by random tall-skinny matrix and check rank
    # J has at most min(n_junk, dim_Ap) rank.
    # If we pick k random vectors and compute J @ R, rank(J@R) = rank(J) w.h.p.
    # when k >= rank(J) + oversampling.

    # First try: probe with increasing k
    t0 = time.time()
    max_rank = min(n_junk, dim_Ap)

    # Use random projection: J @ R where R is dim_Ap x k, k = max_rank + 10
    # Then rank(J@R) = rank(J)
    k_proj = min(max_rank + 10, dim_Ap)
    print(f"  Randomized rank: projecting to k={k_proj}...", flush=True)

    np.random.seed(42)
    R = np.random.randn(dim_Ap, k_proj) + 1j * np.random.randn(dim_Ap, k_proj)
    R /= np.sqrt(2)

    # Compute J @ R (sparse × dense)
    JR = J_sparse @ R  # n_junk × k_proj, complex
    t1 = time.time()
    print(f"  J@R computed: {JR.shape}, {JR.nbytes/(1024**2):.0f}MB ({t1-t0:.1f}s)", flush=True)

    del R
    gc.collect()

    # Now SVD of JR to get rank
    t0 = time.time()
    S_vals = np.linalg.svd(JR, compute_uv=False)
    rank = int(np.sum(S_vals > 1e-8))
    null_dim = dim_Ap - rank
    t1 = time.time()
    print(f"Ω_{pp} = {null_dim} (rank {rank}/{max_rank}, {t1-t0:.1f}s)", flush=True)

    del JR, J_sparse, S_vals
    gc.collect()

# Summary
print("\nKnown Ω dims for P_11 (k=1):", flush=True)
print("  [1, 5, 20, 70, 205, 460, 700, ?, ?, ?, ?]", flush=True)
print("  Alt sum constraint: Ω_7 - Ω_8 + Ω_9 - Ω_10 = 391", flush=True)

print("\nDone.", flush=True)
