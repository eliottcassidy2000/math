#!/usr/bin/env python3
"""
P_11 Ω dims — clever rank computation.

The junk matrix J at dim 7 is 15230×8735 with only 26380 nonzeros.
That's about 3 nonzeros per row on average!

KEY INSIGHT: Most rows have very few nonzeros. We can:
1. Use Gaussian elimination on the sparse matrix
2. Or use the SPQR (sparse QR) from SuiteSparse

Actually, simplest approach: use scipy.sparse.linalg.LinearOperator
with iterative eigenvalue computation to find the rank.

ALTERNATIVE: Since we only need nullity(J), and J is m×n with m>n,
we can compute J^H J (n×n Hermitian positive semidefinite) and find
its rank. J^H J is sparse if J is sparse.

EVEN BETTER: For such a sparse matrix, the null space is likely
high-dimensional. Use randomized null space finder:
Pick random vectors x, compute Jx. If Jx ≈ 0, x is close to null space.
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import svds, eigsh
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
lam = np.exp(2j * np.pi * 1 / p)

print("Enumerating step sequences...", flush=True)
step_seqs = {}
for pp in range(12):
    step_seqs[pp] = enumerate_step_sequences(S, p, pp)
    if len(step_seqs[pp]) > 0:
        print(f"  dim {pp}: {len(step_seqs[pp])}", flush=True)

seq_sets = {pp: set(step_seqs[pp]) for pp in range(12)}

def compute_nullity_sparse(pp):
    """Compute nullity of junk matrix at dimension pp using sparse methods."""
    dim_Ap = len(step_seqs[pp])
    if dim_Ap == 0:
        return 0

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
        return dim_Ap

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

    J = sparse.csr_matrix((vals, (rows, cols)), shape=(n_junk, dim_Ap), dtype=complex)
    nnz = J.nnz
    t1 = time.time()
    print(f"  Sparse J: {n_junk}x{dim_Ap}, {nnz} nonzeros ({t1-t0:.1f}s)", flush=True)

    del face_data, junk_list, junk_idx, rows, cols, vals
    gc.collect()

    # APPROACH: Compute J^H J (Hermitian, n×n), then find its rank
    # J^H J is positive semidefinite, rank(J^H J) = rank(J)
    t0 = time.time()
    JHJ = (J.conj().T @ J).toarray()  # dim_Ap × dim_Ap, dense
    t1 = time.time()
    mem_mb = JHJ.nbytes / (1024**2)
    print(f"  J^H J: {dim_Ap}x{dim_Ap}, {mem_mb:.0f}MB ({t1-t0:.1f}s)", flush=True)

    del J
    gc.collect()

    # SVD of J^H J to find rank
    t0 = time.time()
    eigenvalues = np.linalg.eigvalsh(JHJ)  # eigenvalues of Hermitian matrix
    rank = int(np.sum(eigenvalues > 1e-12))  # eigenvalues, not singular values
    null_dim = dim_Ap - rank
    t1 = time.time()
    print(f"Ω_{pp} = {null_dim} (rank {rank}/{min(n_junk,dim_Ap)}, {t1-t0:.1f}s)", flush=True)

    # Show some eigenvalue info
    sorted_evals = np.sort(eigenvalues)
    print(f"  Smallest eigenvalues: {sorted_evals[:5]}", flush=True)
    print(f"  Largest eigenvalues: {sorted_evals[-3:]}", flush=True)

    del JHJ, eigenvalues
    gc.collect()

    return null_dim

# Process dims 7-10
print("\n--- Computing Ω dims 7-10 ---", flush=True)
results = {}
for pp in [7, 8, 9, 10]:
    t0 = time.time()
    results[pp] = compute_nullity_sparse(pp)
    t1 = time.time()
    print(f"  Total time for dim {pp}: {t1-t0:.1f}s\n", flush=True)

# Summary
known = {0:1, 1:5, 2:20, 3:70, 4:205, 5:460, 6:700}
known.update(results)

print("\n" + "="*70, flush=True)
print("P_11 Ω dims (k=1 eigenspace):", flush=True)
omega_list = [known.get(d, '?') for d in range(11)]
print(f"  {omega_list}", flush=True)

# Alternating sum
alt_sum = sum((-1)**d * v for d, v in enumerate(omega_list) if isinstance(v, int))
print(f"  Alt sum = {alt_sum} (should be 1)", flush=True)

# Palindromic check (inner sequence)
inner = omega_list[1:-1] if omega_list[-1] == 0 else omega_list[1:]
print(f"  Inner: {inner}", flush=True)
if all(isinstance(v, int) for v in inner):
    is_pal = all(inner[i] == inner[-(i+1)] for i in range(len(inner)//2))
    print(f"  Palindromic? {is_pal}", flush=True)

print("\nDone.", flush=True)
