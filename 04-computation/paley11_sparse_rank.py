#!/usr/bin/env python3
"""
P_11 Ω dims — sparse rank computation.

The junk matrix J at dim 7 is 15230×8735 with only 26380 nonzeros
(density 0.0002). Use sparse SVD to find the rank.

APPROACH: The nullity = dim_Ap - rank(J). We need rank(J).
Since J is very sparse, use scipy.sparse.linalg.svds to find
the largest singular values, then estimate rank from them.

Actually better: use the formula nullity = dim(ker J).
Build J^H J (8735×8735 sparse Hermitian), then find its nullity
via sparse eigenvalue computation.

Best: Use randomized low-rank approximation.
Since J has 26380 nonzeros in 133 million entries, it's very low rank.
Actually no — rank could be up to 8735.

Let me try: sparse QR via scipy.sparse.linalg, or use
the fact that J only has about 3 nonzeros per row on average.
With 15230 rows and rank potentially up to 8735, we need care.

ALTERNATIVE: Process J in column blocks. Build J, convert to dense
column-by-column, use incremental SVD.

SIMPLEST: Just build the sparse J and convert to dense, but use
float32 instead of complex128 to save memory.

Wait — actually, let me try: for k=1, the eigenvalue λ = ω = e^{2πi/11}.
The entries are sign * λ^shift. Can we work modulo some prime instead?
For rank computation, we just need to know the rank over C.

Simplest memory reduction: Build J as complex64 (half the memory).
15230 × 8735 × 8 bytes = 1.06 GB (vs 2.1 GB for complex128).
"""
import numpy as np
from scipy import sparse
from scipy.linalg import qr as scipy_qr
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
lam = np.exp(2j * np.pi * 1 / p).astype(np.complex64)

print("Enumerating step sequences...", flush=True)
step_seqs = {}
for pp in range(12):
    step_seqs[pp] = enumerate_step_sequences(S, p, pp)
    if len(step_seqs[pp]) > 0:
        print(f"  dim {pp}: {len(step_seqs[pp])}", flush=True)

seq_sets = {pp: set(step_seqs[pp]) for pp in range(12)}

# Process dims 7-10
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

    # Build dense junk matrix in complex64 to save memory
    junk_idx = {j: i for i, j in enumerate(junk_list)}
    t0 = time.time()

    J = np.zeros((n_junk, dim_Ap), dtype=np.complex64)
    for j_col, classified_faces in enumerate(face_data):
        for face_steps, sign, shift, is_valid in classified_faces:
            if not is_valid and face_steps in junk_idx:
                row = junk_idx[face_steps]
                J[row, j_col] += sign * (lam ** shift)

    t1 = time.time()
    mem_mb = J.nbytes / (1024**2)
    print(f"  Built J: {n_junk}x{dim_Ap}, {mem_mb:.0f}MB ({t1-t0:.1f}s)", flush=True)

    del face_data, junk_list, junk_idx
    gc.collect()

    # Compute rank via SVD on complex64
    t0 = time.time()
    try:
        S_vals = np.linalg.svd(J, compute_uv=False)
        rank = int(np.sum(S_vals > 1e-4))  # lower threshold for float32
        null_dim = dim_Ap - rank
        t1 = time.time()
        print(f"Ω_{pp} = {null_dim} (rank {rank}/{min(n_junk,dim_Ap)}, {t1-t0:.1f}s)", flush=True)
    except MemoryError:
        print(f"  OOM on SVD for dim {pp}. Trying QR on transpose...", flush=True)
        # If J is tall (n_junk > dim_Ap), work with J^H instead
        if n_junk > dim_Ap:
            JH = J.conj().T  # dim_Ap x n_junk
            del J
            gc.collect()
            t0 = time.time()
            S_vals = np.linalg.svd(JH, compute_uv=False)
            rank = int(np.sum(S_vals > 1e-4))
            null_dim = dim_Ap - rank
            t1 = time.time()
            print(f"Ω_{pp} = {null_dim} (rank {rank}, via JH, {t1-t0:.1f}s)", flush=True)
        else:
            print(f"  SKIP dim {pp}", flush=True)

    del J
    gc.collect()

# Summary
known = {0:1, 1:5, 2:20, 3:70, 4:205, 5:460, 6:700}
print(f"\nΩ dims for P_11 (k=1 eigenspace):", flush=True)
for d in range(11):
    val = known.get(d, "?")
    print(f"  Ω_{d} = {val}", flush=True)

# Alt sum check
computed = [known.get(d) for d in range(11)]
known_sum = sum((-1)**d * v for d, v in enumerate(computed) if v is not None)
print(f"\nPartial alt sum (dims 0-6): {known_sum}", flush=True)
print(f"Need dims 7-10 to sum to {known_sum - 1} (alternating)", flush=True)

print("\nDone.", flush=True)
