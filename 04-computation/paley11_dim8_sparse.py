#!/usr/bin/env python3
"""Compute Ω_8 for P_11 (k≠0 eigenspace) using sparse methods.

The previous approach (dense J^H J eigvalsh on 14395×14395) was killed by OOM.
This script uses sparse SVD via scipy.sparse.linalg.svds to find the nullity
of J directly, without forming the dense J^H J matrix.

Key insight: Ω_8 = nullity(J) where J is 35145×14395 with only 50375 nonzeros.
We expect Ω_8 to be in the range 200-700 based on the pattern.
"""
import numpy as np
from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse.linalg import svds, LinearOperator
import sys, time, gc

sys.stdout.reconfigure(line_buffering=True)

# From the previous computation, we know:
# dim 8: 35145 junk types, 14395 seqs, J: 35145x14395, 50375 nonzeros

# Paley tournament P_11: QR = {1,2,3,4,5,9} (mod 11)
p = 11
QR = {pow(x, 2, p) for x in range(1, p)}  # quadratic residues
print(f"P_{p}: QR = {sorted(QR)}")

# Step sequences: length-m+1 sequences of steps from QR
# For dim 8: sequences of 9 steps (8+1 vertices), each step in QR
def enumerate_step_sequences(m, S):
    """Enumerate step sequences of length m from step set S ⊂ Z_p."""
    if m == 0:
        return [()]
    shorter = enumerate_step_sequences(m - 1, S)
    result = []
    for seq in shorter:
        partial_sum = sum(seq) % p
        for s in sorted(S):
            # Check that the new vertex (partial_sum + s) mod p is distinct
            # from all previous vertices
            new_vert = (partial_sum + s) % p
            # Previous vertices: 0, seq[0], seq[0]+seq[1], ...
            verts = [0]
            running = 0
            valid = True
            for step in seq:
                running = (running + step) % p
                verts.append(running)
            if new_vert in verts:
                valid = False
            if valid:
                result.append(seq + (s,))
    return result

print("Enumerating step sequences...")
t0 = time.time()
seqs = enumerate_step_sequences(8, sorted(QR))  # 8 steps for dim 8 (9 vertices)
print(f"  dim 8: {len(seqs)} sequences ({time.time()-t0:.1f}s)")
assert len(seqs) == 14395, f"Expected 14395, got {len(seqs)}"

# Junk matrix: for each way to merge two consecutive steps into one,
# check if the merged step is NOT in QR (= "illegal merge")
# J rows = junk types (illegal merges), J cols = sequences

# For a sequence (s_0, s_1, ..., s_8), merging at position i gives
# step s_i + s_{i+1} (mod p). If this is NOT in QR, we get a junk constraint.
# The junk equation is: sum over all seqs with the same context and merged step = 0.

# Actually, the junk matrix maps sequence to boundary faces:
# ∂(v_0, v_1, ..., v_9) = sum_{i=0}^9 (-1)^i (v_0,...,v̂_i,...,v_9)
# The face with v_i removed corresponds to merging steps s_{i-1} and s_i
# (for interior i) or truncating (for i=0 or i=9).

# For eigenspace decomposition, the junk matrix acts on the step-sequence
# representation. A "junk type" is an 8-step sequence that arises from
# deleting one step position and merging.

# This is complex. Let me just rebuild J from the computation.

# Actually, let me reconstruct J from scratch.
# For dim m, the junk matrix has:
# - Columns indexed by m-step sequences (m+1 vertices, m steps)
# - Rows indexed by "non-allowed" (m-1)-step sequences that appear as faces
# For each m-step sequence s = (s_0,...,s_{m-1}), its boundary faces are:
#   face_i: delete vertex i (for i=0,...,m)
#   face_i has step sequence obtained by:
#     i=0: (s_1, ..., s_{m-1}) — drop first step
#     i=m: (s_0, ..., s_{m-2}) — drop last step
#     0<i<m: (s_0,...,s_{i-2}, s_{i-1}+s_i, s_{i+1},...,s_{m-1}) — merge steps i-1,i
# A face is "non-allowed" if the resulting step sequence contains a merged step
# NOT in QR, OR if the resulting path revisits a vertex.

# First, enumerate allowed 8-step sequences (dim 7 paths)
seqs7 = enumerate_step_sequences(7, sorted(QR))  # 7 steps for dim 7 (8 vertices)
print(f"  dim 7: {len(seqs7)} sequences")
seqs7_set = set(seqs7)

# Build junk matrix: for each 9-step sequence, find non-allowed faces
print("Building junk matrix...")
t0 = time.time()

seq_to_col = {s: i for i, s in enumerate(seqs)}
junk_types = {}  # map non-allowed face tuple -> row index
rows, cols, vals = [], [], []

m = 8  # 8 steps for dim 8 (8+1=9 vertices in path)
for j, seq in enumerate(seqs):
    # Generate all faces
    for i in range(m + 1):
        sign = (-1) ** i
        if i == 0:
            face = seq[1:]
        elif i == m:
            face = seq[:-1]
        else:
            # Merge steps i-1 and i
            merged = (seq[i-1] + seq[i]) % p
            face = seq[:i-1] + (merged,) + seq[i+1:]

        # Check if face is an allowed 8-step sequence
        if face not in seqs7_set:
            # Check vertex distinctness of the face path
            verts = [0]
            running = 0
            distinct = True
            for s in face:
                running = (running + s) % p
                if running in verts:
                    distinct = False
                    break
                verts.append(running)
            if not distinct:
                continue  # degenerate face, skip

            if face not in junk_types:
                junk_types[face] = len(junk_types)
            row = junk_types[face]
            rows.append(row)
            cols.append(j)
            vals.append(sign)

n_junk = len(junk_types)
n_seq = len(seqs)
print(f"  Junk matrix: {n_junk}x{n_seq}, {len(rows)} nonzeros ({time.time()-t0:.1f}s)")

J = csr_matrix((vals, (rows, cols)), shape=(n_junk, n_seq), dtype=np.float64)

# Ω_8 = nullity(J) = n_seq - rank(J)
# Use sparse SVD to estimate rank

# Strategy 1: compute a few of the smallest singular values
# If nullity ≈ 500, we need k ≈ 600 singular values from the bottom
# This is expensive. Better: compute from the top.

# Strategy 2: use randomized rank estimation
# Multiply J by random vectors and check rank
print("\nEstimating rank via randomized projection...")
t0 = time.time()

# Project onto random subspace
np.random.seed(42)
k_proj = 1000  # project to k_proj dimensions
R = np.random.randn(n_seq, k_proj)
JR = J @ R  # n_junk × k_proj
rank_JR = np.linalg.matrix_rank(JR, tol=1e-8)
print(f"  rank(J @ R_{{{n_seq}x{k_proj}}}) = {rank_JR} ({time.time()-t0:.1f}s)")

# If rank_JR < k_proj, then rank(J) = rank_JR
# If rank_JR = k_proj, need more columns

if rank_JR < k_proj:
    print(f"  rank(J) = {rank_JR}")
    print(f"  Ω_8 = {n_seq} - {rank_JR} = {n_seq - rank_JR}")
else:
    print(f"  rank(J) ≥ {k_proj}, trying larger projection...")
    k_proj2 = 2000
    R2 = np.random.randn(n_seq, k_proj2)
    JR2 = J @ R2
    rank_JR2 = np.linalg.matrix_rank(JR2, tol=1e-8)
    print(f"  rank(J @ R_{{{n_seq}x{k_proj2}}}) = {rank_JR2} ({time.time()-t0:.1f}s)")

    if rank_JR2 < k_proj2:
        print(f"  rank(J) = {rank_JR2}")
        print(f"  Ω_8 = {n_seq} - {rank_JR2} = {n_seq - rank_JR2}")
    else:
        # Need to compute more carefully
        # Use iterative approach: keep adding columns until rank saturates
        print("  Rank saturated. Using full J^T J approach with memory mapping...")

        # Alternative: use sparse LU to find nullity
        from scipy.sparse.linalg import splu
        try:
            # Convert to CSC for splu
            J_csc = J.tocsc()
            lu = splu(J_csc[:n_seq, :])  # square part
            print(f"  (sparse LU failed, matrix not square)")
        except:
            pass

        # Try: J^T J as sparse matrix (might not be sparse though)
        JTJ = (J.T @ J).toarray()
        gc.collect()
        print(f"  J^T J formed: {JTJ.shape}, {JTJ.nbytes/(1024**3):.1f}GB")

        # Use eigvalsh with subset_by_index to only compute smallest eigenvalues
        from scipy.linalg import eigvalsh
        t1 = time.time()
        # Compute only the 1000 smallest eigenvalues
        try:
            small_eigs = eigvalsh(JTJ, subset_by_index=[0, 999])
            n_zero = sum(1 for e in small_eigs if abs(e) < 1e-8)
            print(f"  {n_zero} eigenvalues < 1e-8 among smallest 1000 ({time.time()-t1:.1f}s)")
            if n_zero < 1000:
                print(f"  rank(J) = {n_seq - n_zero}")
                print(f"  Ω_8 = {n_zero}")
            else:
                print("  Need to check more eigenvalues!")
        except Exception as e:
            print(f"  eigvalsh failed: {e}")
            # Fallback: matrix_rank
            rank_J = np.linalg.matrix_rank(JTJ, tol=1e-8)
            print(f"  rank(J^T J) = {rank_J}")
            print(f"  Ω_8 = {n_seq - rank_J}")
            del JTJ

print(f"\nTotal time: {time.time()-t0:.1f}s")
print("Done.")
