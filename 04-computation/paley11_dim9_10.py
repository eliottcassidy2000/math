#!/usr/bin/env python3
"""Compute Ω_9 and Ω_10 for P_11 (k≠0 eigenspace).

Known: Ω_0-8 = [1, 5, 20, 70, 205, 460, 700, 690, 450]
Need: Ω_9 (from A_9=15745 seqs, A_8=14395 ref) and Ω_10 (from A_10=8645, A_9=15745 ref).

Use the subset_by_index eigvalsh approach that worked for dim 8.
"""
import numpy as np
from scipy.sparse import csr_matrix
from scipy.linalg import eigvalsh
import sys, time, gc

sys.stdout.reconfigure(line_buffering=True)

p = 11
QR = {pow(x, 2, p) for x in range(1, p)}
print(f"P_{p}: QR = {sorted(QR)}")

def enumerate_step_sequences(m, S, p):
    if m == 0:
        return [()]
    shorter = enumerate_step_sequences(m - 1, S, p)
    result = []
    for seq in shorter:
        partial_sum = sum(seq) % p
        for s in sorted(S):
            new_vert = (partial_sum + s) % p
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

S = sorted(QR)

for dim in [9, 10]:
    print(f"\n{'='*60}")
    print(f"DIM {dim}")
    print(f"{'='*60}")

    t0 = time.time()
    seqs = enumerate_step_sequences(dim, S, p)
    n_seq = len(seqs)
    print(f"  A_{dim}: {n_seq} sequences ({time.time()-t0:.1f}s)")

    seqs_ref = enumerate_step_sequences(dim - 1, S, p)
    n_ref = len(seqs_ref)
    seqs_ref_set = set(seqs_ref)
    print(f"  A_{dim-1}: {n_ref} sequences")

    # Build junk matrix
    t0 = time.time()
    junk_types = {}
    rows, cols, vals = [], [], []

    for j, seq in enumerate(seqs):
        m = dim
        for i in range(m + 1):
            sign = (-1) ** i
            if i == 0:
                face = seq[1:]
            elif i == m:
                face = seq[:-1]
            else:
                merged = (seq[i-1] + seq[i]) % p
                face = seq[:i-1] + (merged,) + seq[i+1:]

            if face not in seqs_ref_set:
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
                    continue

                if face not in junk_types:
                    junk_types[face] = len(junk_types)
                row = junk_types[face]
                rows.append(row)
                cols.append(j)
                vals.append(sign)

    n_junk = len(junk_types)
    print(f"  Junk: {n_junk} types, {len(rows)} nonzeros ({time.time()-t0:.1f}s)")

    J = csr_matrix((vals, (rows, cols)), shape=(n_junk, n_seq), dtype=np.float64)

    # Form J^T J
    t0 = time.time()
    JTJ = (J.T @ J).toarray()
    mem_gb = JTJ.nbytes / (1024**3)
    print(f"  J^T J: {JTJ.shape}, {mem_gb:.1f}GB ({time.time()-t0:.1f}s)")
    gc.collect()

    # Compute smallest eigenvalues
    t0 = time.time()
    # Estimate nullity: try 1000 smallest eigenvalues first
    k_try = min(1000, n_seq)
    small_eigs = eigvalsh(JTJ, subset_by_index=[0, k_try - 1])
    n_zero = sum(1 for e in small_eigs if abs(e) < 1e-8)
    print(f"  {n_zero} zero eigenvalues among smallest {k_try} ({time.time()-t0:.1f}s)")

    if n_zero < k_try:
        omega_dim = n_zero
        print(f"  Ω_{dim} = {omega_dim}")
    else:
        print(f"  Need more eigenvalues!")
        # Compute more
        k_try2 = min(n_zero + 500, n_seq)
        small_eigs2 = eigvalsh(JTJ, subset_by_index=[0, k_try2 - 1])
        n_zero2 = sum(1 for e in small_eigs2 if abs(e) < 1e-8)
        omega_dim = n_zero2
        print(f"  {n_zero2} zero eigenvalues among smallest {k_try2}")
        print(f"  Ω_{dim} = {omega_dim}")

    # Print smallest nonzero eigenvalue for confidence
    nonzero_eigs = [e for e in small_eigs if abs(e) > 1e-8]
    if nonzero_eigs:
        print(f"  Smallest nonzero eigenvalue: {min(nonzero_eigs):.6e}")

    del JTJ
    gc.collect()

# Summary
print(f"\n{'='*60}")
print(f"SUMMARY: P_11 k≠0 eigenspace Ω dims")
print(f"{'='*60}")
known = [1, 5, 20, 70, 205, 460, 700, 690, 450]
print(f"  Ω_0-8 = {known}")
print(f"  (dims 9, 10 computed above)")
print("Done.")
