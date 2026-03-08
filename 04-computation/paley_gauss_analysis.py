#!/usr/bin/env python3
"""
PALEY-GAUSS CONNECTION IN PATH HOMOLOGY

The Paley tournament P_p has path homology controlled by:
1. The Gauss sum g = Σ_{a=1}^{p-1} χ(a) ω^a where χ is the Legendre symbol
2. The character sum Σ_{s ∈ QR} ω^{ks} = (g * χ(k) - 1) / 2

Key observation: For P_7, β = (1,0,0,0,6,0) with χ = 7.
This is 1 + 6 = 7 = p. Is χ = p for all Paley tournaments?

For general circulant tournaments C_p^S (p prime, |S| = (p-1)/2, S ∩ (-S) = ∅):
- Ω_m^(λ) depends on eigenspace λ
- β_p = Σ_λ β_p^(λ)
- The Euler characteristic χ = Σ (-1)^p β_p

ANALYSIS: Compute χ for various circulant tournaments and relate to p.
"""
import numpy as np
from itertools import combinations
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

# Suppress import output
import os
old_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_fourier_v3 import fourier_betti_v3
from path_homology_v2 import path_betti_numbers, circulant_digraph
sys.stdout = old_stdout

def qr(p):
    return set((a*a) % p for a in range(1, p))

def is_prime(n):
    if n < 2: return False
    for d in range(2, int(n**0.5)+1):
        if n % d == 0: return False
    return True

# ===== 1. Euler characteristic of circulant tournaments =====
print("=" * 70)
print("EULER CHARACTERISTIC OF CIRCULANT TOURNAMENTS")
print("=" * 70)

print("\nA circulant tournament on Z_p has S with |S|=(p-1)/2, S∩(-S)=∅.")
print("Question: Is χ(P_p) = p always?\n")

for p in [3, 5, 7]:
    if not is_prime(p): continue
    print(f"\np={p}:", flush=True)

    # All circulant tournaments at this prime
    all_ct = []
    for S_comb in combinations(range(1, p), (p-1)//2):
        S_set = set(S_comb)
        neg_S = {(p - s) % p for s in S_set}
        if S_set & neg_S == set():
            A = circulant_digraph(p, sorted(S_set))
            betti = path_betti_numbers(A, p, max_dim=p-1)
            chi = sum((-1)**k * betti[k] for k in range(len(betti)))
            is_qr = (S_set == qr(p))
            label = " (PALEY)" if is_qr else ""
            all_ct.append((S_set, betti, chi))
            print(f"  S={sorted(S_set)}: β={betti}, χ={chi}{label}")

    # Check if all χ = p
    all_chi_p = all(chi == p for _, _, chi in all_ct)
    print(f"  All χ = p = {p}? {all_chi_p}")

# ===== 2. Gauss sum decomposition =====
print("\n\n" + "=" * 70)
print("GAUSS SUM DECOMPOSITION")
print("=" * 70)

for p in [3, 7, 11, 19, 23]:
    if not is_prime(p) or p % 4 != 3: continue
    S = qr(p)
    omega = np.exp(2j * np.pi / p)

    print(f"\nP_{p}: QR={sorted(S)}")

    # Gauss sum
    g = sum(omega**(a*a) for a in range(1, p))
    print(f"  Gauss sum g = {g.real:.4f} + {g.imag:.4f}i")
    print(f"  |g| = {abs(g):.4f}, sqrt(p) = {p**0.5:.4f}")
    print(f"  g^2 = {(g*g).real:.4f} + {(g*g).imag:.4f}i")
    print(f"  Expected g^2 = {(-1)**((p-1)//2) * p} (= (-1)^{(p-1)//2} * p)")

    # Character sums at each eigenvalue
    print(f"  Eigenspace analysis:")
    for k in range(min(p, 10)):
        lam = omega ** k
        # Σ_{s∈QR} lam^s
        char_sum = sum(lam**s for s in S)
        legendre_k = 1 if k in S else (-1 if k != 0 else 0)
        if k > 0:
            expected = (legendre_k * g - 1) / 2
            print(f"    k={k}: Σ_QR λ^s = {char_sum.real:+.3f}{char_sum.imag:+.3f}i, "
                  f"χ(k)={legendre_k}, (χ(k)g-1)/2 = {expected.real:+.3f}{expected.imag:+.3f}i")
        else:
            print(f"    k=0: Σ_QR λ^s = {char_sum.real:.3f} = (p-1)/2 = {(p-1)/2}")

# ===== 3. Why β_4=6 for P_7? =====
print("\n\n" + "=" * 70)
print("WHY β_4 = 6 FOR P_7?")
print("=" * 70)

p = 7
S = qr(p)
omega = np.exp(2j * np.pi / p)

print(f"\nP_7 = C_7^{{1,2,4}}, S = QR(7) = {{1,2,4}}")
print(f"β = [1,0,0,0,6,0,0], χ = 7 = p")
print(f"\nThe 6 in β_4 must come from 6 eigenspaces each contributing β_4^(λ)=1.")
print(f"Since p=7, there are 7 eigenspaces (k=0,...,6).")
print(f"One eigenspace (k=0, λ=1) contributes β_0=1.")
print(f"So the other 6 eigenspaces (k=1,...,6) each contribute β_4=1?")

# Verify: compute per-eigenspace Betti numbers
# We'll need to modify fourier_betti_v3 to return per-eigenspace data

from path_homology_fourier_v3 import (
    enumerate_step_sequences, compute_faces
)

max_dim = 6
step_seqs = {}
for pp in range(max_dim + 2):
    step_seqs[pp] = enumerate_step_sequences(S, p, pp)

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
        faces = compute_faces(seq, p)
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

print(f"\nStep sequence counts:")
for pp in range(max_dim + 2):
    print(f"  p={pp}: |A_{pp}|={len(step_seqs[pp])}, |junk types|={len(junk_types.get(pp, []))}")

# Per-eigenvalue β
print(f"\nPer-eigenspace Betti numbers:")
roots = [np.exp(2j * np.pi * k / p) for k in range(p)]

for k in range(p):
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

    # Compute boundary ranks (simplified: using full A-space boundaries)
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

    nonzero = [i for i in range(len(betti_k)) if betti_k[i] > 0]
    if nonzero:
        print(f"  k={k}: Ω dims={[omega_dims.get(i,0) for i in range(max_dim+1)]}, "
              f"β^(k)={betti_k}")

print("\nDone.")
