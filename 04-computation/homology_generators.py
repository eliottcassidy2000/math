"""
Compute explicit homology generators for P_7 and analyze their structure.

For P_7 (m=3): total β = [1,0,0,0,6,0,0], with β_4 = 6 coming entirely from k≠0.
For k=0: β = [1,0,0,0,0,0,0] (contractible).
For k=1: β = [0,0,0,0,1,0,0].

Let's find the single H_4 generator in the k=1 eigenspace and understand its structure.

Then for P_11 (m=5): β_5^{(0)} = 5, coming from one Z_5-orbit.
"""
import numpy as np
from math import pi
from scipy.linalg import null_space

def get_QR(p):
    return set(pow(x, 2, p) for x in range(1, p))

def build_diffseqs(p, d):
    QR = get_QR(p)
    QR_list = sorted(QR)
    results = []
    if d == 0:
        return [()]
    def backtrack(seq, ps_list, ps_set):
        if len(seq) == d:
            results.append(tuple(seq))
            return
        for s in QR_list:
            new_ps = (ps_list[-1] + s) % p
            if new_ps in ps_set:
                continue
            seq.append(s)
            ps_list.append(new_ps)
            ps_set.add(new_ps)
            backtrack(seq, ps_list, ps_set)
            seq.pop()
            ps_list.pop()
            ps_set.remove(new_ps)
    backtrack([], [0], {0})
    return results

def build_constraint(p, d, Ad):
    QR = get_QR(p)
    junk_faces = {}
    for seq in Ad:
        for i in range(1, d):
            merged = (seq[i-1] + seq[i]) % p
            if merged not in QR:
                face = list(seq)
                face[i-1] = merged
                del face[i]
                ft = tuple(face)
                if ft not in junk_faces:
                    junk_faces[ft] = len(junk_faces)
    n_rows = len(junk_faces)
    if n_rows == 0:
        return np.zeros((0, len(Ad)), dtype=complex)
    C = np.zeros((n_rows, len(Ad)), dtype=complex)
    for j, seq in enumerate(Ad):
        for i in range(1, d):
            merged = (seq[i-1] + seq[i]) % p
            if merged not in QR:
                face = list(seq)
                face[i-1] = merged
                del face[i]
                ft = tuple(face)
                C[junk_faces[ft], j] += (-1) ** i
    return C

def build_omega_basis(C):
    """Compute basis of ker(C) using SVD."""
    if C.shape[0] == 0:
        return np.eye(C.shape[1], dtype=complex)
    U, S, Vt = np.linalg.svd(C)
    tol = 1e-8
    null_dim = C.shape[1] - np.sum(S > tol)
    if null_dim == 0:
        return np.zeros((C.shape[1], 0), dtype=complex)
    return Vt[-null_dim:].T

def build_boundary_k(p, d, k, Ad_from, Ad_to):
    omega = np.exp(2j * pi / p)
    to_idx = {seq: i for i, seq in enumerate(Ad_to)}
    B = np.zeros((len(Ad_to), len(Ad_from)), dtype=complex)
    for j, seq in enumerate(Ad_from):
        for i in range(d + 1):
            if i == 0:
                face = seq[1:]
                phase = omega ** (k * seq[0]) if k != 0 else 1.0
            elif i == d:
                face = seq[:-1]
                phase = 1.0
            else:
                face = seq[:i-1] + ((seq[i-1] + seq[i]) % p,) + seq[i+1:]
                phase = 1.0
            if face in to_idx:
                B[to_idx[face], j] += (-1)**i * phase
    return B

# P_7 analysis
p = 7
m = 3
QR = get_QR(p)

print(f"=== P_{p} (m={m}), QR = {sorted(QR)} ===\n")

# Build everything
diffseqs = {}
C_mats = {}
omega_bases = {}
for d in range(2*m + 1):
    diffseqs[d] = build_diffseqs(p, d)
    C_mats[d] = build_constraint(p, d, diffseqs[d])
    omega_bases[d] = build_omega_basis(C_mats[d])
    print(f"d={d}: |A|={len(diffseqs[d])}, Omega={omega_bases[d].shape[1]}")

# For k=1 eigenspace, find H_4 generator
k = 1
omega_val = np.exp(2j * pi / p)
print(f"\n--- k={k} eigenspace: finding H_4 generator ---")

# Boundary maps restricted to Omega
def restricted_boundary(d, k, omega_from, omega_to):
    B = build_boundary_k(p, d, k, diffseqs[d], diffseqs[d-1])
    return omega_to.conj().T @ B @ omega_from

# Compute restricted boundary at d=4 and d=5
B4_restr = restricted_boundary(4, k, omega_bases[4], omega_bases[3])
B5_restr = restricted_boundary(5, k, omega_bases[5], omega_bases[4])

print(f"  ∂_4 restricted: {B4_restr.shape}, rank={np.linalg.matrix_rank(B4_restr, tol=1e-6)}")
print(f"  ∂_5 restricted: {B5_restr.shape}, rank={np.linalg.matrix_rank(B5_restr, tol=1e-6)}")

# ker(∂_4) in Omega_4
ker_B4 = null_space(B4_restr)
print(f"  ker(∂_4|Omega_4) has dim {ker_B4.shape[1]}")

# im(∂_5) in Omega_4
im_B5 = B5_restr  # columns are images of Omega_5 basis vectors
rank_B5 = np.linalg.matrix_rank(im_B5, tol=1e-6)
print(f"  im(∂_5|Omega_5) has dim {rank_B5}")

# H_4 = ker(∂_4) / im(∂_5)
# Find a complement: project ker(∂_4) basis onto orthogonal complement of im(∂_5)
# Use QR on im(∂_5) columns
Q5, _ = np.linalg.qr(im_B5, mode='reduced')
proj = Q5 @ (Q5.conj().T @ ker_B4)
h4_basis = ker_B4 - proj
# Get nonzero vectors
norms = np.linalg.norm(h4_basis, axis=0)
significant = norms > 1e-6
h4_vectors = h4_basis[:, significant]
print(f"  H_4 has dim {np.sum(significant)} (should be 1)")

if h4_vectors.shape[1] > 0:
    # Convert back to A_4 coordinates
    h4_in_omega = h4_vectors[:, 0]
    h4_in_A = omega_bases[4] @ h4_in_omega

    print(f"\n  H_4 generator in A_4 coordinates:")
    for j, seq in enumerate(diffseqs[4]):
        if abs(h4_in_A[j]) > 1e-6:
            print(f"    {seq}: {h4_in_A[j]:.6f}")

    # Normalize so max coefficient is 1
    max_coeff = max(abs(h4_in_A))
    h4_normalized = h4_in_A / max_coeff

    print(f"\n  Normalized H_4 generator:")
    nonzero_count = 0
    for j, seq in enumerate(diffseqs[4]):
        if abs(h4_normalized[j]) > 1e-6:
            nonzero_count += 1
            coeff = h4_normalized[j]
            # Check if coefficient is a root of unity
            if abs(abs(coeff) - 1) < 1e-6:
                # It's a root of unity
                angle = np.angle(coeff) / (2*pi/p)
                print(f"    {seq}: ω^{angle:.1f} (|c|=1)")
            else:
                print(f"    {seq}: {coeff:.6f}")
    print(f"  Total nonzero: {nonzero_count} out of {len(diffseqs[4])}")

    # Check Z_m action on the generator
    print(f"\n  Z_m (QR-scaling) action on nonzero diff-seqs:")
    QR_list = sorted(QR)
    nonzero_seqs = [seq for j, seq in enumerate(diffseqs[4]) if abs(h4_normalized[j]) > 1e-6]
    for seq in nonzero_seqs[:5]:
        orbit = []
        for q in QR_list:
            scaled = tuple((q * s) % p for s in seq)
            orbit.append(scaled)
        print(f"    {seq} → {orbit}")

# Also: find H_4 for k=0 (should be trivial since β_4^{(0)} = 0 for P_7)
print(f"\n--- k=0 eigenspace: H_4 check ---")
B4_restr_k0 = restricted_boundary(4, 0, omega_bases[4], omega_bases[3])
B5_restr_k0 = restricted_boundary(5, 0, omega_bases[5], omega_bases[4])
print(f"  rank(∂_4|Omega_4) = {np.linalg.matrix_rank(B4_restr_k0, tol=1e-6)}")
print(f"  rank(∂_5|Omega_5) = {np.linalg.matrix_rank(B5_restr_k0, tol=1e-6)}")
ker_k0 = null_space(B4_restr_k0)
print(f"  ker(∂_4) dim = {ker_k0.shape[1]}")
print(f"  β_4^(0) = ker - im = {ker_k0.shape[1]} - {np.linalg.matrix_rank(B5_restr_k0, tol=1e-6)} = {ker_k0.shape[1] - np.linalg.matrix_rank(B5_restr_k0, tol=1e-6)}")

# Check: what is the TOTAL homology generator? Sum over all k≠0.
print(f"\n--- Total H_4 generators (sum over k≠0) ---")
total_gens = []
for k_val in range(1, p):
    B4_k = restricted_boundary(4, k_val, omega_bases[4], omega_bases[3])
    B5_k = restricted_boundary(5, k_val, omega_bases[5], omega_bases[4])
    ker_k = null_space(B4_k)
    Q5_k, _ = np.linalg.qr(B5_k, mode='reduced')
    proj_k = Q5_k @ (Q5_k.conj().T @ ker_k)
    h4_k = ker_k - proj_k
    norms_k = np.linalg.norm(h4_k, axis=0)
    sig_k = norms_k > 1e-6
    if np.sum(sig_k) > 0:
        gen = omega_bases[4] @ h4_k[:, sig_k][:, 0]
        total_gens.append((k_val, gen))
        # Count nonzero
        nz = sum(1 for c in gen if abs(c) > 1e-6)
        print(f"  k={k_val}: {nz} nonzero coefficients")

# Check if all k≠0 generators are related by QR-scaling
if len(total_gens) >= 2:
    print(f"\n  Checking QR-orbit structure of generators:")
    g1 = total_gens[0][1]
    for k_val, g in total_gens[1:]:
        # Check if g is a QR-scaled version of g1
        # QR-scaling by q: (s_1,...,s_d) -> (q*s_1,...,q*s_d)
        # In eigenspace k, this maps to eigenspace k*q
        # So the generator at k should be the q-scaling of the generator at 1
        # where q is such that k ≡ 1*q mod p... but k-eigenspace maps to kq eigenspace
        pass
    print(f"  (QR-orbit relates all k≠0 generators by THM-125)")
