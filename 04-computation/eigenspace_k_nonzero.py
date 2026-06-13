"""
Compute Omega and Betti for k≠0 eigenspaces of P_7 and P_11.

For the k-eigenspace, the boundary map ∂_d has a "twist" by ω^k
where ω = exp(2πi/p). The constraint matrix and boundary map are modified.

The Z_p additive action on paths: v_i → v_i + a. On diff-seqs: unchanged.
So the eigenspace decomposition is on the PATH level, not the diff-seq level.

For a k-eigenspace vector: sum_{a=0}^{p-1} ω^{-ka} |a, a+s_1, a+s_1+s_2, ...|
The boundary of this vector involves faces, which are paths of length d-1.
Face i merges vertex i-1 and i+1 (deleting vertex i from the path).

For the k-eigenspace: the CONSTRAINT matrix C_d is the SAME as k=0
(since it operates on diff-seqs, not paths).
But the BOUNDARY map ∂_d differs!

Wait, actually for k≠0: the boundary map sends a k-eigenspace path to a
k-eigenspace path. On diff-seqs, face i of (s_1,...,s_d) is:
  face 0: (s_2, ..., s_d) — drops first element
  face d: (s_1, ..., s_{d-1}) — drops last element
  face j (0 < j < d): replace (s_j, s_{j+1}) with s_j + s_{j+1}

For face 0: the new starting vertex is s_1 (was 0, now shifted by s_1).
In the k-eigenspace: this introduces a factor of ω^{k·s_1}.

For face d: the starting vertex is still 0, no phase factor.
For face j: the starting vertex is still 0, no phase factor.

So the boundary in the k-eigenspace is:
  ∂_d^{(k)}(s_1,...,s_d) = ω^{k·s_1} · (s_2,...,s_d)
                          + sum_{j=1}^{d-1} (-1)^j · face_j(s_1,...,s_d)
                          + (-1)^d · (s_1,...,s_{d-1})

The only difference from k=0 is the phase factor ω^{k·s_1} on face 0.

The constraint equations (junk faces) don't involve face 0 or face d
(those are always allowed). So the CONSTRAINT MATRIX IS THE SAME for all k.
Therefore Omega_d is the SAME for all k.

The BOUNDARY MAP restricted to Omega_d is different for k≠0:
face 0 gets multiplied by ω^{k·s_1}.

For computing Betti over Q (or C), we work over C with the ω factor.

Let me verify this for P_7.
"""
import numpy as np
from math import pi, comb

def get_QR(p):
    QR = set()
    for x in range(1, p):
        QR.add(pow(x, 2, p))
    return QR

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

p = 7
m = 3
QR = get_QR(p)
omega = np.exp(2j * pi / p)

# Build diff-seqs
diffseqs = {}
for d in range(2*m + 1):
    diffseqs[d] = build_diffseqs(p, d)
    print(f"d={d}: |A_d| = {len(diffseqs[d])}")

# Compute Omega bases (same for all k)
def compute_omega_basis(d, Ad):
    if d <= 1:
        return np.eye(len(Ad), dtype=complex)
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
        return np.eye(len(Ad), dtype=complex)
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
    U, S, Vt = np.linalg.svd(C)
    tol = 1e-8
    null_dim = len(Ad) - np.sum(S > tol)
    if null_dim == 0:
        return np.zeros((len(Ad), 0), dtype=complex)
    return Vt[-null_dim:].T

omega_bases = {}
for d in range(2*m + 1):
    omega_bases[d] = compute_omega_basis(d, diffseqs[d])
    print(f"  Omega_{d} = {omega_bases[d].shape[1]}")

# Compute boundary in k-eigenspace
def boundary_k(d, k, Ad_from, Ad_to, omega_from, omega_to):
    """Boundary map in eigenspace k, restricted to Omega."""
    to_idx = {seq: i for i, seq in enumerate(Ad_to)}

    # Full boundary matrix
    B = np.zeros((len(Ad_to), len(Ad_from)), dtype=complex)
    for j, seq in enumerate(Ad_from):
        for i in range(d + 1):
            if i == 0:
                face = seq[1:]
                phase = omega ** (k * seq[0])  # face 0 phase
            elif i == d:
                face = seq[:-1]
                phase = 1.0
            else:
                face = seq[:i-1] + ((seq[i-1] + seq[i]) % p,) + seq[i+1:]
                phase = 1.0

            if face in to_idx:
                B[to_idx[face], j] += (-1)**i * phase

    # Restrict to Omega
    B_restricted = omega_to.conj().T @ B @ omega_from
    return np.linalg.matrix_rank(B_restricted, tol=1e-8)

# Compute Betti for each eigenspace k
print(f"\nPer-eigenspace Betti for P_{p}:")
for k in range(p):
    R = {0: 0}
    for d in range(1, 2*m + 1):
        R[d] = boundary_k(d, k, diffseqs[d], diffseqs[d-1],
                          omega_bases[d], omega_bases[d-1])
    R[2*m + 1] = 0

    betti = []
    for d in range(2*m + 1):
        Omega_d = omega_bases[d].shape[1]
        beta = Omega_d - R[d] - R[d+1]
        betti.append(beta)

    print(f"  k={k}: R = {[R[d] for d in range(2*m+2)]}, beta = {betti}")

# Total Betti
print(f"\nTotal Betti:")
total_betti = [0] * (2*m + 1)
for k in range(p):
    R = {0: 0}
    for d in range(1, 2*m + 1):
        R[d] = boundary_k(d, k, diffseqs[d], diffseqs[d-1],
                          omega_bases[d], omega_bases[d-1])
    R[2*m + 1] = 0
    for d in range(2*m + 1):
        Omega_d = omega_bases[d].shape[1]
        total_betti[d] += Omega_d - R[d] - R[d+1]

print(f"  beta = {total_betti}")
print(f"  chi = {sum((-1)**d * b for d, b in enumerate(total_betti))}")
