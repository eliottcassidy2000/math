"""
Compute P_11 k=0 eigenspace Betti numbers using the stacking trick:
  rank(∂_d|Omega_d) = rank([C_d; B_d]) - rank(C_d)

where C_d is the constraint matrix and B_d is the boundary matrix.
This avoids computing the null space of C_d explicitly.

Uses mod-3 Gaussian elimination for speed.

opus-2026-03-13-S71b
"""
import numpy as np
import sys

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

def mod_rank(matrix, prime=3):
    """GF(prime) Gaussian elimination for column rank."""
    M = matrix.copy().astype(np.int64) % prime
    rows, cols = M.shape
    rank = 0
    for col in range(cols):
        pivot = None
        for row in range(rank, rows):
            if M[row, col] % prime != 0:
                pivot = row
                break
        if pivot is None:
            continue
        M[[rank, pivot]] = M[[pivot, rank]]
        inv = pow(int(M[rank, col]), prime - 2, prime)
        M[rank] = (M[rank] * inv) % prime
        for row in range(rows):
            if row != rank and M[row, col] % prime != 0:
                M[row] = (M[row] - M[row, col] * M[rank]) % prime
        rank += 1
    return rank

p = 11
m = 5
QR = get_QR(p)

# Build diff-seqs for all degrees
print(f"P_{p} (m={m}): building diff-seqs...", flush=True)
diffseqs = {}
for d in range(2*m + 1):
    diffseqs[d] = build_diffseqs(p, d)
    print(f"  d={d}: |A_d| = {len(diffseqs[d])}", flush=True)

# Build constraint matrices
def build_constraint(d, Ad):
    """Build constraint matrix C_d."""
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
        return np.zeros((0, len(Ad)), dtype=np.int8)

    C = np.zeros((n_rows, len(Ad)), dtype=np.int8)
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

def build_boundary(d, Ad_from, Ad_to):
    """Build boundary matrix ∂_d: A_d → A_{d-1}."""
    to_idx = {seq: i for i, seq in enumerate(Ad_to)}
    B = np.zeros((len(Ad_to), len(Ad_from)), dtype=np.int8)

    for j, seq in enumerate(Ad_from):
        for i in range(d + 1):
            if i == 0:
                face = seq[1:]
            elif i == d:
                face = seq[:-1]
            else:
                face = seq[:i-1] + ((seq[i-1] + seq[i]) % p,) + seq[i+1:]

            if face in to_idx:
                B[to_idx[face], j] += (-1) ** i

    return B

# Compute Omega dimensions (= |A_d| - rank(C_d))
print("\nComputing Omega dimensions...", flush=True)
Omega = {}
C_rank = {}
for d in range(2*m + 1):
    if d <= 1:
        Omega[d] = len(diffseqs[d])
        C_rank[d] = 0
    else:
        C = build_constraint(d, diffseqs[d])
        r = mod_rank(C) if C.shape[0] > 0 else 0
        C_rank[d] = r
        Omega[d] = len(diffseqs[d]) - r
    print(f"  Omega_{d} = {Omega[d]}", flush=True)

# Compute restricted boundary ranks using stacking trick
# rank(∂_d|Omega_d) = rank([C_d; B_d]) - rank(C_d)
# where B_d sends to A_{d-1} coordinates and C_d acts on A_d columns
print("\nComputing restricted boundary ranks...", flush=True)
R = {0: 0}
for d in range(1, 2*m + 1):
    print(f"  Computing R_{d}...", flush=True)

    Ad_from = diffseqs[d]
    Ad_to = diffseqs[d-1]

    # Constraint matrix for degree d
    C_d = build_constraint(d, Ad_from)

    # Boundary matrix ∂_d: A_d → A_{d-1}
    B_d = build_boundary(d, Ad_from, Ad_to)

    # But we need the boundary in terms of the A_{d-1} basis
    # The boundary of an element of ker(C_d) lands in ker(C_{d-1})
    # So ∂_d|Omega maps to Omega_{d-1}

    # Stack: [C_d; B_d] has shape (junk_rows + |A_{d-1}|) × |A_d|
    stacked = np.vstack([C_d, B_d]) if C_d.shape[0] > 0 else B_d

    # Check memory
    mem_mb = stacked.shape[0] * stacked.shape[1] * 8 / 1e6
    if mem_mb > 4000:
        print(f"    Matrix {stacked.shape} = {mem_mb:.0f} MB — TOO LARGE, skipping")
        R[d] = -1
        continue

    rank_stacked = mod_rank(stacked)
    rank_C = C_rank[d]
    R[d] = rank_stacked - rank_C
    print(f"    rank([C; B]) = {rank_stacked}, rank(C) = {rank_C}, R_{d} = {R[d]}")

# Compute Betti numbers
print(f"\nBetti numbers (k=0 eigenspace):")
betti = []
for d in range(2*m + 1):
    K_d = Omega[d] - R.get(d, 0)
    R_next = R.get(d+1, 0)
    if R.get(d, -1) == -1 or R.get(d+1, -1) == -1:
        beta_d = "?"
    else:
        beta_d = K_d - R_next
    betti.append(beta_d)
    print(f"  beta_{d} = Omega_{d} - R_{d} - R_{d+1} = {Omega[d]} - {R.get(d, '?')} - {R.get(d+1, '?')} = {beta_d}")

print(f"\n  beta^(0) = {betti}")

# Check known values
print(f"\n  Expected: beta_5^(0) = 5, beta_6^(0) = 5")
print(f"  For m>=5: beta_m^(0) = m(m-3)/2 = {m*(m-3)//2}")
