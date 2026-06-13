"""
Compute Omega dimensions for P_19 (m=9) through all degrees using mod-3 GF arithmetic.
"""
import numpy as np
import sys

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

def mod_rank(matrix, prime=3):
    """GF(prime) Gaussian elimination."""
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

def build_constraint_sparse(p, d, Ad):
    """Build constraint matrix C_d as int8 array."""
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

p = 19
m = 9
QR = get_QR(p)
print(f"P_{p} (m={m}), QR = {sorted(QR)}")

# Build diff-seqs and compute Omega for each degree
Omega = {}
A_sizes = {}
for d in range(2*m + 1):
    print(f"\n--- d={d} ---", flush=True)
    Ad = build_diffseqs(p, d)
    A_sizes[d] = len(Ad)
    print(f"  |A_d| = {len(Ad)}", flush=True)

    if d <= 1:
        Omega[d] = len(Ad)
        print(f"  Omega_{d} = {Omega[d]} (no constraints)")
        continue

    # Check memory: constraint matrix size
    # For large d, the constraint matrix could be huge
    C = build_constraint_sparse(p, d, Ad)
    mem_mb = C.shape[0] * C.shape[1] / 1e6
    print(f"  C shape: {C.shape}, ~{mem_mb:.1f} MB", flush=True)

    if mem_mb > 2000:
        print(f"  TOO LARGE, skipping", flush=True)
        Omega[d] = -1
        continue

    r = mod_rank(C) if C.shape[0] > 0 else 0
    Omega[d] = len(Ad) - r
    print(f"  rank(C) = {r}, Omega_{d} = {Omega[d]}", flush=True)

print(f"\n=== SUMMARY ===")
print(f"P_{p} (m={m}):")
print(f"  |A_d| = {[A_sizes.get(d, '?') for d in range(2*m+1)]}")
print(f"  Omega = {[Omega.get(d, '?') for d in range(2*m+1)]}")
print(f"  chi = {sum((-1)**d * Omega[d] for d in range(2*m+1) if Omega.get(d, -1) >= 0)}")
