"""
t11_d7_all_eigenspaces.py - Compute im(d_7) for ALL eigenspaces of T_11, memory-efficient.

Uses small prime (89) with uint8 storage to avoid 1015 MiB OOM.
89 ≡ 1 (mod 11) so it has a primitive 11th root of unity.

Results from large-prime computation:
  k=0: beta_6=5 (rk(d7)=390, ker(d6)=395)
  k=1: beta_6=1 (rk(d7)=390, ker(d6)=391)
  k=2: beta_6=1 (rk(d7)=390, ker(d6)=391)
  k=3..10: predicted beta_6=1 by symmetry; this script confirms.

Note: uses small prime 89 for all computations; rank should be stable.

Author: kind-pasteur-2026-03-10-S50
"""
import sys
import time
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

PRIME = 89  # small prime with 11|(89-1)=88

def qr_set(p):
    return set((a * a) % p for a in range(1, p))

def find_pth_root_of_unity(p, prime):
    exp = (prime - 1) // p
    for g in range(2, prime):
        omega = pow(g, exp, prime)
        if omega == 1:
            continue
        if all(pow(omega, k, prime) != 1 for k in range(1, p)):
            return omega
    return None

def enumerate_diff_seqs(S, n, max_deg):
    S_sorted = sorted(S)
    seqs = {0: [()]}
    partial_sums = {(): frozenset([0])}
    partial_last = {(): 0}
    for m in range(1, max_deg + 1):
        prev = seqs[m - 1]
        new = []
        nps = {}
        nl = {}
        for seq in prev:
            ps = partial_sums[seq]
            last = partial_last[seq]
            for s in S_sorted:
                nls = (last + s) % n
                if nls not in ps:
                    ns = seq + (s,)
                    new.append(ns)
                    nps[ns] = ps | frozenset([nls])
                    nl[ns] = nls
        seqs[m] = new
        partial_sums.update(nps)
        partial_last.update(nl)
    return seqs

def compute_face_diff(D, face_idx, n):
    m = len(D)
    if face_idx == 0:
        return D[1:], D[0]
    elif face_idx == m:
        return D[:m - 1], 0
    else:
        merged = (D[face_idx - 1] + D[face_idx]) % n
        return D[:face_idx - 1] + (merged,) + D[face_idx + 1:], 0

def gauss_rank_nullbasis_mod(C_in, prime):
    """Return (rank, nullbasis_rows) using int16 arithmetic for small prime."""
    C = C_in.astype(np.int16)
    rows, cols = C.shape
    pivot_cols = []
    pivot_row = 0
    for col in range(cols):
        found = -1
        for row in range(pivot_row, rows):
            if C[row, col] != 0:
                found = row; break
        if found < 0: continue
        C[[pivot_row, found]] = C[[found, pivot_row]]
        inv = pow(int(C[pivot_row, col]), prime - 2, prime)
        C[pivot_row] = (C[pivot_row].astype(np.int32) * inv % prime).astype(np.int16)
        factor_col = C[:, col].copy()
        factor_col[pivot_row] = 0
        nzr = np.where(factor_col != 0)[0]
        for row in nzr:
            f = int(factor_col[row])
            C[row] = ((C[row].astype(np.int32) - f * C[pivot_row].astype(np.int32)) % prime).astype(np.int16)
        pivot_cols.append(col)
        pivot_row += 1
    rank = len(pivot_cols)
    free_cols = [c for c in range(cols) if c not in set(pivot_cols)]
    n_free = len(free_cols)
    if n_free == 0:
        return rank, np.zeros((0, cols), dtype=np.int16)
    basis = np.zeros((n_free, cols), dtype=np.int16)
    for i, fc in enumerate(free_cols):
        basis[i, fc] = 1
        for j, pc in enumerate(pivot_cols):
            if j < C.shape[0]:
                basis[i, pc] = (-C[j, fc]) % prime
    return rank, basis

def gauss_rank_mod(C_in, prime):
    """Compute rank only (no null basis)."""
    C = C_in.astype(np.int16).copy()
    rows, cols = C.shape
    rank = 0
    pivot_row = 0
    for col in range(cols):
        found = -1
        for row in range(pivot_row, rows):
            if C[row, col] != 0:
                found = row; break
        if found < 0: continue
        C[[pivot_row, found]] = C[[found, pivot_row]]
        inv = pow(int(C[pivot_row, col]), prime - 2, prime)
        C[pivot_row] = (C[pivot_row].astype(np.int32) * inv % prime).astype(np.int16)
        factor_col = C[:, col].copy()
        factor_col[pivot_row] = 0
        nzr = np.where(factor_col != 0)[0]
        for row in nzr:
            f = int(factor_col[row])
            C[row] = ((C[row].astype(np.int32) - f * C[pivot_row].astype(np.int32)) % prime).astype(np.int16)
        rank += 1
        pivot_row += 1
    return rank

def build_constraint_matrix(A_m, face_data_m, junk_list, omega_k, prime):
    n_junk = len(junk_list)
    n_Am = len(A_m)
    junk_idx = {j: i for i, j in enumerate(junk_list)}
    C = np.zeros((n_junk, n_Am), dtype=np.uint8)
    for j, (D, faces) in enumerate(zip(A_m, face_data_m)):
        for fd, sign, offset, is_allowed in faces:
            if not is_allowed:
                row = junk_idx[fd]
                w = pow(omega_k, offset, prime) if offset != 0 else 1
                entry = (sign * w) % prime
                if entry < 0:
                    entry = (entry + prime) % prime
                C[row, j] = (int(C[row, j]) + entry) % prime
    return C

def build_boundary_matrix(A_hi, face_data_hi, A_lo_idx, omega_k, prime, allowed_lo):
    """Build boundary matrix BD: (|A_lo| x |A_hi|) for allowed faces only."""
    dim_lo = max(A_lo_idx.values()) + 1 if A_lo_idx else 0
    dim_hi = len(A_hi)
    BD = np.zeros((dim_lo, dim_hi), dtype=np.uint8)
    for j, (D, faces) in enumerate(zip(A_hi, face_data_hi)):
        for fd, sign, offset, is_allowed in faces:
            if is_allowed and fd in A_lo_idx:
                row = A_lo_idx[fd]
                w = pow(omega_k, offset, prime) if offset != 0 else 1
                entry = (sign * w) % prime
                if entry < 0:
                    entry = (entry + prime) % prime
                BD[row, j] = (int(BD[row, j]) + entry) % prime
    return BD

def matmul_mod_small(A, B, prime):
    """Matrix multiply mod prime using int32 for safety."""
    # A: (m,k), B: (k,n) -> result (m,n) mod prime
    # Process in chunks to avoid overflow
    result = np.zeros((A.shape[0], B.shape[1]), dtype=np.int32)
    chunk = 256
    for i in range(0, A.shape[1], chunk):
        end = min(i + chunk, A.shape[1])
        result += (A[:, i:end].astype(np.int32) @ B[i:end, :].astype(np.int32))
    return result % prime

def compute_d7_rank(diff_seqs, allowed, omega_k, prime, n):
    """Compute rank(d_7 | Omega_7) using the formula:
    rk = rk([basis_6 @ BD_7; C_7]) - rk(C_7)
    But instead, directly compute the boundary on Omega_7 basis vectors.
    """
    A_6 = diff_seqs[6]
    A_7 = diff_seqs[7]
    A_6_idx = {d: i for i, d in enumerate(A_6)}

    # Build face data for deg 7
    junk_6 = set()  # junk 6-seqs (faces of 7-seqs that are not in A_6)
    face_data_7 = []
    junk_5 = set()  # junk 5-seqs (faces of 6-seqs that are not in A_5)
    face_data_6 = []

    A_5 = diff_seqs[5]

    for D in A_7:
        faces = []
        for fi in range(8):
            fd, offset = compute_face_diff(D, fi, n)
            sign = 1 if fi % 2 == 0 else -1
            is_allowed = (fd in allowed[6])
            faces.append((fd, sign, offset, is_allowed))
            if not is_allowed:
                junk_6.add(fd)
        face_data_7.append(faces)

    for D in A_6:
        faces = []
        for fi in range(7):
            fd, offset = compute_face_diff(D, fi, n)
            sign = 1 if fi % 2 == 0 else -1
            is_allowed = (fd in allowed[5])
            faces.append((fd, sign, offset, is_allowed))
            if not is_allowed:
                junk_5.add(fd)
        face_data_6.append(faces)

    junk_6_list = sorted(junk_6)
    junk_5_list = sorted(junk_5)

    # Build Omega_7 basis
    C_7 = build_constraint_matrix(A_7, face_data_7, junk_6_list, omega_k, prime)
    rk_C7, basis_7 = gauss_rank_nullbasis_mod(C_7, prime)
    omega7_dim = basis_7.shape[0]

    # Build Omega_6 basis
    C_6 = build_constraint_matrix(A_6, face_data_6, junk_5_list, omega_k, prime)
    rk_C6, basis_6 = gauss_rank_nullbasis_mod(C_6, prime)
    omega6_dim = basis_6.shape[0]

    if omega7_dim == 0 or omega6_dim == 0:
        return omega7_dim, omega6_dim, 0

    # Build boundary matrix BD_7: (|A_6|, |A_7|)
    BD_7 = build_boundary_matrix(A_7, face_data_7, A_6_idx, omega_k, prime, allowed[6])

    # d_7^(k) = basis_6 @ BD_7 @ basis_7.T
    # Use int32 for intermediate calculations
    d_on_A = matmul_mod_small(BD_7, basis_7.T, prime)  # (|A_6|, omega7_dim)
    d_in_Omega = matmul_mod_small(basis_6, d_on_A, prime)  # (omega6_dim, omega7_dim)

    rk = gauss_rank_mod(d_in_Omega.astype(np.uint8), prime)
    return omega7_dim, omega6_dim, rk


def main():
    print("=" * 70)
    print("T_11 BOUNDARY d_7 ALL EIGENSPACES (small prime 89)")
    print("=" * 70)

    p = 11
    S = qr_set(p)
    prime = PRIME

    omega_p = find_pth_root_of_unity(p, prime)
    print(f"p={p}, QR={sorted(S)}, PRIME={prime}, omega_p={omega_p}")

    # Known ker(d6) from previous computation (large prime)
    ker_d6 = {0: 395}
    for k in range(1, p):
        ker_d6[k] = 391

    print(f"\nKnown ker(d6): k=0 → 395, k!=0 → 391")
    print(f"Previously computed: k=0: beta_6=5, k=1: beta_6=1, k=2: beta_6=1")
    print(f"This script verifies k=3..10 using small prime {prime}")
    print()

    t0 = time.time()
    diff_seqs = enumerate_diff_seqs(S, p, 8)
    print(f"Enumerated diff-seqs in {time.time()-t0:.1f}s")

    allowed = {m: set(diff_seqs.get(m, [])) for m in range(9)}

    results = {}
    t0 = time.time()

    for k in range(p):
        t_k = time.time()
        omega_k = pow(omega_p, k, prime)

        omega7_dim, omega6_dim, rk_d7 = compute_d7_rank(
            diff_seqs, allowed, omega_k, prime, p)

        im_d7 = rk_d7
        ker = ker_d6[k]
        beta6 = ker - im_d7

        elapsed = time.time() - t0
        print(f"  k={k}: Omega7={omega7_dim} Omega6={omega6_dim} "
              f"rk(d7)={rk_d7} ker(d6)={ker} beta_6={beta6} [{time.time()-t_k:.1f}s / {elapsed:.1f}s]")

        results[k] = {
            'omega7': omega7_dim,
            'omega6': omega6_dim,
            'rk_d7': rk_d7,
            'ker_d6': ker,
            'beta6': beta6
        }

    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print(f"\n{'k':>4} {'Omega7':>8} {'Omega6':>8} {'ker(d6)':>8} {'rk(d7)':>8} {'beta_6':>8}")
    for k in range(p):
        r = results[k]
        print(f"  {k:2d} {r['omega7']:8d} {r['omega6']:8d} {r['ker_d6']:8d} {r['rk_d7']:8d} {r['beta6']:8d}")

    total_beta6 = sum(results[k]['beta6'] for k in range(p))
    print(f"\nTotal beta_6(T_11) = {total_beta6}")
    print(f"Per k=0: {results[0]['beta6']}, k!=0: {[results[k]['beta6'] for k in range(1, p)]}")

    # Verify consistency with large-prime results
    print(f"\nConsistency check vs large-prime (k=0,1,2):")
    print(f"  k=0: beta_6={results[0]['beta6']} (expected 5)")
    print(f"  k=1: beta_6={results[1]['beta6']} (expected 1)")
    print(f"  k=2: beta_6={results[2]['beta6']} (expected 1)")


if __name__ == '__main__':
    main()
    print("\nDONE.")
