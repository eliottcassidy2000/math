"""
t11_beta5_verify.py - Directly verify beta_5(T_11) = 5

Key prediction (HYP-443):
  beta_5^(k=0) = 5, beta_5^(k≠0) = 0 => total beta_5 = 5.

Method: compute rk(d_5^(k)) and rk(d_6^(k)) for k=0,1.
  ker(d_5^(k)) = Omega_5^(k) - rk(d_5^(k))
  beta_5^(k) = ker(d_5^(k)) - rk(d_6^(k))

Known: rk(d_6^(k=0))=305, rk(d_6^(k≠0))=309.
Expected: rk(d_5^(k=0))=150, rk(d_5^(k≠0))=151.

Uses small prime 89 for memory efficiency.

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
    C = C_in.astype(np.int16)
    rows, cols = C.shape
    pivot_cols = []
    pivot_row = 0
    for col in range(cols):
        found = -1
        for row in range(pivot_row, rows):
            if C[row, col] != 0:
                found = row
                break
        if found < 0:
            continue
        C[[pivot_row, found]] = C[[found, pivot_row]]
        inv = pow(int(C[pivot_row, col]), prime - 2, prime)
        C[pivot_row] = (C[pivot_row].astype(np.int32) * inv % prime).astype(np.int16)
        fc = C[:, col].copy()
        fc[pivot_row] = 0
        nzr = np.where(fc != 0)[0]
        for row in nzr:
            f = int(fc[row])
            C[row] = ((C[row].astype(np.int32) - f * C[pivot_row].astype(np.int32)) % prime).astype(np.int16)
        pivot_cols.append(col)
        pivot_row += 1
    rank = len(pivot_cols)
    free_cols = [c for c in range(cols) if c not in set(pivot_cols)]
    n_free = len(free_cols)
    if n_free == 0:
        return rank, np.zeros((0, cols), dtype=np.int16)
    basis = np.zeros((n_free, cols), dtype=np.int16)
    for i, fc_idx in enumerate(free_cols):
        basis[i, fc_idx] = 1
        for j, pc in enumerate(pivot_cols):
            if j < C.shape[0]:
                basis[i, pc] = (-C[j, fc_idx]) % prime
    return rank, basis


def gauss_rank_mod(C_in, prime):
    C = C_in.astype(np.int16).copy()
    rows, cols = C.shape
    rank = 0
    pivot_row = 0
    for col in range(cols):
        found = -1
        for row in range(pivot_row, rows):
            if C[row, col] != 0:
                found = row
                break
        if found < 0:
            continue
        C[[pivot_row, found]] = C[[found, pivot_row]]
        inv = pow(int(C[pivot_row, col]), prime - 2, prime)
        C[pivot_row] = (C[pivot_row].astype(np.int32) * inv % prime).astype(np.int16)
        fc = C[:, col].copy()
        fc[pivot_row] = 0
        nzr = np.where(fc != 0)[0]
        for row in nzr:
            f = int(fc[row])
            C[row] = ((C[row].astype(np.int32) - f * C[pivot_row].astype(np.int32)) % prime).astype(np.int16)
        rank += 1
        pivot_row += 1
    return rank


def matmul_mod_small(A, B, prime):
    result = np.zeros((A.shape[0], B.shape[1]), dtype=np.int32)
    chunk = 256
    for i in range(0, A.shape[1], chunk):
        end = min(i + chunk, A.shape[1])
        result += (A[:, i:end].astype(np.int32) @ B[i:end, :].astype(np.int32))
    return (result % prime).astype(np.int16)


def build_face_data_and_junk(A_m, allowed_lower, n):
    junk = set()
    face_data = []
    for D in A_m:
        faces = []
        for fi in range(len(D) + 1):
            fd, offset = compute_face_diff(D, fi, n)
            sign = 1 if fi % 2 == 0 else -1
            is_allowed = (fd in allowed_lower)
            faces.append((fd, sign, offset, is_allowed))
            if not is_allowed:
                junk.add(fd)
        face_data.append(faces)
    return sorted(junk), face_data


def build_constraint_matrix(A_m, face_data, junk_list, omega_k, prime):
    n_junk = len(junk_list)
    n_Am = len(A_m)
    junk_idx = {j: i for i, j in enumerate(junk_list)}
    C = np.zeros((n_junk, n_Am), dtype=np.uint8)
    for j, (D, faces) in enumerate(zip(A_m, face_data)):
        for fd, sign, offset, is_allowed in faces:
            if not is_allowed:
                row = junk_idx[fd]
                w = pow(omega_k, offset, prime) if offset != 0 else 1
                entry = (sign * w) % prime
                if entry < 0:
                    entry = (entry + prime) % prime
                C[row, j] = (int(C[row, j]) + entry) % prime
    return C


def build_boundary_matrix(A_hi, face_data_hi, A_lo_idx, omega_k, prime):
    dim_lo = len(A_lo_idx)
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


def compute_boundary_rank(diff_seqs, allowed, deg_hi, deg_lo, omega_k, prime, n):
    """Compute rank(d_{deg_hi}: Omega_{deg_hi} -> Omega_{deg_lo})."""
    A_lo = diff_seqs[deg_lo]
    A_hi = diff_seqs[deg_hi]
    A_lo_idx = {d: i for i, d in enumerate(A_lo)}

    # Build Omega_{deg_hi} basis
    junk_hi, face_data_hi = build_face_data_and_junk(A_hi, allowed[deg_hi - 1], n)
    if junk_hi:
        C_hi = build_constraint_matrix(A_hi, face_data_hi, junk_hi, omega_k, prime)
        rk_chi, basis_hi = gauss_rank_nullbasis_mod(C_hi, prime)
    else:
        basis_hi = np.eye(len(A_hi), dtype=np.int16)
        face_data_hi  # reuse

    # Build Omega_{deg_lo} basis
    junk_lo, face_data_lo = build_face_data_and_junk(A_lo, allowed[deg_lo - 1], n)
    if junk_lo:
        C_lo = build_constraint_matrix(A_lo, face_data_lo, junk_lo, omega_k, prime)
        rk_clo, basis_lo = gauss_rank_nullbasis_mod(C_lo, prime)
    else:
        basis_lo = np.eye(len(A_lo), dtype=np.int16)

    omega_hi_dim = basis_hi.shape[0]
    omega_lo_dim = basis_lo.shape[0]

    if omega_hi_dim == 0 or omega_lo_dim == 0:
        return omega_lo_dim, omega_hi_dim, 0

    # Build boundary matrix
    BD = build_boundary_matrix(A_hi, face_data_hi, A_lo_idx, omega_k, prime)

    # d^(k) = basis_lo @ BD @ basis_hi.T
    d_on_A = matmul_mod_small(BD, basis_hi.T, prime)
    d_in_Omega = matmul_mod_small(basis_lo, d_on_A, prime)

    rk = gauss_rank_mod(d_in_Omega.astype(np.uint8), prime)
    return omega_lo_dim, omega_hi_dim, rk


def main():
    print("=" * 70)
    print("T_11 BETA_5 VERIFICATION")
    print("=" * 70)

    p = 11
    S = qr_set(p)
    prime = PRIME

    omega_p = find_pth_root_of_unity(p, prime)
    print(f"p={p}, QR={sorted(S)}, PRIME={prime}, omega_p={omega_p}")
    print()
    print("Prediction (HYP-443):")
    print("  rk(d_5^(k=0))=150 => ker(d_5^(k=0))=310 => beta_5^(k=0)=310-305=5")
    print("  rk(d_5^(k!=0))=151 => ker(d_5^(k!=0))=309 => beta_5^(k!=0)=309-309=0")
    print()

    t0 = time.time()
    diff_seqs = enumerate_diff_seqs(S, p, 7)
    print(f"Enumerated in {time.time()-t0:.1f}s")
    allowed = {m: set(diff_seqs.get(m, [])) for m in range(8)}
    print()

    # Known rk(d_6) from t11_betti6_eigenspace.py
    rk_d6 = {0: 305}
    for k in range(1, p):
        rk_d6[k] = 309

    print(f"Known rk(d_6): k=0 -> 305, k!=0 -> 309")
    print()

    # Compute rk(d_5) for k=0,1 (all eigenspaces should be identical or near-identical)
    print("Computing rk(d_5^(k))...")
    for k in [0, 1]:
        t_k = time.time()
        omega_k = pow(omega_p, k, prime)
        lo_dim, hi_dim, rk_d5 = compute_boundary_rank(
            diff_seqs, allowed, 5, 4, omega_k, prime, p)
        ker_d5 = hi_dim - rk_d5  # ker(d_5) = Omega_5 - rk(d_5)

        # Wait, Omega_5 is hi_dim? No: d_5: Omega_5 -> Omega_4. hi = deg 5, lo = deg 4.
        # lo_dim = Omega_4, hi_dim = Omega_5
        omega4 = lo_dim
        omega5 = hi_dim
        ker_d5 = omega5 - rk_d5

        beta5 = ker_d5 - rk_d6[k]
        elapsed = time.time() - t_k
        print(f"  k={k}: Omega4={omega4} Omega5={omega5} rk(d5)={rk_d5} "
              f"ker(d5)={ker_d5} rk(d6)={rk_d6[k]} beta5={beta5} [{elapsed:.1f}s]")

    print()
    print("Expected: k=0 -> beta5=5, k!=0 -> beta5=0")

    # Also verify by computing rk(d_6) with small prime for cross-check
    print()
    print("Cross-checking rk(d_6) with small prime...")
    for k in [0, 1]:
        t_k = time.time()
        omega_k = pow(omega_p, k, prime)
        lo_dim, hi_dim, rk_d6_check = compute_boundary_rank(
            diff_seqs, allowed, 6, 5, omega_k, prime, p)
        elapsed = time.time() - t_k
        print(f"  k={k}: Omega5={lo_dim} Omega6={hi_dim} rk(d6)={rk_d6_check} "
              f"(expected {rk_d6[k]}) [{elapsed:.1f}s]")


if __name__ == '__main__':
    main()
    print("\nDONE.")
