"""
t11_higher_betti.py - Verify beta_7 through beta_10 = 0 for T_11

Uses the chain complex structure: compute rk(d_m) for m=8,9,10
to check exactness:
  ker(d_7) = 300 (known) -> need rk(d_8) = 300
  ker(d_8) = 150 (if rk=300) -> need rk(d_9) = 150
  ker(d_9) = 30 (if rk=150) -> need rk(d_10) = 30 (full rank, injective)

Per-eigenspace Omega dims: [1,5,20,70,205,460,700,690,450,180,30]

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

    junk_hi, face_data_hi = build_face_data_and_junk(A_hi, allowed[deg_hi - 1], n)
    if junk_hi:
        C_hi = build_constraint_matrix(A_hi, face_data_hi, junk_hi, omega_k, prime)
        _, basis_hi = gauss_rank_nullbasis_mod(C_hi, prime)
    else:
        basis_hi = np.eye(len(A_hi), dtype=np.int16)

    junk_lo, face_data_lo = build_face_data_and_junk(A_lo, allowed[deg_lo - 1], n)
    if junk_lo:
        C_lo = build_constraint_matrix(A_lo, face_data_lo, junk_lo, omega_k, prime)
        _, basis_lo = gauss_rank_nullbasis_mod(C_lo, prime)
    else:
        basis_lo = np.eye(len(A_lo), dtype=np.int16)

    omega_hi_dim = basis_hi.shape[0]
    omega_lo_dim = basis_lo.shape[0]

    if omega_hi_dim == 0 or omega_lo_dim == 0:
        return omega_lo_dim, omega_hi_dim, 0

    BD = build_boundary_matrix(A_hi, face_data_hi, A_lo_idx, omega_k, prime)
    d_on_A = matmul_mod_small(BD, basis_hi.T, prime)
    d_in_Omega = matmul_mod_small(basis_lo, d_on_A, prime)

    rk = gauss_rank_mod(d_in_Omega.astype(np.uint8), prime)
    return omega_lo_dim, omega_hi_dim, rk


def main():
    print("=" * 70)
    print("T_11 HIGHER BETTI VERIFICATION (beta_7 through beta_10 = 0?)")
    print("=" * 70)

    p = 11
    S = qr_set(p)
    prime = PRIME

    omega_p = find_pth_root_of_unity(p, prime)
    print(f"p={p}, QR={sorted(S)}, PRIME={prime}, omega_p={omega_p}")
    print()

    # Known per-eigenspace Omega dims
    omega_dims = [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]
    print(f"Per-eigenspace Omega dims: {omega_dims}")
    print()

    # Known rk(d_7) = 390 for all k -> ker(d_7) = 690 - 390 = 300
    rk_d7 = 390
    ker_d7 = omega_dims[7] - rk_d7  # = 300

    print(f"Known: rk(d_7)=390 for all k -> ker(d_7)={ker_d7}")
    print()
    print("For beta_7=0 need: rk(d_8) = 300")
    print("For beta_8=0 need: rk(d_9) = 150  (= 450-300)")
    print("For beta_9=0 need: rk(d_10) = 30   (= 180-150)")
    print("For beta_10=0: ker(d_10) = 30-rk(d_10) and im(d_11)=0 -> beta_10 = 30-rk(d_10)")
    print("  => need rk(d_10) = 30 (injective)")
    print()

    t0 = time.time()
    diff_seqs = enumerate_diff_seqs(S, p, 11)
    print(f"Enumerated diff-seqs in {time.time()-t0:.1f}s")
    allowed = {m: set(diff_seqs.get(m, [])) for m in range(12)}
    print()

    # Check k=0 and k=1 (by symmetry, all k should give same result)
    for k in [0, 1]:
        omega_k = pow(omega_p, k, prime)
        print(f"--- Eigenspace k={k} ---")
        t_k = time.time()

        # rk(d_8): d_8 maps Omega_8 -> Omega_7
        lo_dim, hi_dim, rk_d8 = compute_boundary_rank(
            diff_seqs, allowed, 8, 7, omega_k, prime, p)
        ker_d8 = hi_dim - rk_d8
        beta7 = ker_d7 - rk_d8
        print(f"  rk(d8): Omega7={lo_dim} Omega8={hi_dim} rk={rk_d8} ker(d8)={ker_d8} "
              f"beta_7={beta7} [{time.time()-t_k:.1f}s]")

        t_m = time.time()
        # rk(d_9): d_9 maps Omega_9 -> Omega_8
        lo_dim9, hi_dim9, rk_d9 = compute_boundary_rank(
            diff_seqs, allowed, 9, 8, omega_k, prime, p)
        ker_d9 = hi_dim9 - rk_d9
        beta8 = ker_d8 - rk_d9
        print(f"  rk(d9): Omega8={lo_dim9} Omega9={hi_dim9} rk={rk_d9} ker(d9)={ker_d9} "
              f"beta_8={beta8} [{time.time()-t_m:.1f}s]")

        t_m = time.time()
        # rk(d_10): d_10 maps Omega_10 -> Omega_9
        lo_dim10, hi_dim10, rk_d10 = compute_boundary_rank(
            diff_seqs, allowed, 10, 9, omega_k, prime, p)
        ker_d10 = hi_dim10 - rk_d10
        beta9 = ker_d9 - rk_d10
        beta10 = ker_d10  # im(d_11) = 0
        print(f"  rk(d10): Omega9={lo_dim10} Omega10={hi_dim10} rk={rk_d10} ker(d10)={ker_d10} "
              f"beta_9={beta9} beta_10={beta10} [{time.time()-t_m:.1f}s]")

        print(f"  Total time for k={k}: {time.time()-t_k:.1f}s")
        print()

    print()
    print("Expected: all beta_7=beta_8=beta_9=beta_10=0")
    print("=> rk(d8)=300, rk(d9)=150, rk(d10)=30")


if __name__ == '__main__':
    main()
    print("\nDONE.")
