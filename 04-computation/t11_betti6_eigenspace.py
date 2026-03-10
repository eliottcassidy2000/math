"""
t11_betti6_eigenspace.py - Compute beta_6 for each eigenspace of T_11

Key question: Is beta_6^(k) = 1 for each non-trivial eigenspace k=1,...,10?

Method:
1. Compute Omega_5^(k) basis (460-dimensional) for each k
2. Compute Omega_6^(k) basis (predicted 460-dimensional) for each k
3. Build boundary d_6: Omega_6^(k) -> Omega_5^(k) as a 460x460 matrix
4. rank(d_6^(k)) = 460 - ker(d_6^(k))
   If rank = 459: ker = 1 => beta_6^(k) >= 1 (need im(d_7) to get exact value)
   If rank = 460: ker = 0 => beta_6^(k) = 0

Also check: Omega_6^(k) dimension (should = 460 = Omega_5^(k) by palindrome)

Author: kind-pasteur-2026-03-10-S50
"""
import sys
import time
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import RANK_PRIME, _gauss_rank_np, _gauss_nullbasis_modp, matmul_mod

PRIME = RANK_PRIME


def qr_set(p):
    return set((a * a) % p for a in range(1, p))


def find_pth_root_of_unity(p, prime):
    exp = (prime - 1) // p
    for g in range(2, 30):
        omega = pow(g, exp, prime)
        if omega == 1:
            continue
        if all(pow(omega, k, prime) != 1 for k in range(1, p)):
            return omega
    return None


def enumerate_diff_seqs_to_deg(S, n, max_deg):
    """Enumerate valid diff seqs with no repeated vertices."""
    S_sorted = sorted(S)
    seqs = {0: [()]}
    partial_sums = {(): frozenset([0])}
    partial_last = {(): 0}

    for m in range(1, max_deg + 1):
        prev = seqs[m - 1]
        new = []
        new_ps = {}
        new_last = {}
        for seq in prev:
            ps = partial_sums[seq]
            last = partial_last[seq]
            for s in S_sorted:
                nls = (last + s) % n
                if nls not in ps:
                    ns = seq + (s,)
                    new.append(ns)
                    new_ps[ns] = ps | frozenset([nls])
                    new_last[ns] = nls
        seqs[m] = new
        partial_sums.update(new_ps)
        partial_last.update(new_last)

    return seqs


def compute_face_diff(D, face_idx, n):
    m = len(D)
    if face_idx == 0:
        return D[1:], D[0]
    elif face_idx == m:
        return D[:m - 1], 0
    else:
        merged = (D[face_idx - 1] + D[face_idx]) % n
        face = D[:face_idx - 1] + (merged,) + D[face_idx + 1:]
        return face, 0


def gauss_rank_mod(C, prime):
    """Gaussian elimination to find rank mod prime."""
    C = C.copy() % prime
    rows, cols = C.shape
    rank = 0
    pivot_row = 0

    for col in range(cols):
        found = -1
        for row in range(pivot_row, rows):
            if C[row, col] != 0:
                found = row
                break
        if found == -1:
            continue
        C[[pivot_row, found]] = C[[found, pivot_row]]
        inv = pow(int(C[pivot_row, col]), prime - 2, prime)
        C[pivot_row] = (C[pivot_row] * inv) % prime
        for row in range(rows):
            if row != pivot_row and C[row, col] != 0:
                factor = C[row, col]
                C[row] = (C[row] - factor * C[pivot_row]) % prime
        rank += 1
        pivot_row += 1

    return rank


def gauss_nullbasis_mod(C, prime):
    """Return null space basis as rows of a matrix, mod prime."""
    C = C.copy() % prime
    rows, cols = C.shape
    # Augment with identity
    aug = np.zeros((rows, cols + cols), dtype=np.int64)
    aug[:, :cols] = C
    for i in range(rows):
        pass  # no augment needed for null space

    # We want null space of C: Cx=0
    # Use reduced row echelon form, then extract free variables
    pivot_cols = []
    pivot_row = 0
    for col in range(cols):
        found = -1
        for row in range(pivot_row, rows):
            if C[row, col] != 0:
                found = row
                break
        if found == -1:
            continue
        C[[pivot_row, found]] = C[[found, pivot_row]]
        inv = pow(int(C[pivot_row, col]), prime - 2, prime)
        C[pivot_row] = (C[pivot_row] * inv) % prime
        for row in range(rows):
            if row != pivot_row and C[row, col] != 0:
                factor = C[row, col]
                C[row] = (C[row] - factor * C[pivot_row]) % prime
        pivot_cols.append(col)
        pivot_row += 1

    free_cols = [c for c in range(cols) if c not in set(pivot_cols)]
    n_free = len(free_cols)
    if n_free == 0:
        return np.zeros((0, cols), dtype=np.int64)

    basis = np.zeros((n_free, cols), dtype=np.int64)
    for i, fc in enumerate(free_cols):
        basis[i, fc] = 1
        for j, pc in enumerate(pivot_cols):
            if j < C.shape[0]:
                basis[i, pc] = (-C[j, fc]) % prime

    return basis


def compute_omega_basis_and_boundary(S, n, deg_lo, deg_hi, omega_p, k):
    """
    Compute Omega bases and boundary map d_{deg_hi}: Omega_{deg_hi} -> Omega_{deg_lo}
    for a single eigenspace k.

    Returns: (basis_lo, basis_hi, rank_d_hi)
    where basis_lo has shape (Omega_lo_dim, |A_lo|)
    and basis_hi has shape (Omega_hi_dim, |A_hi|)
    """
    omega_k = pow(omega_p, k, PRIME)
    diff_seqs = enumerate_diff_seqs_to_deg(S, n, deg_hi + 1)

    allowed = {}
    for m in range(deg_hi + 2):
        allowed[m] = set(diff_seqs.get(m, []))

    def build_omega_basis(m):
        A_m = diff_seqs[m]
        n_Am = len(A_m)
        A_m_idx = {d: i for i, d in enumerate(A_m)}

        if m == 0:
            return np.ones((1, 1), dtype=np.int64)

        # Find all junk faces
        junk = set()
        face_data_m = []
        for D in A_m:
            faces = []
            for fi in range(m + 1):
                fd, offset = compute_face_diff(D, fi, n)
                sign = 1 if fi % 2 == 0 else -1
                is_allowed = (fd in allowed[m - 1])
                faces.append((fd, sign, offset, is_allowed))
                if not is_allowed:
                    junk.add(fd)
            face_data_m.append(faces)

        junk_list = sorted(junk)
        n_junk = len(junk_list)

        if n_junk == 0:
            return np.eye(n_Am, dtype=np.int64)

        junk_idx = {j: i for i, j in enumerate(junk_list)}
        C = np.zeros((n_junk, n_Am), dtype=np.int64)
        for j, (D, faces) in enumerate(zip(A_m, face_data_m)):
            for fd, sign, offset, is_allowed in faces:
                if not is_allowed:
                    row = junk_idx[fd]
                    w = pow(omega_k, offset, PRIME) if offset != 0 else 1
                    entry = (sign * w) % PRIME
                    if entry < 0:
                        entry += PRIME
                    C[row, j] = (C[row, j] + entry) % PRIME

        # Null space = Omega basis
        basis = gauss_nullbasis_mod(C, PRIME)
        return basis, face_data_m

    # Build Omega_lo basis
    result_lo = build_omega_basis(deg_lo)
    if isinstance(result_lo, tuple):
        basis_lo, _ = result_lo
    else:
        basis_lo = result_lo

    # Build Omega_hi basis and face data
    result_hi = build_omega_basis(deg_hi)
    if isinstance(result_hi, tuple):
        basis_hi, face_data_hi = result_hi
    else:
        basis_hi = result_hi
        # Regenerate face_data_hi
        A_hi = diff_seqs[deg_hi]
        face_data_hi = []
        for D in A_hi:
            faces = []
            for fi in range(deg_hi + 1):
                fd, offset = compute_face_diff(D, fi, n)
                sign = 1 if fi % 2 == 0 else -1
                is_allowed = (fd in allowed[deg_hi - 1])
                faces.append((fd, sign, offset, is_allowed))
            face_data_hi.append(faces)

    # Build boundary map d_hi: A_hi -> A_lo (allowed faces only)
    A_lo = diff_seqs[deg_lo]
    A_hi = diff_seqs[deg_hi]
    A_lo_idx = {d: i for i, d in enumerate(A_lo)}

    dim_lo = len(A_lo)
    dim_hi = len(A_hi)

    BD = np.zeros((dim_lo, dim_hi), dtype=np.int64)
    for j, (D, faces) in enumerate(zip(A_hi, face_data_hi)):
        for fd, sign, offset, is_allowed in faces:
            if is_allowed and fd in A_lo_idx:
                row = A_lo_idx[fd]
                w = pow(omega_k, offset, PRIME) if offset != 0 else 1
                entry = (sign * w) % PRIME
                if entry < 0:
                    entry += PRIME
                BD[row, j] = (BD[row, j] + entry) % PRIME

    # Boundary restricted to Omega_hi, projected to Omega_lo:
    # d^(k)_hi = basis_lo^{-1} * BD * basis_hi^T
    # (left multiply by basis_lo from left using least-squares in Omega_lo coords)

    # First: d_on_Omega_hi = BD @ basis_hi.T (express output in A_lo coords)
    if basis_hi.shape[0] == 0 or basis_lo.shape[0] == 0:
        return basis_lo, basis_hi, 0

    # BD: (dim_lo, dim_hi)
    # basis_hi: (Omega_hi_dim, dim_hi)
    # d_on_A = BD @ basis_hi.T: (dim_lo, Omega_hi_dim)
    d_on_A = matmul_mod(BD, basis_hi.T, PRIME)

    # Now express d_on_A in Omega_lo coords: solve basis_lo.T @ x = d_on_A
    # basis_lo: (Omega_lo_dim, dim_lo)
    # We want: basis_lo @ d_on_A = boundary in Omega_lo coords
    # This is: d_in_Omega = basis_lo @ d_on_A: (Omega_lo_dim, Omega_hi_dim)
    d_in_Omega = matmul_mod(basis_lo, d_on_A, PRIME)

    # Rank of this matrix
    rk = int(_gauss_rank_np(d_in_Omega.copy(), PRIME))

    return basis_lo, basis_hi, rk


def main():
    print("=" * 70)
    print("T_11 BETA_6 EIGENSPACE COMPUTATION")
    print("=" * 70)

    p = 11
    S = qr_set(p)
    omega_p = find_pth_root_of_unity(p, PRIME)

    print(f"p={p}, QR={sorted(S)}")
    print(f"omega_p = {omega_p}")
    print()
    print("Computing boundary d_6: Omega_6 -> Omega_5 for each eigenspace k...")
    print("Expected: Omega_6^(k) = 460 = Omega_5^(k) (palindrome prediction)")
    print("Prediction: rank(d_6^(k)) = 459 for k!=0, 460 for k=0 (contractible)")
    print()

    # Pre-enumerate (shared across eigenspaces)
    t_pre = time.time()
    diff_seqs = enumerate_diff_seqs_to_deg(S, p, 7)
    print(f"Diff seqs enumerated in {time.time()-t_pre:.1f}s")
    for m in range(8):
        print(f"  deg {m}: {len(diff_seqs.get(m, []))} seqs")
    print()

    results = {}
    t0 = time.time()

    for k in range(p):
        t_k = time.time()
        basis_lo, basis_hi, rk_d6 = compute_omega_basis_and_boundary(
            S, p, 5, 6, omega_p, k
        )
        omega5_dim = basis_lo.shape[0] if basis_lo is not None else 0
        omega6_dim = basis_hi.shape[0] if basis_hi is not None else 0
        ker_d6 = omega6_dim - rk_d6

        results[k] = {
            'omega5': omega5_dim,
            'omega6': omega6_dim,
            'rk_d6': rk_d6,
            'ker_d6': ker_d6,
        }
        elapsed_k = time.time() - t_k
        elapsed = time.time() - t0
        print(f"  k={k}: Omega5={omega5_dim}, Omega6={omega6_dim}, "
              f"rk(d6)={rk_d6}, ker(d6)={ker_d6}  "
              f"[{elapsed_k:.1f}s / {elapsed:.1f}s total]")

    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print("\nBoundary ranks:")
    print(f"{'k':>4} {'Omega5':>8} {'Omega6':>8} {'rk(d6)':>8} {'ker(d6)':>8}")
    for k in range(p):
        r = results[k]
        print(f"  {k:2d} {r['omega5']:8d} {r['omega6']:8d} {r['rk_d6']:8d} {r['ker_d6']:8d}")

    # Total ker(d6) across all eigenspaces = dim of H_6 components before subtracting im(d7)
    total_ker_d6 = sum(results[k]['ker_d6'] for k in range(p))
    print(f"\nTotal ker(d6) = {total_ker_d6} (= sum of all eigenspace ker dims)")
    print(f"If all im(d7) = 0: total beta_6 = {total_ker_d6}")
    print(f"Prediction (HYP-433): beta_6 = 10 = p-1")

    # Per-eigenspace
    ker_by_k = [results[k]['ker_d6'] for k in range(p)]
    print(f"\nPer-eigenspace ker(d6): {ker_by_k}")
    nonzero_k = [k for k in range(p) if results[k]['ker_d6'] > 0]
    print(f"Eigenspaces with ker(d6)>0: {nonzero_k}")

    # Palindrome check
    omega5s = [results[k]['omega5'] for k in range(p)]
    omega6s = [results[k]['omega6'] for k in range(p)]
    omega5_eq = all(d == omega5s[0] for d in omega5s)
    omega6_eq = all(d == omega6s[0] for d in omega6s)
    palindrome_ok = (omega5s[0] == omega6s[0])
    print(f"\nAll eigenspaces equal for Omega5? {omega5_eq} (dim={omega5s[0]})")
    print(f"All eigenspaces equal for Omega6? {omega6_eq} (dim={omega6s[0]})")
    print(f"Palindrome Omega5=Omega6? {palindrome_ok}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
