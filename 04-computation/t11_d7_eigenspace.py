"""
t11_d7_eigenspace.py - Compute im(d_7) for each eigenspace of T_11

Key question: What is rank(d_7^(k)) for each eigenspace k=0,...,10?

This determines beta_6^(k) = ker(d_6^(k)) - im(d_7^(k)):
  ker(d_6^(k=0))  = 395
  ker(d_6^(k!=0)) = 391

If HYP-433 is correct (beta_6=10, one per non-trivial eigenspace):
  im(d_7^(k!=0)) = 390 (rank 390 out of Omega_7 dim)
  im(d_7^(k=0))  = 395 (trivial eigenspace contractible)

Also computes Omega_7^(k) = dim for each eigenspace.

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


def gauss_nullbasis_mod(C, prime):
    """Return null space basis as rows of a matrix, mod prime."""
    C = C.copy() % prime
    rows, cols = C.shape
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


def compute_omega_basis(diff_seqs, allowed, deg, omega_k, n):
    """Compute Omega_deg^(k) basis and face data."""
    A_m = diff_seqs[deg]
    n_Am = len(A_m)

    if deg == 0:
        return np.ones((1, 1), dtype=np.int64), [[((), 1, 0, True)]]

    # Find junk faces
    junk = set()
    face_data = []
    for D in A_m:
        faces = []
        for fi in range(deg + 1):
            fd, offset = compute_face_diff(D, fi, n)
            sign = 1 if fi % 2 == 0 else -1
            is_allowed = (fd in allowed[deg - 1])
            faces.append((fd, sign, offset, is_allowed))
            if not is_allowed:
                junk.add(fd)
        face_data.append(faces)

    junk_list = sorted(junk)
    n_junk = len(junk_list)

    if n_junk == 0:
        return np.eye(n_Am, dtype=np.int64), face_data

    junk_idx = {j: i for i, j in enumerate(junk_list)}
    C = np.zeros((n_junk, n_Am), dtype=np.int64)
    for j, (D, faces) in enumerate(zip(A_m, face_data)):
        for fd, sign, offset, is_allowed in faces:
            if not is_allowed:
                row = junk_idx[fd]
                w = pow(omega_k, offset, PRIME) if offset != 0 else 1
                entry = (sign * w) % PRIME
                if entry < 0:
                    entry += PRIME
                C[row, j] = (C[row, j] + entry) % PRIME

    basis = gauss_nullbasis_mod(C, PRIME)
    return basis, face_data


def compute_boundary_rank(diff_seqs, allowed, deg_hi, deg_lo, omega_k, n,
                          basis_lo=None, basis_hi=None):
    """
    Compute rank of d_{deg_hi}: Omega_{deg_hi} -> Omega_{deg_lo}.
    Returns (basis_lo, basis_hi, rank).
    If bases are already computed, pass them in.
    """
    A_lo = diff_seqs[deg_lo]
    A_hi = diff_seqs[deg_hi]
    A_lo_idx = {d: i for i, d in enumerate(A_lo)}

    dim_lo = len(A_lo)
    dim_hi = len(A_hi)

    if basis_hi is None:
        basis_hi, face_data_hi = compute_omega_basis(
            diff_seqs, allowed, deg_hi, omega_k, n)
    else:
        # Regenerate face data for deg_hi
        face_data_hi = []
        for D in A_hi:
            faces = []
            for fi in range(deg_hi + 1):
                fd, offset = compute_face_diff(D, fi, n)
                sign = 1 if fi % 2 == 0 else -1
                is_allowed = (fd in allowed[deg_hi - 1])
                faces.append((fd, sign, offset, is_allowed))
            face_data_hi.append(faces)

    if basis_lo is None:
        basis_lo, _ = compute_omega_basis(
            diff_seqs, allowed, deg_lo, omega_k, n)

    if basis_hi.shape[0] == 0 or basis_lo.shape[0] == 0:
        return basis_lo, basis_hi, 0

    # Build boundary map d_{deg_hi}: A_hi -> A_lo (allowed faces only)
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

    # d^(k) = basis_lo @ BD @ basis_hi.T
    d_on_A = matmul_mod(BD, basis_hi.T, PRIME)
    d_in_Omega = matmul_mod(basis_lo, d_on_A, PRIME)
    rk = int(_gauss_rank_np(d_in_Omega.copy(), PRIME))

    return basis_lo, basis_hi, rk


def main():
    print("=" * 70)
    print("T_11 BOUNDARY d_7 EIGENSPACE COMPUTATION")
    print("=" * 70)

    p = 11
    S = qr_set(p)
    omega_p = find_pth_root_of_unity(p, PRIME)

    print(f"p={p}, QR={sorted(S)}")
    print(f"omega_p = {omega_p}")
    print()
    print("Key question: What is rank(d_7^(k)) for each eigenspace k?")
    print("This determines beta_6^(k) = ker(d_6^(k)) - im(d_7^(k))")
    print("  ker(d_6^(k=0)) = 395")
    print("  ker(d_6^(k!=0)) = 391")
    print()

    # Pre-enumerate
    t_pre = time.time()
    diff_seqs = enumerate_diff_seqs_to_deg(S, p, 8)
    print(f"Diff seqs enumerated in {time.time()-t_pre:.1f}s")
    for m in range(9):
        print(f"  deg {m}: {len(diff_seqs.get(m, []))} seqs")
    print()

    # Pre-compute allowed sets
    allowed = {}
    for m in range(9):
        allowed[m] = set(diff_seqs.get(m, []))

    # Known ker(d6) from previous computation
    ker_d6 = {0: 395}
    for k in range(1, p):
        ker_d6[k] = 391

    results = {}
    t0 = time.time()

    for k in range(p):
        t_k = time.time()
        omega_k = pow(omega_p, k, PRIME)

        # Compute Omega_6 basis (known: dim=700)
        basis_6, face_data_6 = compute_omega_basis(
            diff_seqs, allowed, 6, omega_k, p)
        omega6_dim = basis_6.shape[0]

        t_after6 = time.time()
        print(f"  k={k}: Omega_6={omega6_dim} ({t_after6-t_k:.1f}s)...", end='', flush=True)

        # Compute Omega_7 basis
        basis_7, face_data_7 = compute_omega_basis(
            diff_seqs, allowed, 7, omega_k, p)
        omega7_dim = basis_7.shape[0]

        t_after7 = time.time()
        print(f" Omega_7={omega7_dim} ({t_after7-t_after6:.1f}s)...", end='', flush=True)

        # Compute rank(d_7): Omega_7 -> Omega_6
        _, _, rk_d7 = compute_boundary_rank(
            diff_seqs, allowed, 7, 6, omega_k, p,
            basis_lo=basis_6, basis_hi=basis_7)

        im_d7 = rk_d7
        beta6_k = ker_d6[k] - im_d7

        elapsed_k = time.time() - t_k
        elapsed = time.time() - t0
        print(f" rk(d7)={rk_d7} im={im_d7} beta6={beta6_k} [{elapsed_k:.1f}s / {elapsed:.1f}s total]")

        results[k] = {
            'omega6': omega6_dim,
            'omega7': omega7_dim,
            'rk_d7': rk_d7,
            'im_d7': im_d7,
            'ker_d6': ker_d6[k],
            'beta6': beta6_k,
        }

    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print(f"\n{'k':>4} {'Omega6':>8} {'Omega7':>8} {'ker(d6)':>8} {'rk(d7)':>8} {'im(d7)':>8} {'beta_6':>8}")
    for k in range(p):
        r = results[k]
        print(f"  {k:2d} {r['omega6']:8d} {r['omega7']:8d} {r['ker_d6']:8d} {r['rk_d7']:8d} {r['im_d7']:8d} {r['beta6']:8d}")

    total_beta6 = sum(results[k]['beta6'] for k in range(p))
    print(f"\nTotal beta_6(T_11) = {total_beta6}")
    print(f"Prediction (HYP-433): beta_6 = 10")

    # Check per-eigenspace chi
    # Per-eigenspace partial chi through deg 6: 1-5+20-70+205-460+700 = 391
    for k in range(p):
        r = results[k]
        partial_chi = 1 - 5 + 20 - 70 + 205 - 460 + r['omega6']
        print(f"  k={k}: partial_chi(0-6) = {partial_chi}")
        break  # they're all equal

    # Omega_7 dims
    omega7s = [results[k]['omega7'] for k in range(p)]
    omega7_all_equal = all(d == omega7s[0] for d in omega7s)
    print(f"\nAll Omega_7^(k) equal? {omega7_all_equal} (dim={omega7s[0]})")

    # beta_6 per eigenspace
    beta6s = [results[k]['beta6'] for k in range(p)]
    print(f"\nPer-eigenspace beta_6: {beta6s}")
    print(f"Eigenspaces with beta_6=1: {[k for k in range(p) if results[k]['beta6']==1]}")

    # rk(d7) analysis
    rk_d7s = [results[k]['rk_d7'] for k in range(p)]
    rk_all_equal = all(r == rk_d7s[0] for r in rk_d7s)
    print(f"\nAll rk(d_7^(k)) equal? {rk_all_equal}")
    print(f"rk(d_7) per eigenspace: {rk_d7s}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
