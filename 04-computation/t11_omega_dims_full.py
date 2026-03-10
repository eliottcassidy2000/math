"""
t11_omega_dims_full.py - Compute Omega_m^(eig) dimensions for ALL degrees 0-10 for T_11

This extends t11_eigenspace_dims.py to the full range.
For n=11 paths can visit at most 11 vertices, so max degree = 10.

Only computes DIMENSIONS (rank of constraint matrix), not full bases.
This is faster than computing actual null space bases.

Key constraint: chi per eigenspace = 1 requires
  1-5+20-70+205-460+700 - Omega_7 + Omega_8 - Omega_9 + Omega_10 = 1
  => -Omega_7 + Omega_8 - Omega_9 + Omega_10 = -390

Author: kind-pasteur-2026-03-10-S50
"""
import sys
import time
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import RANK_PRIME

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


def enumerate_diff_seqs(S, n, max_deg):
    """Enumerate all valid difference sequences up to max_deg."""
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


def main():
    print("=" * 70)
    print("T_11 FULL OMEGA DIMENSIONS (degrees 0-10, all eigenspaces)")
    print("=" * 70)

    p = 11
    S = qr_set(p)
    max_deg = p - 1  # = 10 (maximum possible for n=11)

    omega_p = find_pth_root_of_unity(p, PRIME)
    print(f"omega_{p} = {omega_p} (mod {PRIME})")
    print(f"QR_{p} = {sorted(S)}")
    print()

    t_pre = time.time()
    # Need deg+1 for face computation
    diff_seqs = enumerate_diff_seqs(S, p, max_deg + 1)
    print(f"Diff seqs enumerated in {time.time()-t_pre:.1f}s")
    for m in range(max_deg + 2):
        n = len(diff_seqs.get(m, []))
        print(f"  deg {m}: {n} seqs")
    print()

    # Precompute face structure and junk for each degree
    allowed = {}
    for m in range(max_deg + 2):
        allowed[m] = set(diff_seqs.get(m, []))

    junk_seqs = {}
    face_data = {}
    print("Precomputing face structure...")
    for m in range(1, max_deg + 2):
        A_m = diff_seqs.get(m, [])
        if not A_m:
            junk_seqs[m] = []
            face_data[m] = []
            continue
        junk = set()
        flist = []
        for D in A_m:
            faces = []
            for fi in range(m + 1):
                fd, offset = compute_face_diff(D, fi, p)
                sign = 1 if fi % 2 == 0 else -1
                is_allowed = (fd in allowed[m - 1])
                faces.append((fd, sign, offset, is_allowed))
                if not is_allowed:
                    junk.add(fd)
            flist.append(faces)
        junk_seqs[m] = sorted(junk)
        face_data[m] = flist
        print(f"  deg {m}: |A|={len(A_m)}, |junk|={len(junk_seqs[m])}")

    print()

    # For each eigenspace k, compute Omega_m dim for all m
    omega_dims = {}
    t0 = time.time()

    for k in range(p):
        omega_k = pow(omega_p, k, PRIME)
        dims_k = [1]  # degree 0

        for m in range(1, max_deg + 1):
            A_m = diff_seqs.get(m, [])
            junk_list = junk_seqs.get(m, [])
            n_junk = len(junk_list)
            n_Am = len(A_m)

            if n_Am == 0:
                dims_k.append(0)
                continue
            if n_junk == 0:
                dims_k.append(n_Am)
                continue

            junk_idx = {j: i for i, j in enumerate(junk_list)}
            C = np.zeros((n_junk, n_Am), dtype=np.int64)

            for j, (D, faces) in enumerate(zip(A_m, face_data[m])):
                for fd, sign, offset, is_allowed in faces:
                    if not is_allowed:
                        row = junk_idx[fd]
                        w = pow(omega_k, offset, PRIME) if offset != 0 else 1
                        entry = (sign * w) % PRIME
                        if entry < 0:
                            entry += PRIME
                        C[row, j] = (C[row, j] + entry) % PRIME

            rk = gauss_rank_mod(C, PRIME)
            dims_k.append(n_Am - rk)

        omega_dims[k] = dims_k
        elapsed = time.time() - t0
        chi_k = sum((-1)**m * dims_k[m] for m in range(len(dims_k)))
        print(f"  k={k}: dims={dims_k}, chi={chi_k} ({elapsed:.1f}s)")

    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)

    k0_dims = omega_dims[0]
    print(f"\nPer-eigenspace dims (k=0): {k0_dims}")
    print(f"Expected palindrome T_7: [1,3,6,9,9,6,3]")

    # Check all eigenspaces equal
    all_equal = all(omega_dims[k] == k0_dims for k in range(1, p))
    print(f"\nAll eigenspaces equal? {all_equal}")

    # Chi per eigenspace
    chi_per = sum((-1)**m * k0_dims[m] for m in range(len(k0_dims)))
    print(f"chi per eigenspace: {chi_per} (should be 1 if chi(T_11)=11)")
    print(f"Total chi(T_11): {p * chi_per} (should be 11)")

    # Palindrome check
    print(f"\nPalindrome check Omega_m = Omega_{{p-1-m}}:")
    for m in range(len(k0_dims) // 2 + 1):
        mirror = (p - 1) - m
        if 0 <= mirror < len(k0_dims):
            ok = "OK" if k0_dims[m] == k0_dims[mirror] else f"FAIL ({k0_dims[m]} != {k0_dims[mirror]})"
            print(f"  Omega_{m}={k0_dims[m]} vs Omega_{mirror}={k0_dims[mirror]}: {ok}")

    # chi constraint check
    alt_sum_7plus = sum((-1)**m * k0_dims[m] for m in range(7, len(k0_dims)))
    print(f"\nchi from deg>=7: {alt_sum_7plus}")
    print(f"chi from deg 0-6: {sum((-1)**m * k0_dims[m] for m in range(7))}")
    print(f"chi constraint: sum must be -390 for chi=1 per eigenspace? Actually need 1-391=-390")

    # Total Omega dims
    print(f"\nTotal Omega dims (=11*per-eigenspace):")
    for m, d in enumerate(k0_dims):
        print(f"  Omega_{m}(T_11) = {11*d}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
