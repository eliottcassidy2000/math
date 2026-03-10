"""
t11_omega_dims_highdeg.py - Compute Omega_m^(eig) for T_11 at high degrees (7-10)

Uses a smaller prime (89, which has 11|(89-1)) to avoid memory issues.
The constraint matrix at deg 8 is (35145 x 14395) which needs 506MB in uint8
vs 3.77 GB in int64.

Checks the chi constraint: for chi(T_11)=11 per eigenspace we need
Omega_8 - Omega_9 + Omega_10 = 300

Author: kind-pasteur-2026-03-10-S50
"""
import sys
import time
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

# Use small prime with 11|(PRIME-1)
# 89: 89-1=88=8*11 => has 11th root of unity
# Verify: 89 is prime ✓
PRIME_SMALL = 89

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


def gauss_rank_mod_small(C_in, prime):
    """
    Gaussian elimination mod prime using int16 (for small prime < 256).
    C_in: 2D array of uint8 values in [0, prime-1]
    Returns rank.
    """
    C = C_in.astype(np.int16)
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
        # Scale pivot row
        C[pivot_row] = (C[pivot_row].astype(np.int32) * inv % prime).astype(np.int16)

        # Eliminate
        factor_col = C[:, col].copy()
        factor_col[pivot_row] = 0
        nonzero_rows = np.where(factor_col != 0)[0]
        for row in nonzero_rows:
            f = int(factor_col[row])
            C[row] = ((C[row].astype(np.int32) - f * C[pivot_row].astype(np.int32)) % prime).astype(np.int16)

        rank += 1
        pivot_row += 1

    return rank


def gauss_rank_mod_chunked(n_junk, n_Am, face_data_m, A_m, junk_idx, omega_k, prime, n):
    """
    Build constraint matrix in chunks and compute rank.
    Uses uint8 for storage (assumes prime < 256).
    """
    # Build full constraint matrix as uint8
    print(f"    Building {n_junk}x{n_Am} matrix (small prime {prime})...", end='', flush=True)
    t0 = time.time()

    C = np.zeros((n_junk, n_Am), dtype=np.uint8)
    for j, (D, faces) in enumerate(zip(A_m, face_data_m)):
        for fd, sign, offset, is_allowed in faces:
            if not is_allowed:
                row = junk_idx[fd]
                if offset != 0:
                    w = pow(omega_k, offset, prime)
                else:
                    w = 1
                entry = (sign * w) % prime
                if entry < 0:
                    entry = (entry + prime) % prime
                C[row, j] = (int(C[row, j]) + entry) % prime

    print(f" {time.time()-t0:.1f}s, running Gauss...", end='', flush=True)
    rk = gauss_rank_mod_small(C, prime)
    print(f" rank={rk}, dim={n_Am-rk} ({time.time()-t0:.1f}s)")
    return rk


def main():
    print("=" * 70)
    print("T_11 HIGH-DEGREE OMEGA DIMENSIONS (using small prime 89)")
    print("=" * 70)

    p = 11
    S = qr_set(p)
    max_deg = p - 1  # = 10

    prime = PRIME_SMALL
    print(f"Using PRIME = {prime} (11 | {prime-1})")

    # Find p-th root of unity mod prime
    omega_p = find_pth_root_of_unity(p, prime)
    print(f"omega_{p} mod {prime} = {omega_p}")
    assert omega_p is not None, f"No {p}-th root of unity mod {prime}"
    assert pow(omega_p, p, prime) == 1
    print(f"QR_{p} = {sorted(S)}")
    print()

    t_pre = time.time()
    diff_seqs = enumerate_diff_seqs(S, p, max_deg + 1)
    print(f"Diff seqs enumerated in {time.time()-t_pre:.1f}s")
    for m in range(max_deg + 2):
        print(f"  deg {m}: {len(diff_seqs.get(m, []))} seqs")
    print()

    # Build allowed sets
    allowed = {}
    for m in range(max_deg + 2):
        allowed[m] = set(diff_seqs.get(m, []))

    # Precompute junk and face data for high degrees
    print("Precomputing face structure for degrees 7-10...")
    junk_seqs = {}
    face_data = {}
    for m in range(8, max_deg + 1):  # degrees 8, 9, 10
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
        print(f"  deg {m}: |A|={len(A_m)}, |junk|={len(junk_seqs[m])}, "
              f"matrix size = {len(junk_seqs[m])}x{len(A_m)}")
    print()

    # Known dims from previous computation (degrees 0-7)
    known_dims = [1, 5, 20, 70, 205, 460, 700, 690]
    print(f"Known dims (deg 0-7): {known_dims}")
    partial_chi_7 = sum((-1)**m * known_dims[m] for m in range(8))
    print(f"Partial chi (deg 0-7) per eigenspace: {partial_chi_7}")
    print()

    # Compute dims for k=0 eigenspace (representative, all should be equal)
    # We'll check k=0 and k=1 for consistency
    dims_full = {}

    for k in [0, 1]:
        print(f"\n--- Eigenspace k={k} ---")
        omega_k = pow(omega_p, k, prime)
        dims_k = list(known_dims)  # copy deg 0-7

        for m in [8, 9, 10]:
            A_m = diff_seqs.get(m, [])
            junk_list = junk_seqs.get(m, [])
            n_junk = len(junk_list)
            n_Am = len(A_m)

            if n_Am == 0:
                dims_k.append(0)
                print(f"  deg {m}: Omega=0 (no seqs)")
                continue

            if n_junk == 0:
                dims_k.append(n_Am)
                print(f"  deg {m}: Omega={n_Am} (no junk)")
                continue

            t0 = time.time()
            print(f"  deg {m}: ", end='', flush=True)

            junk_idx = {j: i for i, j in enumerate(junk_list)}
            rk = gauss_rank_mod_chunked(
                n_junk, n_Am, face_data[m], A_m, junk_idx,
                omega_k, prime, p)
            dims_k.append(n_Am - rk)
            print(f"  deg {m}: Omega_{m}={n_Am-rk} [{time.time()-t0:.1f}s]")

        dims_full[k] = dims_k
        chi_k = sum((-1)**m * d for m, d in enumerate(dims_k))
        print(f"  Full dims (k={k}): {dims_k}")
        print(f"  Chi (k={k}): {chi_k}")

    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)

    k0 = dims_full[0]
    k1 = dims_full.get(1)

    print(f"\nPer-eigenspace Omega dims (k=0): {k0}")
    if k1:
        all_equal = (k0 == k1)
        print(f"Equal to k=1? {all_equal}")

    # Chi
    chi = sum((-1)**m * k0[m] for m in range(len(k0)))
    print(f"Chi per eigenspace: {chi}")
    print(f"Implied chi(T_11) = 11 * {chi} = {11*chi}")

    # Palindrome check Omega_k = Omega_{10-k}
    print(f"\nPalindrome check Omega_k = Omega_{{10-k}}:")
    for k in range(6):
        mirror = 10 - k
        if mirror < len(k0):
            ok = "OK" if k0[k] == k0[mirror] else f"FAIL ({k0[k]} vs {k0[mirror]})"
            print(f"  Omega_{k}={k0[k]} vs Omega_{mirror}={k0[mirror]}: {ok}")

    # Chi constraint check
    hi_chi = sum((-1)**m * k0[m] for m in range(8, len(k0)))
    print(f"\nChi from deg 8-10: {hi_chi}")
    print(f"Need {300} for chi=1 per eigenspace (if chi(T_11)=11)")


if __name__ == '__main__':
    main()
    print("\nDONE.")
