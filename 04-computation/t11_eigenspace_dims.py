"""
t11_eigenspace_dims.py - Per-eigenspace Omega dimensions for T_11 via mod-p arithmetic

Goal: Verify whether all 11 eigenspaces of T_11 have IDENTICAL Omega dims
(as they do for T_7). If true, this proves the Omega palindrome via equal eigenspace structure.

Key question: For T_11, does Omega_m^(k) = Omega_m / 11 for ALL k=0,...,10?

Method: mod-p Gaussian elimination, stopping at max_deg to avoid memory issues.

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
    """Find primitive p-th root of unity mod prime."""
    if (prime - 1) % p != 0:
        return None
    exp = (prime - 1) // p
    for g in range(2, 30):
        omega = pow(g, exp, prime)
        if omega == 1:
            continue
        if all(pow(omega, k, prime) != 1 for k in range(1, p)):
            return omega
    return None


def enumerate_diff_seqs(S, n, max_deg):
    """Enumerate all valid difference sequences up to max_deg.

    A diff sequence (s_1,...,s_m) is valid iff the m+1 partial sums
    {0, s_1, s_1+s_2, ..., s_1+...+s_m} are all DISTINCT mod n.
    This corresponds to the path visiting n distinct vertices.
    """
    S_sorted = sorted(S)
    seqs = {0: [()]}
    # partial_sums[seq] = frozenset of partial sums (including 0)
    partial_sums = {(): frozenset([0])}
    partial_last = {(): 0}  # last partial sum

    for m in range(1, max_deg + 1):
        prev = seqs[m - 1]
        new = []
        new_ps = {}
        new_last = {}
        for seq in prev:
            ps = partial_sums[seq]
            last = partial_last[seq]
            for s in S_sorted:
                new_last_sum = (last + s) % n
                if new_last_sum not in ps:
                    new_seq = seq + (s,)
                    new.append(new_seq)
                    new_ps[new_seq] = ps | frozenset([new_last_sum])
                    new_last[new_seq] = new_last_sum
        seqs[m] = new
        partial_sums.update(new_ps)
        partial_last.update(new_last)

    return seqs


def compute_face_diff(D, face_idx, n):
    """
    Compute face of diff-seq D at position face_idx.
    Returns (face_diff_seq, vertex_offset).

    face_0: drop first vertex, diff = D[1:], offset = D[0]
    face_i (0 < i < m): merge D[i-1]+D[i] (mod n), offset = 0
    face_m: drop last vertex, diff = D[:m-1], offset = 0
    """
    m = len(D)
    if face_idx == 0:
        return D[1:], D[0]
    elif face_idx == m:
        return D[:m - 1], 0
    else:
        merged = (D[face_idx - 1] + D[face_idx]) % n
        face = D[:face_idx - 1] + (merged,) + D[face_idx + 1:]
        return face, 0


def compute_eigenspace_omega_dims(S, n, max_deg, omega_p):
    """
    For each eigenspace k=0,...,n-1, compute Omega_m^(k) dimensions.

    Uses mod-PRIME arithmetic.
    omega_p = primitive n-th root of unity mod PRIME.
    """
    S_set = set(S)
    diff_seqs = enumerate_diff_seqs(S, n, max_deg + 1)

    # Build allowed set and index
    allowed = {}
    for m in range(max_deg + 2):
        ds = diff_seqs[m]
        allowed[m] = set(ds)

    # Precompute face structure for each degree
    face_data = {}
    junk_seqs = {}
    for m in range(1, max_deg + 2):
        A_m = diff_seqs[m]
        junk = set()
        face_list = []
        for D in A_m:
            faces = []
            for fi in range(m + 1):
                fd, offset = compute_face_diff(D, fi, n)
                sign = 1 if fi % 2 == 0 else -1
                is_allowed = (fd in allowed[m - 1])
                faces.append((fd, sign, offset, is_allowed))
                if not is_allowed:
                    junk.add(fd)
            face_list.append(faces)
        face_data[m] = face_list
        junk_seqs[m] = sorted(junk)

    print(f"\nS={sorted(S)}, n={n}")
    for m in range(max_deg + 2):
        n_junk = len(junk_seqs.get(m, []))
        n_A = len(diff_seqs[m])
        print(f"  deg {m}: |A|={n_A}, |junk|={n_junk}")

    # For each eigenspace k, compute rank of constraint matrix
    omega_dims = {}  # k -> list of Omega dims per degree
    t0 = time.time()

    for k in range(n):
        omega_k = pow(omega_p, k, PRIME)
        dims_k = [1]  # degree 0: Omega_0^(k) = 1

        for m in range(1, max_deg + 2):
            A_m = diff_seqs[m]
            junk_list = junk_seqs[m]
            n_junk = len(junk_list)
            n_Am = len(A_m)

            if n_Am == 0:
                dims_k.append(0)
                continue

            if n_junk == 0:
                dims_k.append(n_Am)
                continue

            # Build constraint matrix C (n_junk x n_Am)
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

            # Compute rank via Gaussian elimination mod PRIME
            rk = gauss_rank_mod(C, PRIME)
            dims_k.append(n_Am - rk)

        omega_dims[k] = dims_k[:max_deg + 1]  # trim to max_deg

        elapsed = time.time() - t0
        print(f"  k={k}: Omega_dims={omega_dims[k][:7]}{' ...' if max_deg >= 7 else ''} ({elapsed:.1f}s)")

    return omega_dims


def gauss_rank_mod(C, prime):
    """Gaussian elimination to find rank mod prime."""
    C = C.copy() % prime
    rows, cols = C.shape
    rank = 0
    pivot_row = 0

    for col in range(cols):
        # Find pivot
        found = -1
        for row in range(pivot_row, rows):
            if C[row, col] != 0:
                found = row
                break
        if found == -1:
            continue

        # Swap
        C[[pivot_row, found]] = C[[found, pivot_row]]

        # Eliminate below
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
    print("T_11 PER-EIGENSPACE OMEGA DIMENSIONS")
    print("=" * 70)

    p = 11
    S = qr_set(p)
    max_deg = 5  # Go up to degree 5 first (feasible), degree 6 is ~155MB per eigenspace

    omega_p = find_pth_root_of_unity(p, PRIME)
    print(f"omega_{p} = {omega_p} (mod {PRIME})")
    print(f"QR_{p} = {sorted(S)}")

    dims = compute_eigenspace_omega_dims(S, p, max_deg, omega_p)

    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")

    # Check if all eigenspaces have equal dims
    k0_dims = dims[0]
    all_equal = all(dims[k][:max_deg + 1] == k0_dims[:max_deg + 1] for k in range(1, p))
    print(f"\nAll eigenspaces have equal Omega dims? {all_equal}")

    # Compare with total/p
    print(f"\nPer-eigenspace dims (k=0): {k0_dims}")
    print(f"Expected (total/11): [1, 5, 20, 70, 205, 460, ...]")

    # Chi calculation per eigenspace (up to computed degrees)
    chi = sum((-1) ** m * k0_dims[m] for m in range(len(k0_dims)))
    print(f"\nPartial chi (degrees 0-{max_deg}) per eigenspace: {chi}")
    print(f"Expected total chi = p = {p}")

    # Identify the predicted Betti degree
    pred_deg = (p + 1) // 2
    print(f"\nPredicted Betti degree: {pred_deg}")
    print(f"Omega_{pred_deg}^(k) = {k0_dims[pred_deg] if pred_deg < len(k0_dims) else 'not computed'}")

    # Palindrome check within eigenspace
    print(f"\nPalindrome check Omega_k^(eig) = Omega_{{p-k}}^(eig) for k=1..{min(max_deg, p-1)}:")
    for k in range(1, min(max_deg + 1, p)):
        pk = p - k
        if pk <= max_deg:
            ok = k0_dims[k] == k0_dims[pk]
            print(f"  Omega_{k}={k0_dims[k]} vs Omega_{pk}={k0_dims[pk]}: {'OK' if ok else 'FAIL'}")

    # T_7 for comparison
    print(f"\n{'='*70}")
    print("T_7 PER-EIGENSPACE OMEGA DIMENSIONS (reference)")
    print(f"{'='*70}")
    p7 = 7
    S7 = qr_set(p7)
    omega_7 = find_pth_root_of_unity(p7, PRIME)
    dims7 = compute_eigenspace_omega_dims(S7, p7, p7 - 1, omega_7)

    print(f"\nAll eigenspaces equal? {all(dims7[k] == dims7[0] for k in range(1, p7))}")
    print(f"Per-eigenspace dims: {dims7[0]}")
    chi7 = sum((-1) ** m * dims7[0][m] for m in range(len(dims7[0])))
    print(f"Per-eigenspace chi: {chi7}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
