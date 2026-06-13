"""
paley_fourier_betti.py — Compute Paley tournament Betti numbers via Fourier decomposition

For the Paley tournament T_p (p prime, p = 3 mod 4), the shift tau: v -> v+1
is an automorphism. The chain complex decomposes into p eigenspaces.

In each eigenspace k, the chain groups have dimension = (number of diff-seq orbits).
The boundary maps restrict to these eigenspaces. Total Betti number is the sum
across eigenspaces.

Key advantage: at degree m, the eigenspace chain group has dimension O(d^m / p)
instead of O(d^m), where d = (p-1)/2 = out-degree.

Algorithm:
1. Enumerate allowed m-paths, group by difference sequence -> A_m
2. For each eigenspace k:
   a. Build boundary matrix BD from A_m diffs to ALL (m-1)-diff-seqs (allowed + non-allowed)
   b. Omega_m = ker(rows of BD corresponding to non-allowed faces)
   c. Boundary on Omega: restrict BD to allowed rows, projected onto Omega basis
   d. Betti_k(m) = dim(ker(d_m on Omega_m)) - dim(im(d_{m+1} on Omega_{m+1}))

Author: kind-pasteur-2026-03-10-S50
"""
import sys
import time
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    enumerate_all_allowed,
    _gauss_rank_np, _gauss_nullbasis_modp,
    boundary_faces, RANK_PRIME, matmul_mod
)

PRIME = RANK_PRIME


def paley_adj(p):
    """Construct Paley tournament for prime p = 3 mod 4."""
    qr = set()
    for i in range(1, p):
        qr.add((i * i) % p)
    A = [[0] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                A[i][j] = 1
    return A


def find_primitive_root_of_unity(p, prime):
    """Find a primitive p-th root of unity mod prime."""
    if (prime - 1) % p != 0:
        return None
    exp = (prime - 1) // p
    for g in [2, 3, 5, 6, 7, 10, 11, 13]:
        omega = pow(g, exp, prime)
        if omega == 1:
            continue
        # Check it's a PRIMITIVE p-th root (order exactly p)
        is_primitive = True
        for k in range(1, p):
            if pow(omega, k, prime) == 1:
                is_primitive = False
                break
        if is_primitive:
            return omega
    return None


def compute_face_diff(D, face_idx, m, p):
    """Compute the face diff sequence and eigenspace weight offset.

    Returns (face_diff, offset) where offset is the shift in starting vertex.
    face_0: delete first vertex -> offset = D[0]
    face_j (0 < j < m): merge D[j-1] and D[j] -> offset = 0
    face_m: delete last vertex -> offset = 0
    """
    if face_idx == 0:
        return D[1:], D[0]
    elif face_idx == m:
        return D[:m-1], 0
    else:
        # Merge D[face_idx-1] and D[face_idx]
        new_d = (D[face_idx - 1] + D[face_idx]) % p
        face_diff = D[:face_idx-1] + (new_d,) + D[face_idx+1:]
        return face_diff, 0


def fourier_betti(p, max_deg=None):
    """Compute Betti numbers of T_p using Fourier eigenspace decomposition."""
    if p % 4 != 3:
        print(f"  p={p} not = 3 mod 4, skipping")
        return None

    if max_deg is None:
        max_deg = p - 1

    A = paley_adj(p)

    # Find p-th root of unity mod PRIME
    omega = find_primitive_root_of_unity(p, PRIME)
    if omega is None:
        print(f"  ERROR: cannot find primitive {p}-th root of unity mod {PRIME}")
        return None
    print(f"  omega = {omega}, omega^{p} = {pow(omega, p, PRIME)} (mod {PRIME})")

    # Enumerate all allowed paths and group by difference sequence
    ap = enumerate_all_allowed(A, p, max_deg)

    # Build diff-sequence index for each degree
    diff_to_idx = {}  # degree -> {diff_seq: index}
    omega_bases = {}  # degree -> list of diff seqs

    for m in range(max_deg + 1):
        paths = ap.get(m, [])
        diff_seqs = {}
        for path in paths:
            diffs = tuple((path[i+1] - path[i]) % p for i in range(len(path) - 1))
            if diffs not in diff_seqs:
                diff_seqs[diffs] = len(diff_seqs)
        diff_to_idx[m] = diff_seqs
        omega_bases[m] = list(diff_seqs.keys())

    for m in range(max_deg + 1):
        print(f"  deg {m}: {len(omega_bases[m])} diff-seqs, {len(ap.get(m, []))} paths")

    # For each eigenspace k, compute Betti contribution
    total_betti = [0] * (max_deg + 1)

    for k in range(p):
        omega_k = pow(omega, k, PRIME)

        # For each degree m >= 1, build the full boundary matrix
        # from A_m diff-seqs to ALL (m-1)-diff-seqs that appear as faces
        # Then split into allowed (in A_{m-1}) and non-allowed rows

        # Step 1: Compute all face diff-seqs and build extended boundary matrix
        bd_matrices_allowed = {}  # m -> boundary matrix (A_{m-1} rows only)
        omega_basis_matrices = {}  # m -> Omega_m basis as matrix

        for m in range(max_deg + 1):
            diffs_m = omega_bases[m]
            if not diffs_m:
                omega_basis_matrices[m] = None
                continue

            if m == 0:
                # Degree 0: single empty diff-seq, always in Omega
                omega_basis_matrices[m] = np.eye(1, dtype=np.int64)
                bd_matrices_allowed[m] = None  # no boundary from degree 0
                continue

            # Collect all face diff-seqs (allowed and non-allowed)
            all_face_diffs = set()
            allowed_face_set = set(omega_bases[m - 1])

            for D in diffs_m:
                for face_idx in range(m + 1):
                    fd, _ = compute_face_diff(D, face_idx, m, p)
                    all_face_diffs.add(fd)

            # Separate into allowed and non-allowed
            non_allowed_faces = sorted([f for f in all_face_diffs if f not in allowed_face_set])
            allowed_faces = omega_bases[m - 1]  # keep original ordering

            na_idx = {f: i for i, f in enumerate(non_allowed_faces)}
            al_idx = {f: i for i, f in enumerate(allowed_faces)}

            dim_m = len(diffs_m)

            # Build constraint matrix (non-allowed face rows)
            if non_allowed_faces:
                C = np.zeros((len(non_allowed_faces), dim_m), dtype=np.int64)
                for j, D in enumerate(diffs_m):
                    for face_idx in range(m + 1):
                        fd, offset = compute_face_diff(D, face_idx, m, p)
                        if fd in na_idx:
                            sign = 1 if face_idx % 2 == 0 else PRIME - 1
                            weight = pow(omega_k, offset, PRIME) if (k != 0 and offset != 0) else 1
                            C[na_idx[fd], j] = (C[na_idx[fd], j] + sign * weight) % PRIME

                # Omega_m = ker(C)
                _, nb = _gauss_nullbasis_modp(C, C.shape[0], C.shape[1], PRIME)
                if nb:
                    ob = np.array(nb, dtype=np.int64) % PRIME
                else:
                    ob = np.zeros((0, dim_m), dtype=np.int64)
            else:
                # No non-allowed faces -> all of A_m is Omega_m
                ob = np.eye(dim_m, dtype=np.int64)

            omega_basis_matrices[m] = ob

            # Build boundary matrix (allowed face rows only)
            if allowed_faces:
                BD_al = np.zeros((len(allowed_faces), dim_m), dtype=np.int64)
                for j, D in enumerate(diffs_m):
                    for face_idx in range(m + 1):
                        fd, offset = compute_face_diff(D, face_idx, m, p)
                        if fd in al_idx:
                            sign = 1 if face_idx % 2 == 0 else PRIME - 1
                            weight = pow(omega_k, offset, PRIME) if (k != 0 and offset != 0) else 1
                            BD_al[al_idx[fd], j] = (BD_al[al_idx[fd], j] + sign * weight) % PRIME
                bd_matrices_allowed[m] = BD_al
            else:
                bd_matrices_allowed[m] = np.zeros((0, dim_m), dtype=np.int64)

        # Step 2: Compute Betti numbers for eigenspace k
        betti_k = []
        for m in range(max_deg + 1):
            ob_m = omega_basis_matrices[m]
            dim_omega_m = 0 if ob_m is None else ob_m.shape[0]

            # ker(d_m) in Omega_m
            if m == 0 or ob_m is None or dim_omega_m == 0:
                dim_ker = dim_omega_m
            else:
                BD_al = bd_matrices_allowed.get(m)
                if BD_al is not None and BD_al.size > 0 and ob_m.shape[0] > 0:
                    # d_m on Omega_m: BD_al @ ob_m.T
                    d_on_omega = matmul_mod(BD_al, ob_m.T, PRIME)
                    rk_d = int(_gauss_rank_np(d_on_omega.copy(), PRIME))
                else:
                    rk_d = 0
                dim_ker = dim_omega_m - rk_d

            # im(d_{m+1}) in Omega_m
            if m + 1 > max_deg:
                dim_im = 0
            else:
                ob_mp1 = omega_basis_matrices.get(m + 1)
                dim_omega_mp1 = 0 if ob_mp1 is None else ob_mp1.shape[0]
                if dim_omega_mp1 == 0 or ob_mp1 is None:
                    dim_im = 0
                else:
                    BD_al_mp1 = bd_matrices_allowed.get(m + 1)
                    if BD_al_mp1 is not None and BD_al_mp1.size > 0:
                        d_mp1_on_omega = matmul_mod(BD_al_mp1, ob_mp1.T, PRIME)
                        dim_im = int(_gauss_rank_np(d_mp1_on_omega.copy(), PRIME))
                    else:
                        dim_im = 0

            b_m = dim_ker - dim_im
            betti_k.append(b_m)

        for m in range(max_deg + 1):
            total_betti[m] += betti_k[m]

        # Print nonzero eigenspace contributions
        if any(b > 0 for b in betti_k[1:]):
            nonzero = [(m, betti_k[m]) for m in range(max_deg + 1) if betti_k[m] > 0 and m > 0]
            print(f"  k={k}: betti={betti_k} nonzero={nonzero}")

    return {m: total_betti[m] for m in range(max_deg + 1)}


def main():
    print("=" * 70)
    print("PALEY TOURNAMENT FOURIER BETTI NUMBERS")
    print("=" * 70)

    for p in [3, 7, 11]:
        print(f"\n--- T_{p} ---")
        max_deg = min(p - 1, 8)
        if p >= 19:
            max_deg = 5

        t0 = time.time()
        bettis = fourier_betti(p, max_deg)
        elapsed = time.time() - t0

        if bettis:
            print(f"  RESULT: bettis = {dict(sorted(bettis.items()))}")
            print(f"  Euler char = {sum((-1)**k * v for k, v in bettis.items())}")
            print(f"  Time: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
    print("\nDONE.")
