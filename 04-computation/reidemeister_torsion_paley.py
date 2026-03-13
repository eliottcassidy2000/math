#!/usr/bin/env python3
"""
reidemeister_torsion_paley.py — opus-2026-03-13-S67i

COMPUTING: Reidemeister torsion of the GLMY path homology chain complex
for Paley tournaments P_7 and P_11.

The chain complex: 0 -> C_{p-1} -> ... -> C_1 -> C_0 -> 0
with boundary d_d : C_d -> C_{d-1}.

For an acyclic complex, torsion = alternating product of det'(d_d).
For a complex with homology concentrated at d = 0, m, m+1:
  tau = product over "acyclic chunks"

The key question: does tau relate to F_p (Fibonacci number)?

We use the eigenspace decomposition: each k=0,...,p-1 eigenspace
gives a sub-complex, and the total torsion is the product.

For k != 0: the sub-complex is acyclic except at d=m+1 where
beta_{m+1}^{(k)} = 1. The torsion of each k-eigenspace is computable.

For k = 0: more complex, with beta_m^{(0)} and beta_{m+1}^{(0)} nonzero.
"""

import numpy as np
import math

def build_paley(p):
    """Build Paley tournament adjacency matrix."""
    m = (p - 1) // 2
    QR = set()
    for k in range(1, p):
        QR.add((k * k) % p)

    A = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in QR:
                A[i][j] = 1
    return A

def allowed_paths(A, d):
    """Find all allowed d-paths in tournament A.
    An allowed d-path is a sequence (v_0, v_1, ..., v_d) where
    all v_i are distinct AND v_i -> v_{i+1} is an arc for all i,
    AND the path is 'regular' (i.e., in the GLMY sense:
    v_0 -> v_1 -> ... -> v_d with all arcs present).

    Actually for GLMY path homology, the allowed paths are
    all directed paths of length d (d+1 vertices, d arcs).
    """
    n = A.shape[0]
    if d == 0:
        return [[v] for v in range(n)]

    # Build recursively
    paths_prev = allowed_paths(A, d - 1)
    paths = []
    for path in paths_prev:
        last = path[-1]
        for v in range(n):
            if v not in path and A[last][v] == 1:
                paths.append(path + [v])
    return paths

def boundary_matrix(A, d, paths_d=None, paths_dm1=None):
    """Build the boundary matrix d_d : C_d -> C_{d-1}.

    The GLMY boundary: d_d(v_0,...,v_d) = sum_{i=0}^{d} (-1)^i (v_0,...,hat{v_i},...,v_d)
    where the face (v_0,...,hat{v_i},...,v_d) is included only if it's allowed.
    """
    n = A.shape[0]
    if paths_d is None:
        paths_d = allowed_paths(A, d)
    if paths_dm1 is None:
        paths_dm1 = allowed_paths(A, d - 1)

    # Index the paths
    path_to_idx_d = {tuple(p): i for i, p in enumerate(paths_d)}
    path_to_idx_dm1 = {tuple(p): i for i, p in enumerate(paths_dm1)}

    rows = len(paths_dm1)
    cols = len(paths_d)
    D = np.zeros((rows, cols), dtype=float)

    for j, path in enumerate(paths_d):
        for i_face in range(d + 1):
            face = tuple(path[:i_face] + path[i_face+1:])
            # Check if face is an allowed path
            if face in path_to_idx_dm1:
                # Check that all arcs in face are present
                valid = True
                for k in range(len(face) - 1):
                    if A[face[k]][face[k+1]] != 1:
                        valid = False
                        break
                if valid:
                    sign = (-1) ** i_face
                    row = path_to_idx_dm1[face]
                    D[row][j] += sign

    return D

print("=" * 70)
print("REIDEMEISTER TORSION OF PALEY TOURNAMENTS")
print("=" * 70)

for p in [5, 7]:
    print(f"\n{'='*60}")
    print(f"  P_{p} (m = {(p-1)//2})")
    print(f"{'='*60}")

    A = build_paley(p)
    m = (p - 1) // 2

    # Compute all path spaces and boundary matrices
    all_paths = {}
    for d in range(p):
        all_paths[d] = allowed_paths(A, d)
        if len(all_paths[d]) == 0:
            break
    max_d = max(d for d in all_paths if len(all_paths[d]) > 0)

    print(f"\n  Chain complex dimensions:")
    for d in sorted(all_paths.keys()):
        if len(all_paths[d]) > 0:
            print(f"    C_{d} = {len(all_paths[d])}")

    # Build boundary matrices and compute ranks
    boundaries = {}
    for d in range(1, max_d + 1):
        if d in all_paths and (d-1) in all_paths:
            if len(all_paths[d]) > 0 and len(all_paths[d-1]) > 0:
                D = boundary_matrix(A, d, all_paths[d], all_paths[d-1])
                boundaries[d] = D
                rk = np.linalg.matrix_rank(D)
                print(f"    d_{d}: {D.shape[0]}x{D.shape[1]}, rank={rk}")

    # Compute Betti numbers
    print(f"\n  Betti numbers:")
    for d in range(max_d + 1):
        dim_d = len(all_paths.get(d, []))
        im_dp1 = np.linalg.matrix_rank(boundaries[d+1]) if (d+1) in boundaries else 0
        im_d = np.linalg.matrix_rank(boundaries[d]) if d in boundaries else 0
        # ker(d_d) = dim_d - im_d (since d maps C_d -> C_{d-1})
        # For d=0: ker(d_0) is undefined (no map from C_0), but C_0 has no outgoing boundary
        # Actually: ker(d_d) = dim(C_d) - rank(d_d)
        # beta_d = ker(d_d) - im(d_{d+1}) = dim(C_d) - rank(d_d) - rank(d_{d+1})

        rank_out = np.linalg.matrix_rank(boundaries[d]) if d in boundaries else 0
        rank_in = np.linalg.matrix_rank(boundaries[d+1]) if (d+1) in boundaries else 0
        beta = dim_d - rank_out - rank_in
        if beta != 0 or dim_d > 0:
            print(f"    beta_{d} = {beta}  (dim={dim_d}, rank_out={rank_out}, rank_in={rank_in})")

    # Compute torsion
    # For a chain complex C_0 <- C_1 <- ... <- C_n:
    # Split each C_d = im(d_{d+1}) + complement
    # Torsion = alternating product of |det of d restricted to complement|

    # Actually, for the torsion we need to compute det'(d_d^T * d_d)
    # where det' is the product of nonzero eigenvalues.

    print(f"\n  Torsion computation (det' of d_d^T d_d):")
    log_torsion = 0
    for d in sorted(boundaries.keys()):
        D = boundaries[d]
        DDT = D @ D.T  # d_d * d_d^T
        DTD = D.T @ D  # d_d^T * d_d

        # Nonzero eigenvalues of DTD (same as DDT)
        eigs = np.linalg.eigvalsh(DTD)
        nonzero_eigs = [e for e in eigs if abs(e) > 1e-10]

        if nonzero_eigs:
            log_det_prime = sum(math.log(abs(e)) for e in nonzero_eigs)
            # det'(DTD) = prod of nonzero eigenvalues
            # det'(D) = sqrt(det'(DTD))
            log_det_D = log_det_prime / 2
        else:
            log_det_D = 0

        sign = (-1) ** (d + 1)  # alternating in the torsion formula
        log_torsion += sign * log_det_D

        print(f"    d_{d}: {len(nonzero_eigs)} nonzero eigs, "
              f"log(det')={log_det_D:.6f}, "
              f"sign=(-1)^{d+1}={'+'if sign>0 else '-'}")

    torsion = math.exp(log_torsion) if abs(log_torsion) < 100 else float('inf')

    # Fibonacci number
    fib = [0, 1]
    while fib[-1] < 10**15:
        fib.append(fib[-1] + fib[-2])
    F_p = fib[p]

    print(f"\n  TORSION: tau(P_{p}) = exp({log_torsion:.6f}) = {torsion:.6f}")
    print(f"  F_{p} = {F_p}")
    print(f"  tau / F_p = {torsion / F_p:.6f}" if F_p > 0 else "  F_p = 0")
    print(f"  log(tau) / log(F_p) = {log_torsion / math.log(F_p):.6f}" if F_p > 1 else "")
    print(f"  log(tau) = {log_torsion:.6f}")
    print(f"  log(F_p) = {math.log(F_p):.6f}" if F_p > 0 else "")

    # Also compare with det(I + A)
    det_IpA = abs(np.linalg.det(np.eye(p) + A.astype(float)))
    print(f"  det(I+A) = {det_IpA:.2f}")
    print(f"  log(det(I+A)) = {math.log(det_IpA):.6f}" if det_IpA > 0 else "")

print("\n\nDONE — reidemeister_torsion_paley.py complete")
