#!/usr/bin/env python3
"""
p11_beta8_v4.py - P_11 beta_8 via orbit-based computation with SPARSE junk matrix

Key optimization: build J as a SPARSE matrix instead of forming G = J^H J
via 15230 dense outer products. Use scipy.sparse.linalg.svds for rank.

Since n=11 is prime, all orbits have size 11.
By Aut(P_11) transitivity + conjugation, k=1..10 all give same beta_8.
So compute k=0 and k=1 only.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import svds
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import enumerate_allowed_paths
sys.stdout = _saved

def build_circulant(n, S):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            A[i][(i + s) % n] = 1
    return A

def canon_orbit_rep(path, n):
    best = tuple(path)
    for j in range(1, n):
        rotated = tuple((v + j) % n for v in path)
        if rotated < best:
            best = rotated
    return best

def find_orbit_shift(path, rep, n):
    for j in range(n):
        shifted = tuple((v + j) % n for v in rep)
        if shifted == tuple(path):
            return j
    return None

def find_orbits_fast(paths, n):
    orbit_idx = {}
    reps = []
    path_to_orbit = {}
    path_to_shift = {}
    for p in paths:
        pt = tuple(p)
        rep = canon_orbit_rep(p, n)
        if rep not in orbit_idx:
            orbit_idx[rep] = len(reps)
            reps.append(rep)
        idx = orbit_idx[rep]
        shift = find_orbit_shift(p, rep, n)
        path_to_orbit[pt] = idx
        path_to_shift[pt] = shift
    return reps, path_to_orbit, path_to_shift


def compute_boundary_and_junk_sparse(reps_p, paths_pm1_set, orbit_idx_pm1,
                                      shift_pm1, n, k):
    """Build B (boundary in orbit coords) and J (junk matrix) as SPARSE matrices.

    Returns: (B_sparse, J_sparse, n_junk_orbits)
    B: (n_orbits_pm1 x n_orbits_p) sparse
    J: (n_junk_orbits x n_orbits_p) sparse
    """
    omega = np.exp(2j * np.pi / n) if k != 0 else 1.0
    n_orbits_p = len(reps_p)
    n_orbits_pm1 = max(orbit_idx_pm1.values()) + 1 if orbit_idx_pm1 else 0

    # Collect boundary entries
    B_rows, B_cols, B_vals = [], [], []
    # Collect junk entries: first pass to identify junk orbits
    junk_orbits = {}  # jrep -> junk_orbit_index
    J_rows, J_cols, J_vals = [], [], []

    for r, rep in enumerate(reps_p):
        p_len = len(rep) - 1
        for i in range(p_len + 1):
            face = rep[:i] + rep[i+1:]
            if face in paths_pm1_set:
                q = orbit_idx_pm1[face]
                j = shift_pm1[face]
                if k == 0:
                    coeff = (-1.0)**i
                else:
                    coeff = (-1.0)**i * omega**(-k * j)
                B_rows.append(q)
                B_cols.append(r)
                B_vals.append(coeff)
            else:
                jrep = canon_orbit_rep(face, n)
                if jrep not in junk_orbits:
                    junk_orbits[jrep] = len(junk_orbits)
                jidx = junk_orbits[jrep]
                j = find_orbit_shift(face, jrep, n)
                if k == 0:
                    coeff = (-1.0)**i
                else:
                    coeff = (-1.0)**i * omega**(-k * j)
                J_rows.append(jidx)
                J_cols.append(r)
                J_vals.append(coeff)

    n_junk = len(junk_orbits)

    if k == 0:
        dtype = float
    else:
        dtype = complex

    B = sparse.csr_matrix((B_vals, (B_rows, B_cols)),
                           shape=(n_orbits_pm1, n_orbits_p), dtype=dtype)

    if n_junk > 0:
        # Aggregate duplicate entries (same junk orbit, same column)
        J = sparse.csr_matrix((J_vals, (J_rows, J_cols)),
                               shape=(n_junk, n_orbits_p), dtype=dtype)
    else:
        J = sparse.csr_matrix((n_junk, n_orbits_p), dtype=dtype)

    return B, J, n_junk


def sparse_rank(M, tol=1e-8, max_rank=None):
    """Compute rank of sparse matrix M efficiently.

    Uses iterative SVD starting from largest singular values.
    """
    m, n_cols = M.shape
    if m == 0 or n_cols == 0:
        return 0

    # For small matrices, just use dense SVD
    min_dim = min(m, n_cols)
    if min_dim <= 2000:
        M_dense = M.toarray() if sparse.issparse(M) else M
        S = np.linalg.svd(M_dense, compute_uv=False)
        return int(np.sum(np.abs(S) > tol))

    # For larger matrices, use dense on M^H M which is smaller
    if n_cols <= 5000:
        MHM = (M.conj().T @ M).toarray() if sparse.issparse(M) else M.conj().T @ M
        eigenvals = np.linalg.eigvalsh(MHM)
        return int(np.sum(np.abs(eigenvals) > tol**2))

    # Very large: try sparse SVD iteratively
    if max_rank is None:
        max_rank = min_dim

    # Binary search for rank using sparse svds
    try:
        k_try = min(min_dim - 1, max_rank)
        if k_try <= 0:
            return 0
        U, S, Vh = svds(M.astype(complex), k=k_try)
        return int(np.sum(S > tol))
    except Exception:
        # Fallback to dense
        M_dense = M.toarray() if sparse.issparse(M) else M
        S = np.linalg.svd(M_dense, compute_uv=False)
        return int(np.sum(np.abs(S) > tol))


def compute_beta_orbit_sparse(B_p, J_p, B_pp1, J_pp1, n_orbits_p, n_orbits_pm1, n_orbits_pp1):
    """Compute beta_p using sparse boundary and junk matrices.

    Omega_p = ker(J_p) (paths whose boundary lands in allowed space).
    beta_p = dim(ker(d_p on Om_p)) - dim(im(d_{p+1} from Om_{p+1}))

    d_p restricted to Om_p in orbit coords: B_p restricted to ker(J_p).
    """
    if n_orbits_p == 0:
        return 0, 0

    # Step 1: Compute Omega_p = ker(J_p)
    t1 = time.time()
    if J_p.nnz == 0:
        om_dim = n_orbits_p
        om_basis = np.eye(n_orbits_p, dtype=complex)
    else:
        # J_p is (n_junk x n_orbits_p)
        # We need kernel of J_p
        rk_J = sparse_rank(J_p)
        om_dim = n_orbits_p - rk_J
        if om_dim == 0:
            return 0, 0
        # Get kernel basis via dense SVD of J_p
        J_dense = J_p.toarray()
        U, S, Vh = np.linalg.svd(J_dense, full_matrices=True)
        om_basis = Vh[rk_J:].conj().T  # (n_orbits_p, om_dim)
    print(f"      Om dim={om_dim} (rank J={n_orbits_p - om_dim}, {time.time()-t1:.1f}s)")

    # Step 2: d_p restricted to Om_p: B_p @ om_basis
    t1 = time.time()
    if n_orbits_pm1 == 0:
        rk_dp = 0
    else:
        bd_om = B_p.toarray() @ om_basis if sparse.issparse(B_p) else B_p @ om_basis
        S_vals = np.linalg.svd(bd_om, compute_uv=False)
        rk_dp = int(np.sum(np.abs(S_vals) > 1e-6))
    ker_dp = om_dim - rk_dp
    print(f"      ker(d_p on Om)={ker_dp} (rk d_p={rk_dp}, {time.time()-t1:.1f}s)")

    # Step 3: Omega_{p+1} and im(d_{p+1})
    t1 = time.time()
    if n_orbits_pp1 == 0 or B_pp1 is None:
        im_dpp1 = 0
    else:
        if J_pp1.nnz == 0:
            om_dim_pp1 = n_orbits_pp1
            om_basis_pp1 = np.eye(n_orbits_pp1, dtype=complex)
        else:
            rk_Jpp1 = sparse_rank(J_pp1)
            om_dim_pp1 = n_orbits_pp1 - rk_Jpp1
            if om_dim_pp1 == 0:
                im_dpp1 = 0
                print(f"      Om_{{p+1}} dim=0, im=0 ({time.time()-t1:.1f}s)")
                return om_dim, max(0, ker_dp)
            J_dense_pp1 = J_pp1.toarray()
            U_pp1, S_pp1, Vh_pp1 = np.linalg.svd(J_dense_pp1, full_matrices=True)
            om_basis_pp1 = Vh_pp1[rk_Jpp1:].conj().T

        # im(d_{p+1}) in Om_p: project B_{p+1} @ om_basis_{p+1} onto Om_p
        bd_pp1_om = B_pp1.toarray() @ om_basis_pp1 if sparse.issparse(B_pp1) else B_pp1 @ om_basis_pp1
        S_vals = np.linalg.svd(bd_pp1_om, compute_uv=False)
        im_dpp1 = int(np.sum(np.abs(S_vals) > 1e-6))

    print(f"      im(d_p+1)={im_dpp1} ({time.time()-t1:.1f}s)")

    beta = max(0, int(round(ker_dp - im_dpp1)))
    return om_dim, beta


print("=" * 70)
print("P_11 BETA_8 (v4: SPARSE JUNK MATRIX)")
print("=" * 70)

n = 11
S = {1, 3, 4, 5, 9}
A = build_circulant(n, S)
t0 = time.time()

# Step 1: Enumerate allowed paths
print("\n--- Allowed paths ---")
allowed = {}
for p in range(11):
    t1 = time.time()
    allowed[p] = enumerate_allowed_paths(A, n, p)
    print(f"  A_{p}: {len(allowed[p])} ({time.time()-t1:.1f}s)")
    if not allowed[p]:
        break

# Step 2: Build orbit data
print("\n--- Orbit computation ---")
orbit_reps = {}
orbit_idx = {}
orbit_shift = {}
path_sets = {}

for p in [6, 7, 8, 9]:
    if p not in allowed or not allowed[p]:
        continue
    t1 = time.time()
    reps, oidx, oshift = find_orbits_fast(allowed[p], n)
    orbit_reps[p] = reps
    orbit_idx[p] = oidx
    orbit_shift[p] = oshift
    path_sets[p] = set(tuple(pa) for pa in allowed[p])
    print(f"  dim {p}: {len(reps)} orbits ({time.time()-t1:.1f}s)")

# Step 3: Compute beta_8 for k=0 and k=1
for k in [0, 1]:
    print(f"\n{'='*70}")
    print(f"k={k} {'(TRIVIAL)' if k==0 else '(NON-TRIVIAL)'} EIGENSPACE")
    print(f"{'='*70}")

    # Build sparse boundary and junk for p=7,8,9
    B_data = {}
    J_data = {}
    n_junk = {}

    for p in [7, 8, 9]:
        if p not in orbit_reps:
            continue
        pm1 = p - 1
        t1 = time.time()
        B, J, nj = compute_boundary_and_junk_sparse(
            orbit_reps[p],
            path_sets.get(pm1, set()),
            orbit_idx.get(pm1, {}),
            orbit_shift.get(pm1, {}),
            n, k
        )
        B_data[p] = B
        J_data[p] = J
        n_junk[p] = nj
        print(f"  d_{p}: B={B.shape} nnz={B.nnz}, J={J.shape} nnz={J.nnz} ({time.time()-t1:.1f}s)")

    # Compute beta_8
    print(f"\n  Computing beta_8^({k})...")
    n_orb = {p: len(orbit_reps[p]) for p in [7, 8, 9] if p in orbit_reps}

    om_dim, beta8_k = compute_beta_orbit_sparse(
        B_data.get(8), J_data.get(8),
        B_data.get(9), J_data.get(9),
        n_orb.get(8, 0), n_orb.get(7, 0), n_orb.get(9, 0)
    )

    if k == 0:
        beta8_0 = beta8_k
    else:
        beta8_1 = beta8_k

    print(f"\n  beta_8^({k}) = {beta8_k}")

# Summary
print(f"\n{'='*70}")
print("SUMMARY")
print("=" * 70)
print(f"  beta_8^(0) = {beta8_0} (trivial eigenspace)")
print(f"  beta_8^(1) = {beta8_1} (non-trivial, same for k=1..10)")
print(f"  TOTAL beta_8 = {beta8_0} + 10 * {beta8_1} = {beta8_0 + 10 * beta8_1}")
print(f"\n  HYP-212 prediction: beta_8 = 10 (= p-1 + 0 for prime p)")
print(f"  Matches: {beta8_0 + 10 * beta8_1 == 10}")

print(f"\nTotal time: {time.time()-t0:.1f}s")
