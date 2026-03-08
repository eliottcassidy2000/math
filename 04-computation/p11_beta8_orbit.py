#!/usr/bin/env python3
"""
p11_beta8_orbit.py - P_11 beta_8 via orbit-based computation

Since n=11 prime, all orbits have size 11, per-eigenspace dim = |A_p|/11.
By Aut(P_11) transitivity + conjugation, all k=1..10 have same beta_8.
So compute k=0 and k=1 only.

Strategy: work in orbit coordinates. For each orbit rep, compute boundary
faces and classify as allowed (in A_{p-1} orbit) or junk.
Build junk matrix incrementally using JHJ approach for memory efficiency.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
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
    """Return canonical orbit representative = lex-min rotation."""
    best = tuple(path)
    for j in range(1, n):
        rotated = tuple((v + j) % n for v in path)
        if rotated < best:
            best = rotated
    return best

def find_orbit_shift(path, rep, n):
    """Find shift j such that sigma^j(rep) = path."""
    for j in range(n):
        shifted = tuple((v + j) % n for v in rep)
        if shifted == tuple(path):
            return j
    return None

def find_orbits_fast(paths, n):
    """Find orbits and return (reps, path_to_orbit, path_to_shift)."""
    orbit_idx = {}  # rep -> orbit index
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


def compute_junk_and_boundary_orbit(reps_p, paths_pm1_set, orbit_idx_pm1,
                                     shift_pm1, n, k):
    """Compute junk matrix J and boundary matrix B in orbit coordinates.

    For eigenvalue k: coefficient contribution from face at shift j is
    omega^{-k*j} * (-1)^i.

    Returns: (J_dict, B_matrix, junk_reps)
    J_dict: {junk_orbit_rep: column_vector (length |reps_p|)}
    B_matrix: dense matrix (|orbits_pm1| x |reps_p|), complex
    """
    omega = np.exp(2j * np.pi / n) if k != 0 else 1.0

    n_orbits_p = len(reps_p)
    n_orbits_pm1 = max(orbit_idx_pm1.values()) + 1 if orbit_idx_pm1 else 0

    B = np.zeros((n_orbits_pm1, n_orbits_p), dtype=complex)
    junk_data = {}  # junk_orbit_rep -> dict{col_idx: coeff}

    for r, rep in enumerate(reps_p):
        p_len = len(rep) - 1  # path dimension
        for i in range(p_len + 1):
            face = rep[:i] + rep[i+1:]
            if face in paths_pm1_set:
                # Allowed face: goes to orbit of A_{p-1}
                q = orbit_idx_pm1[face]
                j = shift_pm1[face]
                if k == 0:
                    coeff = (-1.0)**i
                else:
                    coeff = (-1.0)**i * omega**(-k * j)
                B[q, r] += coeff
            else:
                # Junk face: classify by orbit
                jrep = canon_orbit_rep(face, n)
                if jrep not in junk_data:
                    junk_data[jrep] = {}
                j = find_orbit_shift(face, jrep, n)
                if k == 0:
                    coeff = (-1.0)**i
                else:
                    coeff = (-1.0)**i * omega**(-k * j)
                junk_data[jrep][r] = junk_data[jrep].get(r, 0) + coeff

    return junk_data, B


def compute_omega_via_junk(junk_data, n_orbits_p, tol=1e-8):
    """Compute Omega dimension = n_orbits_p - rank(J).

    Uses JHJ approach: G = J^H J = sum of outer products.
    """
    if not junk_data:
        return n_orbits_p  # no junk -> everything is Omega

    # Build G = J^H J incrementally
    G = np.zeros((n_orbits_p, n_orbits_p), dtype=complex)
    for jrep, col_dict in junk_data.items():
        # j_q is a sparse vector of length n_orbits_p
        j_q = np.zeros(n_orbits_p, dtype=complex)
        for col, val in col_dict.items():
            j_q[col] = val
        # G += j_q.conj() @ j_q.T (outer product)
        G += np.outer(j_q.conj(), j_q)

    # rank(J) = rank(G) since G = J^H J is PSD
    eigenvals = np.linalg.eigvalsh(G)
    rank = sum(abs(e) > tol for e in eigenvals)
    return n_orbits_p - rank


def compute_beta_p_orbit(B_p, omega_dim_p, omega_basis_p,
                          B_pp1, omega_dim_pp1, omega_basis_pp1,
                          tol=1e-6):
    """Compute beta_p in orbit coordinates.

    B_p: boundary d_p in orbit coords (orbits_{p-1} x orbits_p)
    omega_basis_p: basis for Omega_p in orbit coords (orbits_p x omega_dim_p)
    B_pp1: boundary d_{p+1} in orbit coords
    omega_basis_pp1: basis for Omega_{p+1}
    """
    if omega_dim_p == 0:
        return 0

    # d_p restricted to Omega_p
    bd_om = B_p @ omega_basis_p
    S_vals = np.linalg.svd(bd_om, compute_uv=False)
    rk_dp = sum(s > tol for s in S_vals)
    ker_dp = omega_dim_p - rk_dp

    # im(d_{p+1}) from Omega_{p+1}
    if omega_dim_pp1 == 0 or B_pp1 is None:
        im_dpp1 = 0
    else:
        bd_om_pp1 = B_p @ omega_basis_pp1  # Wait, B_pp1 goes to orbits_p, not orbits_{p-1}
        # Actually: d_{p+1} maps A_{p+1} -> A_p, and im(d_{p+1}) lands in ker(d_p)
        # In orbit coords: B_pp1: (orbits_p x orbits_{p+1})
        # Restricted to Omega_{p+1}: B_pp1 @ omega_basis_{p+1}
        bd_om_pp1 = B_pp1 @ omega_basis_pp1
        S_vals = np.linalg.svd(bd_om_pp1, compute_uv=False)
        im_dpp1 = sum(s > tol for s in S_vals)

    return max(0, int(round(ker_dp - im_dpp1)))


print("=" * 70)
print("P_11 BETA_8 (ORBIT-BASED)")
print("=" * 70)

n = 11
S = {1, 3, 4, 5, 9}
A = build_circulant(n, S)
t0 = time.time()

# Step 1: Enumerate allowed paths
print("\n--- Allowed paths ---")
allowed = {}
for p in range(10):
    t1 = time.time()
    allowed[p] = enumerate_allowed_paths(A, n, p)
    print(f"  A_{p}: {len(allowed[p])} ({time.time()-t1:.1f}s)")
    if not allowed[p]:
        break

# Step 2: Build orbit data for dims 6-9
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

# Step 3: Compute for k=0 (trivial eigenspace)
print("\n" + "=" * 70)
print("k=0 (TRIVIAL EIGENSPACE)")
print("=" * 70)

for p_target in [7, 8, 9]:
    if p_target not in orbit_reps:
        continue
    t1 = time.time()
    pm1 = p_target - 1
    junk_data, B = compute_junk_and_boundary_orbit(
        orbit_reps[p_target],
        path_sets.get(pm1, set()),
        orbit_idx.get(pm1, {}),
        orbit_shift.get(pm1, {}),
        n, k=0
    )
    n_orbits = len(orbit_reps[p_target])
    n_junk_orbits = len(junk_data)
    print(f"  d_{p_target}: B shape {B.shape}, junk orbits: {n_junk_orbits} ({time.time()-t1:.1f}s)")

    # Compute Omega dimension
    t1 = time.time()
    om_dim = compute_omega_via_junk(junk_data, n_orbits)
    print(f"  Om_{p_target}^(0): dim={om_dim} ({time.time()-t1:.1f}s)")

# Actually compute beta_8 for k=0 properly
print("\n  Computing beta_8^(0)...")

# Need Om_7, Om_8, Om_9 and boundaries d_8, d_9

# Compute Om and boundaries for p=7,8,9
omega_data = {}
boundary_data = {}

for p in [7, 8, 9]:
    pm1 = p - 1
    junk_data, B = compute_junk_and_boundary_orbit(
        orbit_reps[p], path_sets.get(pm1, set()),
        orbit_idx.get(pm1, {}), orbit_shift.get(pm1, {}),
        n, k=0
    )
    n_orb = len(orbit_reps[p])
    om_dim = compute_omega_via_junk(junk_data, n_orb)

    # Get Omega basis
    if not junk_data:
        om_basis = np.eye(n_orb, dtype=complex)
    else:
        G = np.zeros((n_orb, n_orb), dtype=complex)
        for jrep, col_dict in junk_data.items():
            j_q = np.zeros(n_orb, dtype=complex)
            for col, val in col_dict.items():
                j_q[col] = val
            G += np.outer(j_q.conj(), j_q)
        eigenvals, eigvecs = np.linalg.eigh(G)
        # Kernel = eigenvectors with eigenvalue ~0
        mask = np.abs(eigenvals) < 1e-8
        om_basis = eigvecs[:, mask]

    omega_data[p] = (om_dim, om_basis)
    boundary_data[p] = B
    print(f"  Om_{p}^(0) = {om_dim}")

# Compute beta_8^(0)
om8_dim, om8_basis = omega_data[8]
om9_dim, om9_basis = omega_data[9]
B_8 = boundary_data[8]  # orbits_7 x orbits_8
B_9 = boundary_data[9]  # orbits_8 x orbits_9

if om8_dim == 0:
    beta8_0 = 0
else:
    # ker(d_8) on Om_8
    bd8_om = B_8 @ om8_basis
    S_vals = np.linalg.svd(bd8_om, compute_uv=False)
    rk8 = sum(s > 1e-6 for s in S_vals)
    ker8 = om8_dim - rk8

    # im(d_9) from Om_9
    if om9_dim == 0:
        im9 = 0
    else:
        bd9_om = B_9 @ om9_basis
        S_vals = np.linalg.svd(bd9_om, compute_uv=False)
        im9 = sum(s > 1e-6 for s in S_vals)

    beta8_0 = max(0, int(round(ker8 - im9)))

print(f"\n  beta_8^(0) = {beta8_0}")
print(f"  (Expected: 0 if HYP-212 correct for prime n)")

# Step 4: Compute for k=1 (non-trivial eigenspace)
print("\n" + "=" * 70)
print("k=1 (NON-TRIVIAL EIGENSPACE)")
print("=" * 70)

omega_data_k1 = {}
boundary_data_k1 = {}

for p in [7, 8, 9]:
    pm1 = p - 1
    t1 = time.time()
    junk_data, B = compute_junk_and_boundary_orbit(
        orbit_reps[p], path_sets.get(pm1, set()),
        orbit_idx.get(pm1, {}), orbit_shift.get(pm1, {}),
        n, k=1
    )
    n_orb = len(orbit_reps[p])
    print(f"  d_{p}^(1): junk orbits={len(junk_data)} ({time.time()-t1:.1f}s)")

    t1 = time.time()
    om_dim = compute_omega_via_junk(junk_data, n_orb)

    # Get Omega basis
    if not junk_data:
        om_basis = np.eye(n_orb, dtype=complex)
    else:
        G = np.zeros((n_orb, n_orb), dtype=complex)
        for jrep, col_dict in junk_data.items():
            j_q = np.zeros(n_orb, dtype=complex)
            for col, val in col_dict.items():
                j_q[col] = val
            G += np.outer(j_q.conj(), j_q)
        eigenvals, eigvecs = np.linalg.eigh(G)
        mask = np.abs(eigenvals) < 1e-8
        om_basis = eigvecs[:, mask]

    omega_data_k1[p] = (om_dim, om_basis)
    boundary_data_k1[p] = B
    print(f"  Om_{p}^(1) = {om_dim} ({time.time()-t1:.1f}s)")

# Compute beta_8^(1)
om8_dim_k1, om8_basis_k1 = omega_data_k1[8]
om9_dim_k1, om9_basis_k1 = omega_data_k1[9]
B_8_k1 = boundary_data_k1[8]
B_9_k1 = boundary_data_k1[9]

if om8_dim_k1 == 0:
    beta8_1 = 0
else:
    bd8_om = B_8_k1 @ om8_basis_k1
    S_vals = np.linalg.svd(bd8_om, compute_uv=False)
    rk8 = sum(s > 1e-6 for s in S_vals)
    ker8 = om8_dim_k1 - rk8

    if om9_dim_k1 == 0:
        im9 = 0
    else:
        bd9_om = B_9_k1 @ om9_basis_k1
        S_vals = np.linalg.svd(bd9_om, compute_uv=False)
        im9 = sum(s > 1e-6 for s in S_vals)

    beta8_1 = max(0, int(round(ker8 - im9)))

print(f"\n  beta_8^(1) = {beta8_1}")

# Summary
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"  beta_8^(0) = {beta8_0} (trivial eigenspace)")
print(f"  beta_8^(1) = {beta8_1} (non-trivial, same for k=1..10)")
print(f"  TOTAL beta_8 = {beta8_0} + 10 * {beta8_1} = {beta8_0 + 10 * beta8_1}")
print(f"\n  HYP-212 prediction: beta_8 = 10 (= p-1 + 0 for prime p)")
print(f"  Matches: {beta8_0 + 10 * beta8_1 == 10}")

print(f"\nTotal time: {time.time()-t0:.1f}s")
