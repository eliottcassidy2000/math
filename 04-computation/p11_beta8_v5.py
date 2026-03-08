#!/usr/bin/env python3
"""
p11_beta8_v5.py - P_11 beta_8 via RANDOMIZED rank estimation

The orbit-based approach still needs rank of J (35145 x 14395 for k=0).
Dense SVD of this is infeasible.

Strategy: use random projection to estimate rank.
- Random projection: Phi = random (d x m) matrix, compute rank(Phi @ J^T)
  where d ~ 2*expected_rank. If rank = expected rank, we're confident.
- For iterative nullspace: use randomized range finder to get approximate
  basis for range(J^T), then nullspace = complement.

Key insight: for k=0, all matrices are REAL (no complex phases).

Alternative: exploit the FULL automorphism group of P_11.
Aut(P_p) = {x -> ax+b : a in QR(p), b in Z_p}, order p(p-1)/2 = 55.
This gives orbits of size 55 (not just 11), reducing dimensions by 5x more.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from scipy import sparse
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

def apply_aut(path, a, b, n):
    """Apply affine map x -> a*x + b mod n."""
    return tuple((a * v + b) % n for v in path)

def canon_aut_rep(path, n, QR):
    """Return canonical orbit representative under full Aut(P_p)."""
    best = tuple(path)
    for a in QR:
        for b in range(n):
            img = apply_aut(path, a, b, n)
            if img < best:
                best = img
    return best

def find_aut_orbits(paths, n, QR):
    """Find orbits under full automorphism group.

    Returns: (reps, path_to_orbit_idx)
    """
    orbit_idx = {}  # rep -> index
    reps = []
    path_to_orbit = {}

    for p in paths:
        pt = tuple(p)
        rep = canon_aut_rep(p, n, QR)
        if rep not in orbit_idx:
            orbit_idx[rep] = len(reps)
            reps.append(rep)
        path_to_orbit[pt] = orbit_idx[rep]

    return reps, path_to_orbit

def find_aut_shift(path, rep, n, QR):
    """Find (a,b) such that x->ax+b maps rep to path."""
    for a in QR:
        for b in range(n):
            img = apply_aut(rep, a, b, n)
            if img == tuple(path):
                return (a, b)
    return None


print("=" * 70)
print("P_11 BETA_8 (v5: FULL AUTOMORPHISM GROUP)")
print("=" * 70)

n = 11
QR = [1, 3, 4, 5, 9]  # quadratic residues mod 11
S = set(QR)
aut_size = n * len(QR)  # = 55
A = build_circulant(n, S)
t0 = time.time()

print(f"\nAut(P_11) = AGL(1, F_11) = {{x -> ax+b : a in QR}}, |Aut| = {aut_size}")

# Step 1: Enumerate allowed paths
print("\n--- Allowed paths ---")
allowed = {}
for p in range(11):
    t1 = time.time()
    allowed[p] = enumerate_allowed_paths(A, n, p)
    print(f"  A_{p}: {len(allowed[p])} ({time.time()-t1:.1f}s)")
    if not allowed[p]:
        break

# Step 2: Find orbits under full automorphism group
print("\n--- Full automorphism orbits ---")
aut_reps = {}
aut_orb_map = {}
for p in [6, 7, 8, 9]:
    if p not in allowed or not allowed[p]:
        continue
    t1 = time.time()
    reps, orb_map = find_aut_orbits(allowed[p], n, QR)
    aut_reps[p] = reps
    aut_orb_map[p] = orb_map
    exp_orbits = len(allowed[p]) // aut_size
    print(f"  dim {p}: {len(reps)} orbits (expect ~{exp_orbits}, ratio={len(allowed[p])/len(reps):.1f}) ({time.time()-t1:.1f}s)")

# The orbit sizes may vary (not all orbits have full size 55)
# Check orbit size distribution
for p in [7, 8, 9]:
    if p not in aut_reps:
        continue
    orbit_sizes = {}
    for rep in aut_reps[p]:
        # Count orbit size
        orbit = set()
        for a in QR:
            for b in range(n):
                img = apply_aut(rep, a, b, n)
                if img in aut_orb_map[p]:
                    orbit.add(img)
        sz = len(orbit)
        orbit_sizes[sz] = orbit_sizes.get(sz, 0) + 1
    print(f"  dim {p} orbit sizes: {dict(sorted(orbit_sizes.items()))}")

# Step 3: For the trivial representation (k=0 for ALL elements of Aut),
# the invariant subspace consists of orbit-sum vectors.
# Boundary in orbit coords: for each orbit rep, compute boundary faces,
# classify as allowed (which orbit?) or junk (which junk orbit?).
# Coefficient = sum over orbit elements of boundary coefficient / |orbit|

# For k=0 under FULL Aut, all orbit contributions are uniform (coefficient 1).
# The boundary map d_p sends orbit_rep r to:
#   d(r) = sum_i (-1)^i face_i(r)
# Each face lands in some orbit (allowed or junk).
# In orbit coords, the contribution to orbit q from face_i(r) is:
#   (-1)^i if face_i(r) is in orbit q

print("\n--- Building orbit boundary matrices ---")

path_sets = {}
for p in [6, 7, 8, 9]:
    if p in allowed:
        path_sets[p] = set(tuple(pa) for pa in allowed[p])

# For trivial rep: boundary in orbit coordinates
def build_orbit_boundary_trivial(reps_p, path_set_pm1, orb_map_pm1, n, QR):
    """Build boundary matrix in orbit coordinates for trivial representation.

    For trivial rep, each orbit rep contributes equally (all aut images
    have same boundary structure up to relabeling). So we just need to
    compute the boundary of the rep and classify faces.

    Returns: (B, junk_reps, J)
    B: dense matrix (n_orbits_pm1 x n_orbits_p), real
    J: dense matrix (n_junk_orbits x n_orbits_p), real
    """
    n_orbits_p = len(reps_p)
    n_orbits_pm1 = max(orb_map_pm1.values()) + 1 if orb_map_pm1 else 0

    junk_orbit_map = {}  # junk_rep -> index
    junk_reps = []

    # First pass: identify all junk orbits
    for r, rep in enumerate(reps_p):
        p_len = len(rep) - 1
        for i in range(p_len + 1):
            face = rep[:i] + rep[i+1:]
            if face not in path_set_pm1:
                jrep = canon_aut_rep(face, n, QR)
                if jrep not in junk_orbit_map:
                    junk_orbit_map[jrep] = len(junk_reps)
                    junk_reps.append(jrep)

    n_junk = len(junk_reps)
    B = np.zeros((n_orbits_pm1, n_orbits_p), dtype=float)
    J = np.zeros((n_junk, n_orbits_p), dtype=float)

    for r, rep in enumerate(reps_p):
        p_len = len(rep) - 1
        for i in range(p_len + 1):
            face = rep[:i] + rep[i+1:]
            sign = (-1.0)**i
            if face in path_set_pm1:
                q = orb_map_pm1[face]
                B[q, r] += sign
            else:
                jrep = canon_aut_rep(face, n, QR)
                jidx = junk_orbit_map[jrep]
                J[jidx, r] += sign

    return B, J, n_junk

B_data = {}
J_data = {}
n_junk = {}

for p in [7, 8, 9]:
    if p not in aut_reps:
        continue
    pm1 = p - 1
    t1 = time.time()
    B, J, nj = build_orbit_boundary_trivial(
        aut_reps[p],
        path_sets.get(pm1, set()),
        aut_orb_map.get(pm1, {}),
        n, QR
    )
    B_data[p] = B
    J_data[p] = J
    n_junk[p] = nj
    print(f"  d_{p}: B=({B.shape[0]}x{B.shape[1]}), J=({J.shape[0]}x{J.shape[1]}) ({time.time()-t1:.1f}s)")

# Step 4: Compute beta_8 for trivial representation (k=0 under full Aut)
print(f"\n{'='*70}")
print("TRIVIAL REPRESENTATION (k=0, full Aut)")
print(f"{'='*70}")

for p in [7, 8, 9]:
    if p not in J_data:
        continue
    t1 = time.time()
    J = J_data[p]
    n_orb = len(aut_reps[p])
    if J.shape[0] == 0:
        om_dim = n_orb
    else:
        S_vals = np.linalg.svd(J, compute_uv=False)
        rk_J = int(np.sum(np.abs(S_vals) > 1e-8))
        om_dim = n_orb - rk_J
    print(f"  Om_{p}^(triv): dim={om_dim} (rk J={n_orb-om_dim}) ({time.time()-t1:.1f}s)")

# Now compute beta_8^(triv)
J8 = J_data[8]
n_orb8 = len(aut_reps[8])
S_vals = np.linalg.svd(J8, compute_uv=False)
rk_J8 = int(np.sum(np.abs(S_vals) > 1e-8))
om8_dim = n_orb8 - rk_J8

if om8_dim == 0:
    beta8_triv = 0
    print(f"\n  beta_8^(triv) = 0 (Om_8 = 0)")
else:
    # Get Om_8 basis via J^H J eigenvectors (kernel = eigenvalue ~0)
    JHJ8 = J8.T @ J8
    eigenvals8, eigvecs8 = np.linalg.eigh(JHJ8)
    mask8 = np.abs(eigenvals8) < 1e-10
    om8_basis = eigvecs8[:, mask8]
    print(f"  Om_8 basis shape: {om8_basis.shape}")

    # ker(d_8 on Om_8): d_8 restricted to Om_8
    bd8_om = B_data[8] @ om8_basis
    S_vals_bd = np.linalg.svd(bd8_om, compute_uv=False)
    rk_d8 = int(np.sum(np.abs(S_vals_bd) > 1e-6))
    ker_d8 = om8_dim - rk_d8

    # im(d_9 from Om_9)
    J9 = J_data[9]
    n_orb9 = len(aut_reps[9])
    JHJ9 = J9.T @ J9
    eigenvals9, eigvecs9 = np.linalg.eigh(JHJ9)
    rk_J9 = int(np.sum(np.abs(eigenvals9) > 1e-10))
    om9_dim = n_orb9 - rk_J9

    if om9_dim == 0:
        im_d9 = 0
    else:
        mask9 = np.abs(eigenvals9) < 1e-10
        om9_basis = eigvecs9[:, mask9]
        bd9_om = B_data[9] @ om9_basis
        S_vals_im = np.linalg.svd(bd9_om, compute_uv=False)
        im_d9 = int(np.sum(np.abs(S_vals_im) > 1e-6))

    beta8_triv = max(0, ker_d8 - im_d9)
    print(f"\n  Om_8^(triv) = {om8_dim}, ker(d_8) = {ker_d8}, im(d_9) = {im_d9}")
    print(f"  beta_8^(triv) = {beta8_triv}")

# Step 5: Now for non-trivial eigenspaces
# Under full Aut = AGL(1, F_11), the character group is more complex.
# The cyclic part Z/11 gives eigenvalues omega^k (k=0..10).
# The QR scaling part permutes these: k -> a*k mod 11 for a in QR.
# The QR orbits on {1..10} are:
# {1,3,4,5,9} and {2,6,7,8,10} (the two QR/QNR cosets)
# Wait: QR = {1,3,4,5,9}. QR acts on Z/11* by multiplication.
# k=1 orbit: {1*1, 3*1, 4*1, 5*1, 9*1} = {1,3,4,5,9} = QR itself
# k=2 orbit: {1*2, 3*2, 4*2, 5*2, 9*2} = {2,6,8,10,7} = QNR
# So under full Aut, the 10 non-trivial cyclic eigenspaces split into
# 2 groups of 5: QR-type (k in QR) and QNR-type (k in QNR).

# For Aut = Z/p x (Z/p)* acting on eigenspace k, the scaling part
# maps eigenspace k to eigenspace ak. So eigenspaces within each orbit
# are isomorphic.

# Therefore: total beta_8 = beta_8^(triv) + 5 * beta_8^(QR) + 5 * beta_8^(QNR)
# And since conjugation k -> -k = (p-1)*k maps QR to itself (since -1 in QR for p=11),
# actually QR eigenspaces have conjugation symmetry within themselves.

# For eigenspace k=1 (QR representative):
print(f"\n{'='*70}")
print("k=1 (QR-TYPE) EIGENSPACE")
print(f"{'='*70}")

# For eigenspace k under full Aut, the dimension per orbit depends on
# the orbit size under the cyclic part AND the QR scaling.
# Under just Z/11, each orbit has size 11, giving per-eigenspace dim = orbits/11.
# Under full Aut, orbits may have size 55, giving per-eigenspace dim = orbits/55.
# But we need to be careful: the eigenspace projection is only over Z/11.

# Actually, for Paley P_p, the full Aut orbits on paths have size divisible by p(p-1)/2.
# The eigenspace k has dimension |A_p|/p for each k.
# The QR scaling then acts WITHIN each eigenspace k, giving:
#   per-QR-eigenspace dim = |A_p| / (p * |QR|) = |A_p| / (p*(p-1)/2)

# Let's just compute the cyclic eigenspace k=1 directly using the orbit approach.
# We need B and J for k=1 in cyclic orbit coordinates.

def build_orbit_boundary_eigenspace(reps_p, paths_pm1_set, cyclic_orb_idx_pm1,
                                     cyclic_orb_shift_pm1, n, k):
    """Build boundary in CYCLIC orbit coordinates for eigenspace k.

    In cyclic orbit coordinates, for eigenvalue omega^k:
    - Each orbit rep r has boundary faces
    - Face at cyclic shift j from its orbit rep contributes omega^{-kj} * (-1)^i

    Returns: (B, J, n_junk)
    """
    omega = np.exp(2j * np.pi / n)
    n_orbits_p = len(reps_p)

    # Need cyclic orbit data for p-1 level
    n_orbits_pm1 = max(cyclic_orb_idx_pm1.values()) + 1 if cyclic_orb_idx_pm1 else 0

    junk_map = {}
    junk_reps = []

    # First pass: identify junk orbits (cyclic orbits)
    for r, rep in enumerate(reps_p):
        p_len = len(rep) - 1
        for i in range(p_len + 1):
            face = rep[:i] + rep[i+1:]
            if face not in paths_pm1_set:
                # Cyclic canon rep
                from p11_beta8_orbit import canon_orbit_rep as cyclic_canon
                jrep = cyclic_canon(face, n)
                if jrep not in junk_map:
                    junk_map[jrep] = len(junk_reps)
                    junk_reps.append(jrep)

    n_junk = len(junk_map)
    B = np.zeros((n_orbits_pm1, n_orbits_p), dtype=complex)
    J = np.zeros((n_junk, n_orbits_p), dtype=complex)

    for r, rep in enumerate(reps_p):
        p_len = len(rep) - 1
        for i in range(p_len + 1):
            face = rep[:i] + rep[i+1:]
            sign = (-1.0)**i
            if face in paths_pm1_set:
                q = cyclic_orb_idx_pm1[face]
                j = cyclic_orb_shift_pm1[face]
                B[q, r] += sign * omega**(-k * j)
            else:
                from p11_beta8_orbit import canon_orbit_rep as cyclic_canon
                from p11_beta8_orbit import find_orbit_shift as cyclic_shift
                jrep = cyclic_canon(face, n)
                jidx = junk_map[jrep]
                j = cyclic_shift(face, jrep, n)
                if j is not None:
                    J[jidx, r] += sign * omega**(-k * j)

    return B, J, n_junk


# We need cyclic orbit data. Let's compute it.
print("\n--- Cyclic orbit computation ---")

def cyclic_canon_orbit_rep(path, n):
    best = tuple(path)
    for j in range(1, n):
        rotated = tuple((v + j) % n for v in path)
        if rotated < best:
            best = rotated
    return best

def cyclic_find_orbit_shift(path, rep, n):
    for j in range(n):
        shifted = tuple((v + j) % n for v in rep)
        if shifted == tuple(path):
            return j
    return None

cyclic_reps = {}
cyclic_orb_idx = {}
cyclic_orb_shift = {}

for p in [6, 7, 8, 9]:
    if p not in allowed or not allowed[p]:
        continue
    t1 = time.time()
    reps = []
    oidx = {}
    oshift = {}
    orbit_map = {}

    for pa in allowed[p]:
        pt = tuple(pa)
        rep = cyclic_canon_orbit_rep(pa, n)
        if rep not in orbit_map:
            orbit_map[rep] = len(reps)
            reps.append(rep)
        idx = orbit_map[rep]
        shift = cyclic_find_orbit_shift(pa, rep, n)
        oidx[pt] = idx
        oshift[pt] = shift

    cyclic_reps[p] = reps
    cyclic_orb_idx[p] = oidx
    cyclic_orb_shift[p] = oshift
    print(f"  dim {p}: {len(reps)} cyclic orbits ({time.time()-t1:.1f}s)")

# Build boundary for k=1 in cyclic orbit coords
print("\n--- Building cyclic-orbit boundary for k=1 ---")
omega = np.exp(2j * np.pi / n)

B_k1 = {}
J_k1 = {}
n_junk_k1 = {}

for p in [7, 8, 9]:
    if p not in cyclic_reps:
        continue
    pm1 = p - 1
    t1 = time.time()
    n_orbits_p = len(cyclic_reps[p])
    n_orbits_pm1 = len(cyclic_reps[pm1]) if pm1 in cyclic_reps else 0

    junk_map = {}
    junk_list = []

    B_rows, B_cols, B_vals = [], [], []
    J_entries = {}  # (jidx, col) -> val

    for r, rep in enumerate(cyclic_reps[p]):
        p_len = len(rep) - 1
        for i in range(p_len + 1):
            face = rep[:i] + rep[i+1:]
            sign = (-1.0)**i
            if face in path_sets.get(pm1, set()):
                q = cyclic_orb_idx[pm1][face]
                j = cyclic_orb_shift[pm1][face]
                coeff = sign * omega**(-1 * j)
                B_rows.append(q)
                B_cols.append(r)
                B_vals.append(coeff)
            else:
                jrep = cyclic_canon_orbit_rep(face, n)
                if jrep not in junk_map:
                    junk_map[jrep] = len(junk_list)
                    junk_list.append(jrep)
                jidx = junk_map[jrep]
                j = cyclic_find_orbit_shift(face, jrep, n)
                if j is not None:
                    coeff = sign * omega**(-1 * j)
                    key = (jidx, r)
                    J_entries[key] = J_entries.get(key, 0) + coeff

    nj = len(junk_list)

    # Build sparse B
    B = sparse.csr_matrix((B_vals, (B_rows, B_cols)),
                           shape=(n_orbits_pm1, n_orbits_p), dtype=complex)

    # Build sparse J
    if nj > 0:
        jr, jc, jv = [], [], []
        for (ji, ci), val in J_entries.items():
            jr.append(ji)
            jc.append(ci)
            jv.append(val)
        J = sparse.csr_matrix((jv, (jr, jc)),
                               shape=(nj, n_orbits_p), dtype=complex)
    else:
        J = sparse.csr_matrix((nj, n_orbits_p), dtype=complex)

    B_k1[p] = B
    J_k1[p] = J
    n_junk_k1[p] = nj
    print(f"  d_{p}: B=({B.shape}), J=({J.shape}) ({time.time()-t1:.1f}s)")

# Compute beta_8 for k=1
# The cyclic orbit matrices are still large: 14395 columns for dim 8.
# But with sparse J, computing rank should be tractable via
# J^H J (which is n_orbits_p x n_orbits_p = 14395 x 14395 complex)
# OR via random projection if needed.

print(f"\n--- Computing beta_8^(k=1) ---")

# Use J^H @ J approach with eigvalsh for rank
t1 = time.time()
J8 = J_k1[8]
n_orb8 = len(cyclic_reps[8])  # 14395

print(f"  J_8 shape: {J8.shape}, nnz={J8.nnz}")
print(f"  Computing J_8^H @ J_8 ({n_orb8}x{n_orb8})...")

# Use sparse J^H @ J
JHJ = (J8.conj().T @ J8).toarray()
print(f"  JHJ computed, shape {JHJ.shape} ({time.time()-t1:.1f}s)")

t1 = time.time()
print(f"  Computing eigenvalues of JHJ...")
eigenvals = np.linalg.eigvalsh(JHJ)
rk_J8 = int(np.sum(np.abs(eigenvals) > 1e-12))
om8_dim = n_orb8 - rk_J8
print(f"  Om_8^(k=1) = {om8_dim} (rk J = {rk_J8}) ({time.time()-t1:.1f}s)")

if om8_dim == 0:
    beta8_k1 = 0
    print(f"\n  beta_8^(k=1) = 0")
else:
    # Get Om_8 basis from JHJ kernel
    t1 = time.time()
    eigenvals_full, eigvecs = np.linalg.eigh(JHJ)
    mask = np.abs(eigenvals_full) < 1e-10
    om8_basis = eigvecs[:, mask]
    print(f"  Om_8 basis: {om8_basis.shape} ({time.time()-t1:.1f}s)")

    # d_8 restricted to Om_8
    t1 = time.time()
    B8 = B_k1[8].toarray()
    bd8_om = B8 @ om8_basis
    S_vals_bd = np.linalg.svd(bd8_om, compute_uv=False)
    rk_d8 = int(np.sum(np.abs(S_vals_bd) > 1e-6))
    ker_d8 = om8_dim - rk_d8
    print(f"  ker(d_8 on Om_8) = {ker_d8} ({time.time()-t1:.1f}s)")

    # Om_9 and im(d_9)
    t1 = time.time()
    J9 = J_k1[9]
    n_orb9 = len(cyclic_reps[9])
    print(f"  J_9 shape: {J9.shape}, nnz={J9.nnz}")

    JHJ9 = (J9.conj().T @ J9).toarray()
    eigenvals9 = np.linalg.eigvalsh(JHJ9)
    rk_J9 = int(np.sum(np.abs(eigenvals9) > 1e-12))
    om9_dim = n_orb9 - rk_J9
    print(f"  Om_9^(k=1) = {om9_dim} (rk J = {rk_J9}) ({time.time()-t1:.1f}s)")

    if om9_dim == 0:
        im_d9 = 0
    else:
        eigenvals9_full, eigvecs9 = np.linalg.eigh(JHJ9)
        mask9 = np.abs(eigenvals9_full) < 1e-10
        om9_basis = eigvecs9[:, mask9]
        B9 = B_k1[9].toarray()
        bd9_om = B9 @ om9_basis
        S_vals_im = np.linalg.svd(bd9_om, compute_uv=False)
        im_d9 = int(np.sum(np.abs(S_vals_im) > 1e-6))

    beta8_k1 = max(0, ker_d8 - im_d9)
    print(f"\n  im(d_9) = {im_d9}")
    print(f"  beta_8^(k=1) = {beta8_k1}")

# Compute for k=2 (QNR representative) similarly
print(f"\n--- Computing beta_8^(k=2) ---")

B_k2 = {}
J_k2 = {}

for p in [8, 9]:
    if p not in cyclic_reps:
        continue
    pm1 = p - 1
    t1 = time.time()
    n_orbits_p = len(cyclic_reps[p])
    n_orbits_pm1 = len(cyclic_reps[pm1]) if pm1 in cyclic_reps else 0

    junk_map = {}
    junk_list = []
    B_rows, B_cols, B_vals = [], [], []
    J_entries = {}

    for r, rep in enumerate(cyclic_reps[p]):
        p_len = len(rep) - 1
        for i in range(p_len + 1):
            face = rep[:i] + rep[i+1:]
            sign = (-1.0)**i
            if face in path_sets.get(pm1, set()):
                q = cyclic_orb_idx[pm1][face]
                j = cyclic_orb_shift[pm1][face]
                coeff = sign * omega**(-2 * j)
                B_rows.append(q)
                B_cols.append(r)
                B_vals.append(coeff)
            else:
                jrep = cyclic_canon_orbit_rep(face, n)
                if jrep not in junk_map:
                    junk_map[jrep] = len(junk_list)
                    junk_list.append(jrep)
                jidx = junk_map[jrep]
                j = cyclic_find_orbit_shift(face, jrep, n)
                if j is not None:
                    coeff = sign * omega**(-2 * j)
                    key = (jidx, r)
                    J_entries[key] = J_entries.get(key, 0) + coeff

    nj = len(junk_list)
    B = sparse.csr_matrix((B_vals, (B_rows, B_cols)),
                           shape=(n_orbits_pm1, n_orbits_p), dtype=complex)
    if nj > 0:
        jr, jc, jv = [], [], []
        for (ji, ci), val in J_entries.items():
            jr.append(ji)
            jc.append(ci)
            jv.append(val)
        J = sparse.csr_matrix((jv, (jr, jc)),
                               shape=(nj, n_orbits_p), dtype=complex)
    else:
        J = sparse.csr_matrix((nj, n_orbits_p), dtype=complex)

    B_k2[p] = B
    J_k2[p] = J
    print(f"  d_{p}: B=({B.shape}), J=({J.shape}) ({time.time()-t1:.1f}s)")

# Compute beta_8^(k=2)
t1 = time.time()
J8_k2 = J_k2[8]
n_orb8_k2 = len(cyclic_reps[8])
JHJ_k2 = (J8_k2.conj().T @ J8_k2).toarray()
eigenvals_k2 = np.linalg.eigvalsh(JHJ_k2)
rk_J8_k2 = int(np.sum(np.abs(eigenvals_k2) > 1e-12))
om8_dim_k2 = n_orb8_k2 - rk_J8_k2
print(f"  Om_8^(k=2) = {om8_dim_k2} ({time.time()-t1:.1f}s)")

if om8_dim_k2 == 0:
    beta8_k2 = 0
else:
    eigenvals_k2_f, eigvecs_k2 = np.linalg.eigh(JHJ_k2)
    mask_k2 = np.abs(eigenvals_k2_f) < 1e-10
    om8_basis_k2 = eigvecs_k2[:, mask_k2]
    B8_k2 = B_k2[8].toarray()
    bd8_om_k2 = B8_k2 @ om8_basis_k2
    S_vals_bd_k2 = np.linalg.svd(bd8_om_k2, compute_uv=False)
    rk_d8_k2 = int(np.sum(np.abs(S_vals_bd_k2) > 1e-6))
    ker_d8_k2 = om8_dim_k2 - rk_d8_k2

    J9_k2 = J_k2[9]
    JHJ9_k2 = (J9_k2.conj().T @ J9_k2).toarray()
    eigenvals9_k2 = np.linalg.eigvalsh(JHJ9_k2)
    rk_J9_k2 = int(np.sum(np.abs(eigenvals9_k2) > 1e-12))
    om9_dim_k2 = n_orb9 - rk_J9_k2

    if om9_dim_k2 == 0:
        im_d9_k2 = 0
    else:
        eigenvals9_k2_f, eigvecs9_k2 = np.linalg.eigh(JHJ9_k2)
        mask9_k2 = np.abs(eigenvals9_k2_f) < 1e-10
        om9_basis_k2 = eigvecs9_k2[:, mask9_k2]
        B9_k2 = B_k2[9].toarray()
        bd9_om_k2 = B9_k2 @ om9_basis_k2
        S_vals_im_k2 = np.linalg.svd(bd9_om_k2, compute_uv=False)
        im_d9_k2 = int(np.sum(np.abs(S_vals_im_k2) > 1e-6))

    beta8_k2 = max(0, ker_d8_k2 - im_d9_k2)

print(f"  beta_8^(k=2) = {beta8_k2}")

# Summary
print(f"\n{'='*70}")
print("SUMMARY")
print("=" * 70)
print(f"  beta_8^(triv, k=0) = {beta8_triv}")
print(f"  beta_8^(QR, k=1) = {beta8_k1} (same for k in {{1,3,4,5,9}})")
print(f"  beta_8^(QNR, k=2) = {beta8_k2} (same for k in {{2,6,7,8,10}})")
total = beta8_triv + 5 * beta8_k1 + 5 * beta8_k2
print(f"  TOTAL beta_8 = {beta8_triv} + 5*{beta8_k1} + 5*{beta8_k2} = {total}")
print(f"\n  HYP-212 prediction: beta_8 = 10 (= p-1 + 0 for prime p)")
print(f"  Matches: {total == 10}")
if beta8_k1 == beta8_k2:
    print(f"  QR = QNR symmetry: YES (both give {beta8_k1})")
else:
    print(f"  QR = QNR symmetry: NO ({beta8_k1} vs {beta8_k2})")

print(f"\nTotal time: {time.time()-t0:.1f}s")
