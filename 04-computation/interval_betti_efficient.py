#!/usr/bin/env python3
"""
interval_betti_efficient.py — Memory-efficient Betti for Interval T_11.

Avoids constructing the full BD matrix (14937 x 8457 = 1GB).
Instead, computes BD @ omega_m_basis.T column-by-column.

Author: opus-2026-03-12-S68
"""
import sys, time, gc
import numpy as np
sys.path.insert(0, '04-computation')
from circulant_homology import CirculantHomology, find_nth_root_of_unity, _gauss_rank

p = 11
h = CirculantHomology(n=p, S={1,2,3,4,5})
h._ensure_enumerated(p)
prime = h.prime
omega_k = 1  # k=0

def boundary_rank_efficient(h, m, omega_k, prime):
    """Compute boundary rank without full BD matrix."""
    A_m = h._diff_seqs.get(m, [])
    A_m1 = h._diff_seqs.get(m - 1, [])
    if not A_m or not A_m1:
        return 0

    # Get omega bases
    omega_m_basis, face_data_m = h._omega_basis_k(m, omega_k)
    if omega_m_basis.shape[0] == 0:
        return 0

    dim_Om = omega_m_basis.shape[0]  # = Omega_m

    if m == 1:
        omega_m1_basis = np.ones((1, 1), dtype=np.int64)
    else:
        omega_m1_basis, _ = h._omega_basis_k(m - 1, omega_k)

    dim_Om1 = omega_m1_basis.shape[0] if omega_m1_basis.ndim == 2 else 0
    if dim_Om1 == 0:
        return dim_Om

    A_m1_idx = {d: i for i, d in enumerate(A_m1)}
    n_Am1 = len(A_m1)

    # Compute BD @ omega_m_basis.T WITHOUT storing BD.
    # Result shape: (n_Am1, dim_Om)
    # Process in chunks to manage memory
    bd_omega = np.zeros((n_Am1, dim_Om), dtype=np.int64)

    for j, (D, faces) in enumerate(zip(A_m, face_data_m)):
        omega_j = omega_m_basis[j]  # row j, shape (dim_Om,)
        for fd, sign, offset, is_allowed in faces:
            if is_allowed and fd in A_m1_idx:
                row = A_m1_idx[fd]
                w = pow(omega_k, offset, prime) if offset != 0 else 1
                entry = (sign * w) % prime
                bd_omega[row] = (bd_omega[row] + entry * omega_j) % prime

    del omega_m_basis
    gc.collect()

    # Now compute d_restricted = omega_m1_basis @ bd_omega
    # Shape: (dim_Om1, dim_Om) — small!
    d_restricted = omega_m1_basis @ bd_omega % prime % prime

    del bd_omega, omega_m1_basis
    gc.collect()

    return _gauss_rank(d_restricted, prime)


# Already known: r_0..r_7 = [0, 0, 4, 16, 58, 166, 356, 570]
known_ranks = {0: 0, 1: 0, 2: 4, 3: 16, 4: 58, 5: 166, 6: 356, 7: 570}

# Compute r_8, r_9, r_10
for m in [8, 9, 10]:
    A_m = h._diff_seqs.get(m, [])
    A_m1 = h._diff_seqs.get(m - 1, [])
    omega_m = h.omega_dims(max_degree=m, verbose=False)[-1]
    omega_m1 = h.omega_dims(max_degree=m-1, verbose=False)[-1]

    print(f'm={m}: |A_m|={len(A_m)}, |A_{m-1}|={len(A_m1)}, '
          f'Omega_m={omega_m}, Omega_{m-1}={omega_m1}', flush=True)

    t0 = time.time()
    rk = boundary_rank_efficient(h, m, omega_k, prime)
    dt = time.time() - t0
    known_ranks[m] = rk
    print(f'  r_{m}^(0) = {rk}  ({dt:.1f}s)', flush=True)
    gc.collect()

known_ranks[11] = 0  # No degree-11 paths

# Compute Betti
omega_dims = h.omega_dims(max_degree=p-1, verbose=False)
print(f'\nOmega: {omega_dims}')
print(f'Ranks: {[known_ranks[m] for m in range(p+1)]}')

betti = []
for m in range(p):
    if m == 0:
        b = 1  # β_0 = 1 always
    elif m == 1:
        b = 1  # β_1^(k=0) = 5 - 0 - 4 = 1, β_1^(k>=1) = 0
    else:
        per_k = omega_dims[m] - known_ranks[m] - known_ranks[m + 1]
        b = p * per_k
    betti.append(b)

print(f'\nInterval T_{p} Betti: {betti}')
chi = sum((-1)**i * b for i, b in enumerate(betti))
print(f'chi = {chi}')
print(f'\nPaley  T_{p} Betti:  [1, 0, 0, 0, 0, 5, 15, 0, 0, 0, 0]')

# Save to file
with open('05-knowledge/results/interval_betti_p11.out', 'w') as f:
    f.write(f'Interval T_{p}\n')
    f.write(f'Omega: {omega_dims}\n')
    f.write(f'Ranks (k=0): {[known_ranks[m] for m in range(p+1)]}\n')
    f.write(f'Betti: {betti}\n')
    f.write(f'chi = {chi}\n')
    f.write(f'\nPaley T_{p}: [1, 0, 0, 0, 0, 5, 15, 0, 0, 0, 0], chi=11\n')
print('\nDONE.')
