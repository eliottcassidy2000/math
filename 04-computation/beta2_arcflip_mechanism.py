#!/usr/bin/env python3
"""
beta2_arcflip_mechanism.py - Understand WHY delta(Z_2) = delta(rk d_3)

Strategy: For each arc flip u->v to v->u, compute:
1. The exact change in Omega_2 (TT triples gained/lost) -- KNOWN formula
2. The exact change in Z_2 -- need to understand
3. The exact change in rk(d_2) -- since dim(Z_2) = dim(Om_2) - rk(d_2)
4. The exact change in rk(d_3)

Key identity we want to prove:
  delta(rk d_2) + delta(rk d_3) = delta(dim Omega_2)

Since delta(dim Z_2) = delta(dim Om_2) - delta(rk d_2),
and delta(rk d_3) = delta(dim Z_2) (our observation),
this is equivalent to: delta(rk d_2) + delta(dim Z_2) = delta(dim Om_2)
which is just the definition. So the nontrivial content is:
  delta(rk d_3) = delta(dim Om_2) - delta(rk d_2)

Can we find closed forms for delta(rk d_2)?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def compute_full_data(A, n):
    """Compute full exactness data including beta_1."""
    paths = {}
    omega = {}
    for p in range(5):
        paths[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            omega[p] = np.eye(n)
        elif len(paths[p]) > 0 and len(paths[p-1]) > 0:
            omega[p] = compute_omega_basis(A, n, p, paths[p], paths[p-1])
        else:
            omega[p] = np.zeros((max(1, len(paths[p])), 0))

    dims = {}
    rks = {}
    for p in range(5):
        dims[p] = omega[p].shape[1] if omega[p].ndim == 2 else 0

    for p in range(1, 5):
        if dims[p] > 0 and dims[p-1] > 0:
            bd = build_full_boundary_matrix(paths[p], paths[p-1])
            bd_om = bd @ omega[p]
            if p >= 2:
                im_coords, _, _, _ = np.linalg.lstsq(omega[p-1], bd_om, rcond=None)
                rks[p] = np.linalg.matrix_rank(im_coords, tol=1e-8)
            else:
                rks[p] = np.linalg.matrix_rank(bd_om, tol=1e-8)
        else:
            rks[p] = 0

    betas = {}
    for p in range(4):
        z_p = dims[p] - rks.get(p, 0)
        b_p = z_p - rks.get(p+1, 0)
        betas[p] = b_p

    return dims, rks, betas


def flip_arc(bits, i, j, n):
    idx = 0
    for a in range(n):
        for b in range(a+1, n):
            if a == i and b == j:
                return bits ^ (1 << idx)
            idx += 1
    return bits


print("=" * 70)
print("ARC-FLIP MECHANISM ANALYSIS")
print("=" * 70)

n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

print(f"\nPrecomputing n={n}...")
t0 = time.time()
all_data = {}
for bits in range(total):
    A = build_adj(n, bits)
    all_data[bits] = compute_full_data(A, n)
dt = time.time() - t0
print(f"  Done in {dt:.1f}s")

# Detailed delta analysis
print(f"\n--- Detailed delta analysis ---")

delta_tuples = Counter()
delta_beta1_dist = Counter()

for bits in range(total):
    A = build_adj(n, bits)
    d, r, b = all_data[bits]

    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc(bits, i, j, n)
            d2, r2, b2 = all_data[bits2]

            dOm2 = d2[2] - d[2]
            dOm3 = d2[3] - d[3]
            drk2 = r2[2] - r[2]
            drk3 = r2[3] - r[3]
            dZ2 = dOm2 - drk2
            db1 = b2[1] - b[1]
            db2 = b2[2] - b[2]

            delta_tuples[(dOm2, drk2, dZ2, dOm3, drk3, db1)] += 1
            delta_beta1_dist[db1] += 1

print("\nDelta tuple distribution (dOm2, drk2, dZ2, dOm3, drk3, db1):")
for key, count in sorted(delta_tuples.items()):
    dOm2, drk2, dZ2, dOm3, drk3, db1 = key
    b2_check = dZ2 - drk3
    print(f"  ({dOm2:+d}, {drk2:+d}, {dZ2:+d}, {dOm3:+d}, {drk3:+d}, {db1:+d}): "
          f"count={count}, db2={b2_check}")

print(f"\ndelta(beta_1) distribution: {dict(sorted(delta_beta1_dist.items()))}")

# KEY: dim(Z_1) is constant for all tournaments on n vertices
print(f"\n{'='*70}")
print("KEY INSIGHT: dim(Z_1) is constant for tournaments")
print("=" * 70)
print("dim(Omega_1) = n(n-1)/2 (all arcs)")
print("rk(d_1) = n-1 (tournament is strongly connected <=> beta_0=1, always)")
print("dim(Z_1) = n(n-1)/2 - (n-1) = (n-1)(n-2)/2")
print()
print("Since rk(d_2) = dim(Z_1) - beta_1:")
print("  delta(rk d_2) = -delta(beta_1)")

# Verify
errors_rk2 = 0
for bits in range(total):
    d, r, b = all_data[bits]
    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc(bits, i, j, n)
            d2, r2, b2 = all_data[bits2]
            drk2 = r2[2] - r[2]
            db1 = b2[1] - b[1]
            if drk2 != -db1:
                errors_rk2 += 1

print(f"Verification: delta(rk d_2) = -delta(beta_1) errors: {errors_rk2}")

# So beta_2 = 0 invariance <=> delta(rk d_3) = delta(dim Om_2) + delta(beta_1)
print(f"\n{'='*70}")
print("REFORMULATION: beta_2 invariance <=> delta(rk d_3) = delta(dim Om_2) + delta(beta_1)")
print("=" * 70)

errors_formula = 0
for bits in range(total):
    d, r, b = all_data[bits]
    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc(bits, i, j, n)
            d2, r2, b2 = all_data[bits2]
            drk3 = r2[3] - r[3]
            dOm2 = d2[2] - d[2]
            db1 = b2[1] - b[1]
            if drk3 != dOm2 + db1:
                errors_formula += 1

print(f"  Errors: {errors_formula}")

# Now check: is delta(beta_1) determined by local data?
print(f"\n{'='*70}")
print("Is delta(beta_1) determined by local data?")
print("=" * 70)

db1_by_local = defaultdict(list)
for bits in range(total):
    A = build_adj(n, bits)
    d, r, b = all_data[bits]
    scores = [sum(row) for row in A]

    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc(bits, i, j, n)
            d2, r2, b2 = all_data[bits2]
            db1 = b2[1] - b[1]

            if A[i][j] == 1:
                u, v = i, j
            else:
                u, v = j, i

            cout = sum(1 for w in range(n) if w != u and w != v and A[u][w] and A[v][w])
            cin = sum(1 for w in range(n) if w != u and w != v and A[w][u] and A[w][v])
            p_uv = sum(1 for w in range(n) if w != u and w != v and A[u][w] and A[w][v])
            p_vu = sum(1 for w in range(n) if w != u and w != v and A[v][w] and A[w][u])

            db1_by_local[(scores[u], scores[v], cout, cin, p_uv, p_vu)].append(db1)

non_det = 0
for key in sorted(db1_by_local.keys()):
    vals = db1_by_local[key]
    val_set = set(vals)
    if len(val_set) > 1:
        non_det += 1
        if non_det <= 5:
            print(f"  {key}: values={dict(Counter(vals))}")

print(f"\n  Non-deterministic: {non_det}/{len(db1_by_local)}")

# Euler characteristic
print(f"\n{'='*70}")
print("EULER CHARACTERISTIC")
print("=" * 70)
chi_dist = Counter()
for bits in range(total):
    d, r, b = all_data[bits]
    chi = sum((-1)**p * b[p] for p in range(4))
    chi_dist[chi] += 1
print(f"  chi distribution: {dict(sorted(chi_dist.items()))}")

if len(chi_dist) == 1:
    chi_val = list(chi_dist.keys())[0]
    print(f"  chi is CONSTANT = {chi_val} for all tournaments!")
    print(f"  beta_0 - beta_1 + beta_2 - beta_3 = {chi_val}")
    print(f"  If beta_2 = 0: 1 - beta_1 - beta_3 = {chi_val}")
    print(f"  So: beta_1 + beta_3 = {1 - chi_val}")
else:
    print("  chi is NOT constant")

# Check at n=6
print(f"\n{'='*70}")
print("EULER CHARACTERISTIC at n=6,7 (sampled)")
print("=" * 70)
import random
random.seed(42)

for n_test in [6, 7]:
    n_arcs = n_test*(n_test-1)//2
    samples = 2000 if n_test == 6 else 500
    chi_dist = Counter()

    for _ in range(samples):
        bits = random.randint(0, (1 << n_arcs) - 1)
        A = build_adj(n_test, bits)
        _, _, b = compute_full_data(A, n_test)
        chi = sum((-1)**p * b.get(p, 0) for p in range(n_test))
        chi_dist[chi] += 1

    print(f"  n={n_test}: chi distribution ({samples} samples): {dict(sorted(chi_dist.items()))}")

print("\nDone.")
