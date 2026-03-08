#!/usr/bin/env python3
"""
beta2_z2_identity.py - Derive and verify the Z2 dimension formula

Chain complex: Om3 -bd3-> Om2 -bd2-> Om1 -bd1-> Om0

Known:
- Om1 = A1 = C(n,2), rk(bd1) = n-1
- Z1 = ker(bd1) dim = C(n-1,2)
- beta1 = c3, so rk(bd2|_Om2) = C(n-1,2) - c3
- Om2 = |A2| - J2

Therefore: Z2 = Om2 - rk(bd2) = |A2| - J2 - C(n-1,2) + c3
For beta2=0: rk(bd3|_Om3) = |A2| - J2 - C(n-1,2) + c3

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, time
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


def count_c3(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]):
                    c3 += 1
    return c3


def compute_J2(A, n):
    J2 = 0
    for a in range(n):
        for c in range(n):
            if a == c or not A[c][a]:
                continue
            for b in range(n):
                if b != a and b != c and A[a][b] and A[b][c]:
                    J2 += 1
                    break
    return J2


# ============================================================
# Verify Z2 formula and beta2=0 at n=5 and n=6
# ============================================================
for n in [5, 6]:
    print("=" * 70)
    print(f"Z2 FORMULA AND BETA2=0 VERIFICATION: n={n}")
    print("=" * 70)

    total = 1 << (n*(n-1)//2)
    cn2 = (n-1)*(n-2)//2

    z2_mismatches = 0
    bd3_mismatches = 0
    t0 = time.time()

    for bits in range(total):
        if n == 6 and bits % 10000 == 0 and bits > 0:
            dt = time.time() - t0
            print(f"  ... {bits}/{total} ({dt:.0f}s)")

        A = build_adj(n, bits)
        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)

        om2 = compute_omega_basis(A, n, 2, a2, a1)
        om3 = compute_omega_basis(A, n, 3, a3, a2)
        d_om2 = om2.shape[1] if om2.ndim == 2 else 0
        d_om3 = om3.shape[1] if om3.ndim == 2 else 0

        bd2 = build_full_boundary_matrix(a2, a1)
        bd3 = build_full_boundary_matrix(a3, a2)

        rk_bd2 = 0
        if d_om2 > 0:
            S2 = np.linalg.svd(bd2 @ om2, compute_uv=False)
            rk_bd2 = sum(s > 1e-8 for s in S2)
        dim_z2 = d_om2 - rk_bd2

        rk_bd3 = 0
        if d_om3 > 0:
            S3 = np.linalg.svd(bd3 @ om3, compute_uv=False)
            rk_bd3 = sum(s > 1e-8 for s in S3)

        c3 = count_c3(A, n)
        j2 = compute_J2(A, n)
        pred_z2 = len(a2) - j2 - cn2 + c3

        if pred_z2 != dim_z2:
            z2_mismatches += 1
        if rk_bd3 != dim_z2:
            bd3_mismatches += 1

    dt = time.time() - t0
    print(f"\nn={n}: {total} tournaments in {dt:.0f}s")
    print(f"Z2 = |A2| - J2 - C(n-1,2) + c3 mismatches: {z2_mismatches}/{total}")
    print(f"beta2=0 (rk_bd3 = dim_Z2) mismatches: {bd3_mismatches}/{total}")
    if z2_mismatches == 0 and bd3_mismatches == 0:
        print(f"  ** BOTH CONFIRMED at n={n}! **")


# ============================================================
# n=7 sampled
# ============================================================
print(f"\n{'='*70}")
print("VERIFICATION: n=7 (sampled)")
print("=" * 70)

import random
random.seed(42)
n = 7
cn2 = (n-1)*(n-2)//2
n_samples = 3000

z2_m = 0
bd3_m = 0
t0 = time.time()
for trial in range(n_samples):
    bits = random.randint(0, (1 << (n*(n-1)//2)) - 1)
    A = build_adj(n, bits)
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    bd2 = build_full_boundary_matrix(a2, a1)
    bd3 = build_full_boundary_matrix(a3, a2)

    rk_bd2 = 0
    if d_om2 > 0:
        S2 = np.linalg.svd(bd2 @ om2, compute_uv=False)
        rk_bd2 = sum(s > 1e-8 for s in S2)

    rk_bd3 = 0
    if d_om3 > 0:
        S3 = np.linalg.svd(bd3 @ om3, compute_uv=False)
        rk_bd3 = sum(s > 1e-8 for s in S3)

    dim_z2 = d_om2 - rk_bd2
    c3 = count_c3(A, n)
    j2 = compute_J2(A, n)
    pred_z2 = len(a2) - j2 - cn2 + c3

    if pred_z2 != dim_z2:
        z2_m += 1
    if rk_bd3 != dim_z2:
        bd3_m += 1

    if trial % 1000 == 999:
        dt = time.time() - t0
        print(f"  ... {trial+1}/{n_samples} ({dt:.0f}s)")

dt = time.time() - t0
print(f"\nn=7: {n_samples} samples in {dt:.0f}s")
print(f"Z2 formula mismatches: {z2_m}/{n_samples}")
print(f"beta2=0 mismatches: {bd3_m}/{n_samples}")


# ============================================================
# Now examine what |A2| - J2 actually counts
# ============================================================
print(f"\n{'='*70}")
print("COMBINATORIAL MEANING OF |A2| - J2")
print("=" * 70)

print("""
|A2| = # of allowed 2-paths (a,b,c) with a->b->c, all distinct

J2 = # of pairs (a,c) with c->a AND some b with a->b->c
   = # of "junk pairs" where the backward edge c->a is bridged

A 2-path (a,b,c) has boundary:
  bd(a,b,c) = (b,c) - (a,c) + (a,b)

The face (a,c) is NOT an allowed 1-path iff c->a (backward).
So Om2 = {u in A2 : for each backward pair (a,c), sum of coefficients = 0}

Om2 = ker(P) where P is the junk projection matrix.
Each junk pair (a,c) gives one row: P[(a,c), (a,b,c)] = -1 for each b.
(Actually the sign is (-1)^1 = -1 from the d_1 face.)

dim(Om2) = |A2| - rk(P) = |A2| - J2
(Since each junk pair has at least one path, and the rows are independent.)
""")

# Let's understand what rk(bd2|_Om2) = C(n-1,2) - c3 means
print("rk(bd2|_Om2) counts:")
print("  = dim(im bd2 in Z1) = dim(B1)")
print("  = dim(Z1) - beta1 = C(n-1,2) - c3")
print()
print("Z1 has a nice basis: for each pair {i,j}, the 'cycle difference'")
print("In a tournament, Z1 = ker(bd1) where bd1: edge (i,j) -> j - i")
print("dim(Z1) = C(n,2) - (n-1) = C(n-1,2)")
print()
print("beta1 = c3 means: each 3-cycle contributes one independent")
print("element of Z1 that is NOT a boundary of any Om2 element.")


# ============================================================
# Euler characteristic decomposition
# ============================================================
print(f"\n{'='*70}")
print("EULER CHARACTERISTIC")
print("=" * 70)

n = 5
print(f"\nn={n}:")
for bits in [0, 1, 10, 100, 500]:
    A = build_adj(n, bits)
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a4 = enumerate_allowed_paths(A, n, 4)

    c3 = count_c3(A, n)
    j2 = compute_J2(A, n)
    scores = tuple(sorted([sum(row) for row in A]))

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    om4 = compute_omega_basis(A, n, 4, a4, a3)
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0
    d_om4 = om4.shape[1] if om4.ndim == 2 else 0

    # chi = beta0 - beta1 + beta2 - ...
    # If beta2=0: chi = 1 - c3 + 0 - beta3 + ...
    # Also chi = sum (-1)^p dim(Om_p) - ... no that's wrong

    print(f"  bits={bits}, scores={scores}, c3={c3}")
    print(f"    Om: {n}, {len(a1)}, {d_om2}, {d_om3}, {d_om4}")
    print(f"    A:  {n}, {len(a1)}, {len(a2)}, {len(a3)}, {len(a4)}")


print("\n\nDone.")
