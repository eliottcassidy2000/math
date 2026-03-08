#!/usr/bin/env python3
"""
beta2_rank_identity.py - Study the rank identity rk(d_2) + rk(d_3) = dim(Omega_2)

For beta_2 = 0 iff rk(d_2|O2) + rk(d_3|O3) = dim(O2).

This is the key identity to prove. Let's understand what determines rk(d_2) and rk(d_3)
separately, and why their sum equals dim(O2).

rk(d_2) = dim(O2) - dim(Z_2) = dim(im(d_2|O2))
rk(d_3) = dim(im(d_3|O3)) = dim(O3) - dim(ker(d_3|O3)) = dim(O3) - beta_3

So: rk(d_2) + rk(d_3) = dim(O2) - Z_2 + dim(O3) - beta_3
   = dim(O2) + surplus - beta_3  (since surplus = O3 - Z_2)

For this to equal dim(O2): surplus = beta_3.

So beta_2 = 0 iff surplus = beta_3!

Surplus = dim(O3) - dim(Z_2) = dim(O3) - Z_2
Beta_3 = dim(ker d_3|O3) = dim(O3) - rk(d_3)

surplus = beta_3 iff dim(O3) - Z_2 = dim(O3) - rk(d_3) iff Z_2 = rk(d_3) = im(d_3).
Which is just beta_2 = Z_2 - im(d_3) = 0. So this is circular.

But we learn: surplus = beta_3 + beta_2. So surplus >= beta_3 always,
and surplus >= 0 iff beta_2 + beta_3 >= 0 (which is trivial since dimensions are non-negative).

Wait: surplus = O3 - Z_2. And Z_2 = rk(d_3) + beta_2. And rk(d_3) = O3 - beta_3.
So surplus = O3 - (O3 - beta_3 + beta_2) = beta_3 - beta_2.

So SURPLUS = BETA_3 - BETA_2.

This means:
- surplus >= 0 iff beta_3 >= beta_2
- beta_2 = 0 iff surplus = beta_3
- surplus = 0 iff beta_3 = beta_2

This is a HUGE insight! The question "is beta_2 = 0 for all tournaments?"
is equivalent to "is surplus = beta_3?"

And "is surplus >= 0?" is equivalent to "is beta_3 >= beta_2?"

Let me verify this identity computationally.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from collections import Counter
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

def compute_bettis(A, n, max_dim=4):
    allowed = {}
    for p in range(max_dim + 2):
        allowed[p] = enumerate_allowed_paths(A, n, p)
        if not allowed[p]:
            break

    omega_basis = {}
    for p in range(max_dim + 1):
        if p not in allowed or not allowed[p]:
            omega_basis[p] = np.zeros((0, 0))
            continue
        basis = compute_omega_basis(A, n, p, allowed[p],
                                     allowed[p-1] if p >= 1 and p-1 in allowed else [])
        omega_basis[p] = basis

    dims = {p: omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0
            for p in range(max_dim + 1)}

    betti = {}
    rk = {}
    for p in range(1, max_dim + 1):
        if dims[p] == 0:
            rk[p] = 0
            continue
        bd = build_full_boundary_matrix(allowed[p], allowed.get(p-1, []))
        bd_om = bd @ omega_basis[p]
        S_v = np.linalg.svd(bd_om, compute_uv=False)
        rk[p] = int(np.sum(np.abs(S_v) > 1e-8))

    # ker(d_p) = dims[p] - rk[p]
    # betti[p] = ker(d_p) - rk[p+1] = dims[p] - rk[p] - rk.get(p+1, 0)
    for p in range(max_dim + 1):
        ker_p = dims[p] - rk.get(p, 0) if p >= 1 else dims[p]
        im_p1 = rk.get(p+1, 0)
        betti[p] = ker_p - im_p1

    return betti, dims, rk


print("=" * 70)
print("SURPLUS = BETA_3 - BETA_2 VERIFICATION")
print("=" * 70)

for n in [5, 6]:
    print(f"\n  n={n}:")
    n_arcs = n*(n-1)//2
    total = 1 << n_arcs
    violations = 0

    for bits in range(total):
        if n == 6 and bits % 10000 == 0 and bits > 0:
            print(f"    ... {bits}/{total}")

        A = build_adj(n, bits)
        betti, dims, rk = compute_bettis(A, n, max_dim=4)

        surplus = dims.get(3, 0) - (dims.get(2, 0) - rk.get(2, 0))
        expected = betti.get(3, 0) - betti.get(2, 0)

        if surplus != expected:
            violations += 1
            print(f"    VIOLATION at bits={bits}: surplus={surplus}, beta3-beta2={expected}")
            print(f"      dims={dims}, rk={rk}, betti={betti}")

    print(f"    surplus = beta_3 - beta_2: {violations} violations out of {total}")

    # Distribution of beta_2, beta_3, surplus
    beta_dist = Counter()
    for bits in range(total):
        A = build_adj(n, bits)
        betti, dims, rk = compute_bettis(A, n, max_dim=4)
        beta_dist[(betti.get(2, 0), betti.get(3, 0))] += 1

    print(f"\n    (beta_2, beta_3) distribution:")
    for key in sorted(beta_dist.keys()):
        surplus = key[1] - key[0]
        print(f"      beta_2={key[0]}, beta_3={key[1]}: {beta_dist[key]} (surplus={surplus})")

print("\nDone.")
