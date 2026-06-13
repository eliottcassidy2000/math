#!/usr/bin/env python3
"""
beta2_filling_structure.py - Explicit Z_2 cycle filling analysis

For each tournament T at n=5, compute:
1. Explicit basis of Z_2 = ker(d_2) in Omega_2
2. For each Z_2 basis vector, find the Omega_3 chain that fills it
3. Analyze the structure of fillings to find patterns

Key identity to prove: rk(d_3) = dim(Z_2) for all tournaments.

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


def get_chain_data(A, n, max_p=4):
    """Get full chain complex data."""
    paths = {}
    omega = {}
    for p in range(max_p+1):
        paths[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            omega[p] = np.eye(n)
        elif len(paths[p]) > 0 and len(paths[p-1]) > 0:
            omega[p] = compute_omega_basis(A, n, p, paths[p], paths[p-1])
        else:
            omega[p] = np.zeros((max(1, len(paths[p])), 0))

    dims = {p: (omega[p].shape[1] if omega[p].ndim == 2 else 0) for p in range(max_p+1)}
    bds = {}
    bd_om = {}
    rks = {}

    for p in range(1, max_p+1):
        if len(paths.get(p, [])) > 0 and len(paths.get(p-1, [])) > 0:
            bds[p] = build_full_boundary_matrix(paths[p], paths[p-1])
        else:
            bds[p] = None

    for p in range(2, max_p+1):
        if dims[p] > 0 and bds.get(p) is not None and dims[p-1] > 0:
            raw = bds[p] @ omega[p]
            coords, _, _, _ = np.linalg.lstsq(omega[p-1], raw, rcond=None)
            bd_om[p] = coords
            rks[p] = np.linalg.matrix_rank(coords, tol=1e-8)
        else:
            bd_om[p] = None
            rks[p] = 0

    if bds.get(1) is not None and dims[1] > 0:
        rks[1] = np.linalg.matrix_rank(bds[1] @ omega[1], tol=1e-8)
    else:
        rks[1] = 0

    betas = {}
    for p in range(max_p):
        z_p = dims[p] - rks.get(p, 0)
        b_p = z_p - rks.get(p+1, 0)
        betas[p] = b_p

    return paths, omega, dims, bds, bd_om, rks, betas


print("=" * 70)
print("Z_2 CYCLE FILLING STRUCTURE ANALYSIS")
print("=" * 70)

n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

# Statistics
dimZ2_dist = Counter()
dimZ3_dist = Counter()
cycle_support_sizes = Counter()
filling_support_sizes = Counter()

# Find examples with different dim(Z_2)
examples_by_dimZ2 = {}

for bits in range(total):
    A = build_adj(n, bits)
    paths, omega, dims, bds, bd_om, rks, betas = get_chain_data(A, n)

    dimZ2 = dims[2] - rks.get(2, 0)
    dimZ3 = dims[3] - rks.get(3, 0)
    dimZ2_dist[dimZ2] += 1
    dimZ3_dist[dimZ3] += 1

    if dimZ2 not in examples_by_dimZ2:
        examples_by_dimZ2[dimZ2] = bits

print(f"dim(Z_2) distribution: {dict(sorted(dimZ2_dist.items()))}")
print(f"dim(Z_3) distribution: {dict(sorted(dimZ3_dist.items()))}")

# Detailed analysis for each dim(Z_2) type
for dimZ2_val in sorted(examples_by_dimZ2.keys()):
    bits = examples_by_dimZ2[dimZ2_val]
    A = build_adj(n, bits)
    paths, omega, dims, bds, bd_om, rks, betas = get_chain_data(A, n)

    dimZ2 = dims[2] - rks.get(2, 0)
    print(f"\n{'='*50}")
    print(f"Example: bits={bits}, dim(Z_2)={dimZ2}")
    print(f"  Scores: {[sum(row) for row in A]}")
    print(f"  dims: Om0={dims[0]}, Om1={dims[1]}, Om2={dims[2]}, Om3={dims[3]}, Om4={dims.get(4,0)}")
    print(f"  rks: d1={rks.get(1,0)}, d2={rks.get(2,0)}, d3={rks.get(3,0)}, d4={rks.get(4,0)}")
    print(f"  betas: {[betas.get(p,0) for p in range(4)]}")

    if dimZ2 == 0:
        print("  Z_2 is trivial.")
        continue

    if bd_om.get(2) is None:
        continue

    # Compute Z_2 basis
    U, S, Vt = np.linalg.svd(bd_om[2], full_matrices=True)
    rk2 = int(np.sum(np.abs(S) > 1e-8))
    Z2_coords = Vt[rk2:].T  # columns = Z_2 basis in Omega_2 coords

    # Show each Z_2 cycle
    for j in range(Z2_coords.shape[1]):
        z_om2 = Z2_coords[:, j]
        # Convert to path coordinates
        z_paths = omega[2] @ z_om2
        nonzero = [(i, z_paths[i]) for i in range(len(z_paths)) if abs(z_paths[i]) > 1e-8]

        print(f"\n  Z_2 cycle z_{j} ({len(nonzero)} nonzero paths):")
        vertex_set = set()
        for idx, coeff in sorted(nonzero, key=lambda x: abs(x[1]), reverse=True)[:8]:
            a, b, c = paths[2][idx]
            tt = "TT" if A[a][c] else "NT"
            print(f"    {coeff:+.4f} * ({a},{b},{c}) [{tt}]")
            vertex_set.update([a, b, c])
        if len(nonzero) > 8:
            print(f"    ... and {len(nonzero)-8} more")
        print(f"    Uses vertices: {sorted(vertex_set)}")

        # Find filling in Omega_3
        if bd_om.get(3) is not None and dims[3] > 0:
            d3_mat = bd_om[3]
            x, _, _, _ = np.linalg.lstsq(d3_mat, z_om2, rcond=None)
            check = np.max(np.abs(d3_mat @ x - z_om2))
            if check < 1e-8:
                # Convert to 3-path coordinates
                x_paths = omega[3] @ x
                nonzero3 = [(i, x_paths[i]) for i in range(len(x_paths)) if abs(x_paths[i]) > 1e-8]
                print(f"    Filled by Omega_3 chain ({len(nonzero3)} nonzero 3-paths):")
                vertex_set3 = set()
                for idx, coeff in sorted(nonzero3, key=lambda x: abs(x[1]), reverse=True)[:6]:
                    a, b, c, d = paths[3][idx]
                    face02 = "a->c" if A[a][c] else "c->a"
                    face13 = "b->d" if A[b][d] else "d->b"
                    print(f"      {coeff:+.4f} * ({a},{b},{c},{d}) [{face02}, {face13}]")
                    vertex_set3.update([a, b, c, d])
                if len(nonzero3) > 6:
                    print(f"      ... and {len(nonzero3)-6} more")
                print(f"      Uses vertices: {sorted(vertex_set3)}")
            else:
                print(f"    NOT fillable (err={check:.6f}) -- beta_2 > 0??")


# KEY ANALYSIS: Z_2 cycle types by vertex set
print(f"\n{'='*70}")
print("ANALYSIS: Z_2 cycle support structure")
print("=" * 70)

# For each tournament, decompose Z_2 into cycles supported on
# different vertex subsets
support_analysis = Counter()

for bits in range(total):
    A = build_adj(n, bits)
    paths, omega, dims, bds, bd_om, rks, betas = get_chain_data(A, n)

    dimZ2 = dims[2] - rks.get(2, 0)
    if dimZ2 == 0 or bd_om.get(2) is None:
        continue

    U, S, Vt = np.linalg.svd(bd_om[2], full_matrices=True)
    rk2 = int(np.sum(np.abs(S) > 1e-8))
    Z2_coords = Vt[rk2:].T

    # Check: are Z_2 cycles "localized" to vertex subsets?
    for j in range(Z2_coords.shape[1]):
        z_paths = omega[2] @ Z2_coords[:, j]
        vertex_set = set()
        for i, c in enumerate(z_paths):
            if abs(c) > 1e-8:
                vertex_set.update(paths[2][i])
        support_analysis[len(vertex_set)] += 1

print(f"Z_2 cycle vertex support sizes: {dict(sorted(support_analysis.items()))}")


# CRITICAL TEST: Is there a LOCALIZED filling?
# I.e., can each Z_2 cycle be filled by 3-paths on the SAME vertex set?
print(f"\n{'='*70}")
print("CRITICAL: Can each Z_2 cycle be filled within its vertex support?")
print("=" * 70)

local_fill_possible = 0
non_local_fill = 0

for bits in range(total):
    A = build_adj(n, bits)
    paths, omega, dims, bds, bd_om, rks, betas = get_chain_data(A, n)

    dimZ2 = dims[2] - rks.get(2, 0)
    if dimZ2 == 0 or bd_om.get(2) is None or bd_om.get(3) is None:
        continue

    U, S, Vt = np.linalg.svd(bd_om[2], full_matrices=True)
    rk2 = int(np.sum(np.abs(S) > 1e-8))
    Z2_coords = Vt[rk2:].T

    d3_mat = bd_om[3]

    for j in range(Z2_coords.shape[1]):
        z_om2 = Z2_coords[:, j]
        z_paths = omega[2] @ z_om2
        cycle_verts = set()
        for i, c in enumerate(z_paths):
            if abs(c) > 1e-8:
                cycle_verts.update(paths[2][i])

        # Find minimum-support filling
        x, _, _, _ = np.linalg.lstsq(d3_mat, z_om2, rcond=None)
        check = np.max(np.abs(d3_mat @ x - z_om2))
        if check > 1e-8:
            continue

        x_paths = omega[3] @ x
        fill_verts = set()
        for i, c in enumerate(x_paths):
            if abs(c) > 1e-8:
                fill_verts.update(paths[3][i])

        if fill_verts <= cycle_verts:
            local_fill_possible += 1
        else:
            non_local_fill += 1

print(f"  Locally fillable: {local_fill_possible}")
print(f"  Non-locally fillable: {non_local_fill}")
total_cycles_checked = local_fill_possible + non_local_fill
if total_cycles_checked > 0:
    print(f"  Fraction locally fillable: {local_fill_possible/total_cycles_checked:.4f}")


# ANALYSIS: dim(Omega_3) - dim(Z_2) = "surplus" dimension
print(f"\n{'='*70}")
print("SURPLUS ANALYSIS: dim(Omega_3) - dim(Z_2) for each tournament")
print("=" * 70)

surplus_dist = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    paths, omega, dims, bds, bd_om, rks, betas = get_chain_data(A, n)
    dimZ2 = dims[2] - rks.get(2, 0)
    surplus = dims[3] - dimZ2
    surplus_dist[surplus] += 1

print(f"Surplus distribution: {dict(sorted(surplus_dist.items()))}")

# The surplus tells us how much "room" there is for d_3 to fill Z_2.
# If surplus >= 0, there are at least as many Omega_3 generators as Z_2 cycles.
# For beta_2 = 0, we need rk(d_3) = dim(Z_2), which requires dim(Omega_3) >= dim(Z_2).
if all(s >= 0 for s in surplus_dist.keys()):
    print("  ALL surpluses >= 0: dim(Omega_3) >= dim(Z_2) always")
    print("  This is a NECESSARY condition for beta_2 = 0")
else:
    print("  WARNING: Some tournaments have dim(Omega_3) < dim(Z_2)!")


print("\nDone.")
