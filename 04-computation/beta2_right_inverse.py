#!/usr/bin/env python3
"""
beta2_right_inverse.py — Construct explicit right inverse for ∂₃: Ω₃ → Z₂

KEY IDEA: β₂=0 ⟺ ∃ linear map F: Z₂ → Ω₃ such that ∂₃∘F = id on Z₂.
F is a "right inverse" (section) of ∂₃.

Can we write F explicitly using tournament-local data?

APPROACH: For each A₂ path (a,b,c), define F(a,b,c) as a linear combination
of A₃ paths involving a,b,c plus one more vertex. The coefficients should
depend ONLY on the local tournament structure.

Test: F(a,b,c) = Σ_d w(a,b,c,d) · (a,b,c,d) where the weights w depend
on the sub-tournament T[{a,b,c,d}].

If such a "local" right inverse exists, it proves β₂=0.

Author: opus-2026-03-08-S49
"""
import sys
import numpy as np
from collections import Counter, defaultdict
from itertools import combinations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved


def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0


print("=" * 70)
print("RIGHT INVERSE STRUCTURE FOR ∂₃: Ω₃ → Z₂")
print("=" * 70)

# Strategy: For each tournament, compute the actual right inverse F
# and analyze its structure to find a pattern.

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

# For each tournament, extract the right inverse F in A₂ → A₃ form
# F: Z₂ → Ω₃, so F maps each Z₂ basis to an Ω₃ element.
# In coordinates: F_matrix has columns = Z₂ basis, rows = Ω₃ basis.
# The filling is Ω₃_basis @ F_matrix.

# Key question: is there a UNIVERSAL pattern for F that depends only on
# the local tournament structure?

# First: analyze the structure of the right inverse for a few tournaments.

print("\nAnalyzing right inverse structure...")
print("(Looking at F expressed in A₃ coordinates)")

# For each (T, z), find the MINIMUM NORM filling
# (i.e., the unique filling w with ||w|| minimal, which is the pseudoinverse)

inverse_data = []  # Store filling data for analysis

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2); d3 = dim_om(om3)
    if d2 == 0 or d3 == 0: continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk2
    if z2_dim == 0: continue

    z2_basis_om2 = np.linalg.svd(coords2, full_matrices=True)[2][rk2:].T

    bd3 = build_full_boundary_matrix(ap3, ap2)
    bd3_om = bd3 @ om3
    coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]

    ap2_list = [tuple(x) for x in ap2]
    ap3_list = [tuple(x) for x in ap3]

    # For each Z₂ basis vector, find the min-norm Ω₃ preimage
    for z_idx in range(z2_dim):
        z_om2 = z2_basis_om2[:, z_idx]
        # Solve: coords3 @ x = z_om2 (min norm)
        x, _, _, _ = np.linalg.lstsq(coords3, z_om2, rcond=None)

        # Convert to A₃ coordinates
        filling_a3 = om3 @ x
        cycle_a2 = om2 @ z_om2

        # For each nonzero term in the filling, record the 4-path info
        for i in range(len(ap3_list)):
            if abs(filling_a3[i]) < 1e-10: continue
            path = ap3_list[i]
            a,b,c,d = path
            # Classify this extension
            # Which face of the 4-path IS the target cycle element?
            # ∂₃(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
            # The face -(a,b,c) is a "right extension by d"
            # The face (b,c,d) is a "left extension by a" (removing a)
            inverse_data.append({
                'bits': bits, 'z_idx': z_idx,
                'path4': path, 'coeff': filling_a3[i],
                'a_beats_c': A[a][c], 'b_beats_d': A[b][d],
                'a_beats_d': A[a][d], 'c_beats_a': A[c][a],
                'd_beats_b': A[d][b],
            })

    if bits % 200 == 0 and bits > 0:
        print(f"  ... {bits}/{1 << m}")

# Analyze: what cross-arc patterns appear in fillings?
print(f"\nTotal filling terms: {len(inverse_data)}")

# Cross-arc classification of 4-paths used in fillings
cross_patterns = Counter()
for d in inverse_data:
    pattern = (d['a_beats_c'], d['b_beats_d'], d['a_beats_d'])
    cross_patterns[pattern] += 1

print(f"\nCross-arc patterns (a→c, b→d, a→d) in filling 4-paths:")
for pattern, count in sorted(cross_patterns.items()):
    label = f"{'a→c' if pattern[0] else 'c→a'}, {'b→d' if pattern[1] else 'd→b'}, {'a→d' if pattern[2] else 'd→a'}"
    print(f"  ({label}): {count}")

# KEY: Are DT paths (a→c AND b→d) always sufficient for filling?
dt_count = sum(c for (ac,bd,ad), c in cross_patterns.items() if ac and bd)
non_dt = sum(c for (ac,bd,ad), c in cross_patterns.items() if not (ac and bd))
print(f"\nDT paths in fillings: {dt_count}")
print(f"Non-DT paths in fillings: {non_dt}")

# Now: test if the filling coefficient depends ONLY on the 4-vertex sub-tournament type
print(f"\n{'='*70}")
print("LOCALITY TEST: Does filling depend only on 4-vertex sub-tournament?")
print("=" * 70)

# For each 4-path (a,b,c,d), the coefficient in the filling depends on
# the full tournament T. But does it depend only on T[{a,b,c,d}]?

# Group fillings by (4-vertex sub-tournament type, position in 4-vertex sub-tournament)
from_4sub = defaultdict(list)

for bits in range(min(200, 1 << m)):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2); d3 = dim_om(om3)
    if d2 == 0 or d3 == 0: continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk2
    if z2_dim == 0: continue

    z2_basis_om2 = np.linalg.svd(coords2, full_matrices=True)[2][rk2:].T

    bd3 = build_full_boundary_matrix(ap3, ap2)
    bd3_om = bd3 @ om3
    coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]

    ap2_list = [tuple(x) for x in ap2]
    ap3_list = [tuple(x) for x in ap3]

    for z_idx in range(z2_dim):
        z_om2 = z2_basis_om2[:, z_idx]
        x, _, _, _ = np.linalg.lstsq(coords3, z_om2, rcond=None)
        filling_a3 = om3 @ x
        cycle_a2 = om2 @ z_om2

        for i in range(len(ap3_list)):
            if abs(filling_a3[i]) < 1e-10: continue
            a,b,c,d = ap3_list[i]
            # Get the 4-vertex sub-tournament type (normalized)
            sub = [[0]*4 for _ in range(4)]
            verts = [a,b,c,d]
            for p in range(4):
                for q in range(4):
                    if p != q:
                        sub[p][q] = A[verts[p]][verts[q]]

            # Encode as a tuple
            sub_key = tuple(sub[p][q] for p in range(4) for q in range(4) if p != q)

            # Also record the ratio filling_coeff / cycle_coeff for each face
            face_coeff = {}
            faces = [(1, (b,c,d)), (-1, (a,c,d)), (1, (a,b,d)), (-1, (a,b,c))]
            for sign, face in faces:
                if face in ap2_list:
                    f_idx = ap2_list.index(face)
                    if abs(cycle_a2[f_idx]) > 1e-10:
                        face_coeff[face] = filling_a3[i] / cycle_a2[f_idx]

            from_4sub[(sub_key, z_idx)].append({
                'bits': bits, 'path': (a,b,c,d),
                'fill_coeff': filling_a3[i],
                'face_ratios': face_coeff
            })

# Check: for same sub-tournament type and z_idx, are the ratios consistent?
print(f"\nTotal 4-sub groups: {len(from_4sub)}")
# This is too many groups. Let me instead check a simpler question:
# Is the filling always "uniform" within each 4-vertex subset?

print(f"\n{'='*70}")
print("ALTERNATIVE: RATIO ANALYSIS")
print("=" * 70)

# For each (T, z), and each A₃ path in the filling, compute the ratio
# fill_coeff(a,b,c,d) / cycle_coeff(a,b,c)
# If this ratio is always ±1/(n-3), it's a simple averaging formula.

ratio_dist = Counter()

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2); d3 = dim_om(om3)
    if d2 == 0 or d3 == 0: continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk2
    if z2_dim == 0: continue

    z2_basis_om2 = np.linalg.svd(coords2, full_matrices=True)[2][rk2:].T

    bd3 = build_full_boundary_matrix(ap3, ap2)
    bd3_om = bd3 @ om3
    coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]

    ap2_list = [tuple(x) for x in ap2]
    ap3_list = [tuple(x) for x in ap3]

    for z_idx in range(z2_dim):
        z_om2 = z2_basis_om2[:, z_idx]
        x, _, _, _ = np.linalg.lstsq(coords3, z_om2, rcond=None)
        filling_a3 = om3 @ x
        cycle_a2 = om2 @ z_om2

        for i in range(len(ap3_list)):
            if abs(filling_a3[i]) < 1e-10: continue
            a,b,c,d = ap3_list[i]

            # The face -(a,b,c) has sign -1 in ∂₃
            # So ∂₃(a,b,c,d) contributes -fill_coeff to (a,b,c)
            abc_idx = ap2_list.index((a,b,c)) if (a,b,c) in ap2_list else -1
            if abc_idx >= 0 and abs(cycle_a2[abc_idx]) > 1e-10:
                ratio = -filling_a3[i] / cycle_a2[abc_idx]
                ratio_dist[round(ratio, 4)] += 1

    if bits % 200 == 0 and bits > 0:
        print(f"  ... {bits}/{1 << m}")

print(f"\nRatio -fill_coeff(a,b,c,d) / cycle_coeff(a,b,c):")
for r, count in sorted(ratio_dist.items()):
    print(f"  {r}: {count}")

print("\nDone.")
