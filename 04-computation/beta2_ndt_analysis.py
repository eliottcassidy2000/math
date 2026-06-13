#!/usr/bin/env python3
"""
beta2_ndt_analysis.py — Analyze non-DT Ω₃ elements at n=6

At n=6, DT boundaries alone span Z₂ in 31808/32768 tournaments.
In 960 cases (all score [1,2,2,3,3,4]), DT rank = 9 but Z₂ = 10.
Full Ω₃ still gives B₂ = Z₂ = 10 (β₂=0).

Question: What do the non-DT Ω₃ elements look like?
A 4-path (a,b,c,d) is in Ω₃ but NOT DT when a→/c or b→/d.
These paths cancel in ∂ despite NOT having the simplicial structure.

Author: opus-2026-03-08-S45
"""
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved


n = 6
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

print("=" * 70)
print("NON-DT Ω₃ ANALYSIS AT n=6")
print("=" * 70)

# Examine one failure case in detail
bits = 154
A = [[0]*n for _ in range(n)]
for idx, (i,j) in enumerate(pairs):
    if (bits >> idx) & 1: A[i][j] = 1
    else: A[j][i] = 1

print(f"\nT#{bits} adjacency:")
for i in range(n):
    out = [j for j in range(n) if A[i][j]]
    print(f"  {i} → {out}")

scores = [sum(A[i]) for i in range(n)]
print(f"Scores: {scores}")

ap0 = enumerate_allowed_paths(A, n, 0)
ap1 = enumerate_allowed_paths(A, n, 1)
ap2 = enumerate_allowed_paths(A, n, 2)
ap3 = enumerate_allowed_paths(A, n, 3)

om1 = compute_omega_basis(A, n, 1, ap1, ap0)
om2 = compute_omega_basis(A, n, 2, ap2, ap1)
om3 = compute_omega_basis(A, n, 3, ap3, ap2)

d2 = om2.shape[1] if om2.ndim == 2 else 0
d3 = om3.shape[1] if om3.ndim == 2 else 0

print(f"\n|A₃| = {len(ap3)}, dim Ω₃ = {d3}")

# Classify all A₃ paths
ap3_tuples = [tuple(p) for p in ap3]
dt_paths = []
ndt_paths = []
for i, (a,b,c,d) in enumerate(ap3_tuples):
    is_dt = A[a][c] and A[b][d]
    if is_dt:
        dt_paths.append(i)
    else:
        ndt_paths.append(i)
        # What condition fails?
        ac = A[a][c]
        bd = A[b][d]
        print(f"  Non-DT path ({a},{b},{c},{d}): a→c={ac}, b→d={bd}")

print(f"\n#DT = {len(dt_paths)}, #non-DT = {len(ndt_paths)}")

# Show Ω₃ basis elements
print(f"\nΩ₃ basis elements:")
for col in range(d3):
    v = om3[:, col]
    terms = [(ap3_tuples[i], round(v[i],4)) for i in range(len(v)) if abs(v[i]) > 1e-8]
    dt_count = sum(1 for p, _ in terms if A[p[0]][p[2]] and A[p[1]][p[3]])
    ndt_count = len(terms) - dt_count
    print(f"  e{col}: {len(terms)} terms ({dt_count} DT, {ndt_count} non-DT)")
    for path, coeff in terms:
        a,b,c,d = path
        ac = A[a][c]
        bd = A[b][d]
        tag = "DT" if ac and bd else f"nDT(ac={ac},bd={bd})"
        print(f"    {coeff:+.4f} * ({a},{b},{c},{d}) [{tag}]")

# Z₂ computation
bd2 = build_full_boundary_matrix(ap2, ap1)
bd2_om = bd2 @ om2
coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
z2 = d2 - rk2
print(f"\nZ₂ = {z2}, dim Ω₂ = {d2}, rk(∂₂) = {rk2}")

# DT boundary rank
bd3 = build_full_boundary_matrix(ap3, ap2)
dt_boundaries = bd3[:, dt_paths]
dt_in_om2 = np.linalg.lstsq(om2, dt_boundaries, rcond=None)[0]
rk_dt = np.linalg.matrix_rank(dt_in_om2, tol=1e-8)
print(f"DT boundary rank in Ω₂: {rk_dt}")

# Full Ω₃ boundary rank
bd3_om = bd3 @ om3
bd3_in_om2 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]
rk_full = np.linalg.matrix_rank(bd3_in_om2, tol=1e-8)
print(f"Full Ω₃ boundary rank: {rk_full}")
print(f"Non-DT contribution: {rk_full - rk_dt} extra dimension(s)")

# What does the missing direction look like?
# Get the kernel direction in Ω₂ that DT misses
print(f"\n{'='*70}")
print("MISSING DIRECTION ANALYSIS")
print("=" * 70)

# DT image in Ω₂
U_dt, S_dt, Vt_dt = np.linalg.svd(dt_in_om2, full_matrices=True)
rk = int(sum(s > 1e-8 for s in S_dt))
# Complement of DT image in Z₂ space

# First get Z₂ basis
U2, S2, Vt2 = np.linalg.svd(coords2, full_matrices=True)
rk2 = int(sum(s > 1e-8 for s in S2))
z2_basis = Vt2[rk2:]  # rows are Z₂ basis vectors in Ω₂ coords

print(f"Z₂ basis: {z2_basis.shape[0]} vectors")

# Project DT image into Z₂
# DT image vectors in Ω₂ coords
dt_image_om2 = bd3_in_om2  # This is the full Ω₃ image... need DT only
dt_image_z2 = z2_basis @ dt_in_om2  # Z₂ coords of DT image
rk_dt_z2 = np.linalg.matrix_rank(dt_image_z2, tol=1e-8)
print(f"DT rank within Z₂: {rk_dt_z2}")

# Full Ω₃ image in Z₂
full_image_z2 = z2_basis @ bd3_in_om2
rk_full_z2 = np.linalg.matrix_rank(full_image_z2, tol=1e-8)
print(f"Full Ω₃ rank within Z₂: {rk_full_z2}")

# Find the missing vector
# SVD of dt_image_z2 to find null complement
U3, S3, Vt3 = np.linalg.svd(dt_image_z2.T, full_matrices=True)
rk3 = int(sum(s > 1e-8 for s in S3))
if rk3 < dt_image_z2.shape[0]:
    missing_z2 = U3[:, rk3:]  # in Z₂ coords
    missing_om2 = z2_basis.T @ missing_z2  # back to Ω₂ coords
    # Convert to A₂ coords
    missing_a2 = om2 @ missing_om2
    print(f"\nMissing direction (in A₂ coords):")
    for col in range(missing_a2.shape[1]):
        v = missing_a2[:, col]
        terms = [(tuple(ap2[i]), round(v[i],4)) for i in range(len(v)) if abs(v[i]) > 1e-8]
        print(f"  Missing Z₂ vector:")
        for path, coeff in terms:
            a,b,c = path
            tt = A[a][c]
            print(f"    {coeff:+.4f} * ({a},{b},{c}) [{'TT' if tt else 'NT'}]")

# Check: are ALL 960 failures isomorphic? (same score seq suggests yes)
print(f"\n{'='*70}")
print("STRUCTURE OF ALL 960 FAILURE CASES")
print("=" * 70)

fail_3cycles = Counter()
fail_scores = Counter()
for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    ap3_l = enumerate_allowed_paths(A, n, 3)
    if not ap3_l:
        continue
    ap3_t = [tuple(p) for p in ap3_l]
    dt_idx = [i for i, (a,b,c,d) in enumerate(ap3_t) if A[a][c] and A[b][d]]

    # Quick check: is this a failure?
    ap2_l = enumerate_allowed_paths(A, n, 2)
    ap1_l = enumerate_allowed_paths(A, n, 1)
    ap0_l = enumerate_allowed_paths(A, n, 0)
    if not ap2_l:
        continue
    om1_l = compute_omega_basis(A, n, 1, ap1_l, ap0_l)
    om2_l = compute_omega_basis(A, n, 2, ap2_l, ap1_l)
    d2_l = om2_l.shape[1] if om2_l.ndim == 2 and om2_l.shape[0] > 0 else 0
    if d2_l == 0:
        continue

    bd2_l = build_full_boundary_matrix(ap2_l, ap1_l)
    bd2_om_l = bd2_l @ om2_l
    coords2_l = np.linalg.lstsq(om1_l, bd2_om_l, rcond=None)[0]
    rk2_l = np.linalg.matrix_rank(coords2_l, tol=1e-8)
    z2_l = d2_l - rk2_l

    if z2_l == 0 or not dt_idx:
        continue

    bd3_l = build_full_boundary_matrix(ap3_l, ap2_l)
    dt_bd = bd3_l[:, dt_idx]
    dt_om2 = np.linalg.lstsq(om2_l, dt_bd, rcond=None)[0]
    rk_dt_l = np.linalg.matrix_rank(dt_om2, tol=1e-8)

    if rk_dt_l < z2_l:
        c3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                 if (A[i][j]+A[j][i])*(A[j][k]+A[k][j])*(A[i][k]+A[k][i]) > 0
                 and ((A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k])))
        fail_3cycles[c3] += 1
        scores = tuple(sorted(sum(A[i]) for i in range(n)))
        fail_scores[scores] += 1

print(f"Failure score distributions: {dict(fail_scores)}")
print(f"Failure 3-cycle counts: {dict(fail_3cycles)}")

# Count ALL 4-path types in failures
print(f"\n{'='*70}")
print("4-PATH TYPES IN Ω₃ FOR FAILURES")
print("=" * 70)

# Detailed look at first failure
bits = 154
A = [[0]*n for _ in range(n)]
for idx, (i,j) in enumerate(pairs):
    if (bits >> idx) & 1: A[i][j] = 1
    else: A[j][i] = 1

ap3 = enumerate_allowed_paths(A, n, 3)
ap3_t = [tuple(p) for p in ap3]
print(f"\nT#{bits}: {len(ap3)} allowed 4-paths")

# Classify more finely
for i, (a,b,c,d) in enumerate(ap3_t):
    ac = A[a][c]
    bd = A[b][d]
    ad = A[a][d]
    # Inner edges: a→b (always), b→c (always), c→d (always)
    # Crossing edges: a→c, b→d, a→d
    print(f"  ({a},{b},{c},{d}): a→c={ac} b→d={bd} a→d={ad} "
          f"{'DT+' if ac and bd and ad else 'DT-' if ac and bd else 'nDT'}")

print("\nDone.")
