#!/usr/bin/env python3
"""
beta2_explicit_filling.py — Find explicit Ω₃ filling for each Z₂ cycle

Strategy: For each tournament T, find a CANONICAL construction that maps
each Z₂ basis element to an Ω₃ element whose boundary is that cycle.

If we can describe this construction purely in terms of tournament structure,
that would give an algebraic proof of β₂=0.

Key question: What is the structure of the filling? Does it decompose
nicely in terms of 4-vertex sub-tournaments?

Author: opus-2026-03-08-S49
"""
import sys
import numpy as np
from collections import Counter, defaultdict
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


def tournament_type(A, verts):
    """Classify a 4-vertex sub-tournament by its score sequence."""
    scores = []
    for v in verts:
        s = sum(A[v][w] for w in verts if w != v)
        scores.append(s)
    return tuple(sorted(scores))


print("=" * 70)
print("EXPLICIT Ω₃ FILLING OF Z₂ CYCLES")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

# Analyze a few representative tournaments
test_bits = [0, 1, 5, 15, 31, 100, 200, 500, 1000]

for bits in test_bits:
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    scores = tuple(sorted(sum(A[i][j] for j in range(n) if j != i) for i in range(n)))

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2)
    d3 = dim_om(om3)

    if d2 == 0:
        continue

    # Compute Z₂ = ker(∂₂) in Ω₂ coords
    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    U2, S2, Vt2 = np.linalg.svd(coords2, full_matrices=True)
    rk2 = int(sum(s > 1e-8 for s in S2))
    z2_dim = d2 - rk2

    if z2_dim == 0:
        continue

    z2_basis = Vt2[rk2:].T  # d2 × z2_dim

    # Compute B₂ = im(∂₃) in Ω₂ coords
    if d3 > 0:
        bd3 = build_full_boundary_matrix(ap3, ap2)
        bd3_om = bd3 @ om3
        coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]
        rk3 = np.linalg.matrix_rank(coords3, tol=1e-8)
    else:
        coords3 = np.zeros((d2, 0))
        rk3 = 0

    beta2 = z2_dim - rk3

    print(f"\nT#{bits} scores={scores}: dim(Ω₂)={d2}, dim(Ω₃)={d3}, "
          f"dim(Z₂)={z2_dim}, rk(∂₃)={rk3}, β₂={beta2}")

    if beta2 != 0:
        print("  WARNING: β₂ ≠ 0!")
        continue

    # For each Z₂ basis element, find which Ω₃ element fills it
    # Solve: coords3 @ x = z2_basis for each column
    ap2_list = [tuple(x) for x in ap2]
    ap3_list = [tuple(x) for x in ap3]

    for z_idx in range(z2_dim):
        z_vec = z2_basis[:, z_idx]
        # Find x such that coords3 @ x = z_vec
        x, residual, _, _ = np.linalg.lstsq(coords3, z_vec, rcond=None)
        reconstruction = coords3 @ x
        err = np.max(np.abs(reconstruction - z_vec))

        if err > 1e-6:
            print(f"  Z₂ basis {z_idx}: CANNOT FILL (error={err:.2e})")
            continue

        # x is the Ω₃ coords of the filling element
        # Convert to A₃ coords
        filling_A3 = om3 @ x

        # Which A₃ paths have nonzero coefficient?
        nonzero = [(ap3_list[i], round(filling_A3[i], 4))
                   for i in range(len(filling_A3)) if abs(filling_A3[i]) > 1e-8]

        # Also express the Z₂ cycle in A₂ coords
        cycle_A2 = om2 @ z_vec
        cycle_terms = [(ap2_list[i], round(cycle_A2[i], 4))
                      for i in range(len(cycle_A2)) if abs(cycle_A2[i]) > 1e-8]

        print(f"  Z₂[{z_idx}] cycle ({len(cycle_terms)} terms):")
        for path, coeff in sorted(cycle_terms):
            print(f"    {coeff:+.4f} * {path}")

        print(f"  Filled by Ω₃ element ({len(nonzero)} terms):")
        for path, coeff in sorted(nonzero):
            # Classify each 4-path
            a,b,c,d = path
            # Check: which arcs exist beyond the path?
            cross = []
            if A[a][c]: cross.append('a→c')
            if A[c][a]: cross.append('c→a')
            if A[b][d]: cross.append('b→d')
            if A[d][b]: cross.append('d→b')
            if A[a][d]: cross.append('a→d')
            if A[d][a]: cross.append('d→a')
            print(f"    {coeff:+.4f} * {path}  cross={cross}")

        # What 4-vertex sub-tournament types appear?
        sub_types = Counter()
        for path, coeff in nonzero:
            verts = list(path)
            sub_types[tournament_type(A, verts)] += 1
        print(f"  Sub-tournament types: {dict(sub_types)}")

print("\n" + "=" * 70)
print("CANONICAL FILLING STRUCTURE ANALYSIS")
print("=" * 70)

# For ALL n=5 tournaments, characterize the filling
filling_stats = Counter()
filling_by_score = defaultdict(list)

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    scores = tuple(sorted(sum(A[i][j] for j in range(n) if j != i) for i in range(n)))

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2)
    d3 = dim_om(om3)
    if d2 == 0:
        continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk2

    if z2_dim == 0:
        continue

    z2_basis = np.linalg.svd(coords2, full_matrices=True)[2][rk2:].T

    if d3 > 0:
        bd3 = build_full_boundary_matrix(ap3, ap2)
        bd3_om = bd3 @ om3
        coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]
        rk3 = np.linalg.matrix_rank(coords3, tol=1e-8)
    else:
        coords3 = np.zeros((d2, 0))
        rk3 = 0

    beta2 = z2_dim - rk3
    ap3_list = [tuple(x) for x in ap3]

    # For each Z₂ basis, find minimum-support filling
    for z_idx in range(z2_dim):
        z_vec = z2_basis[:, z_idx]
        x, _, _, _ = np.linalg.lstsq(coords3, z_vec, rcond=None)
        filling_A3 = om3 @ x

        n_nonzero = sum(abs(filling_A3[i]) > 1e-8 for i in range(len(filling_A3)))
        filling_stats[n_nonzero] += 1

        # Count how many distinct 4-vertex subsets used
        used_quads = set()
        for i in range(len(filling_A3)):
            if abs(filling_A3[i]) > 1e-8:
                path = ap3_list[i]
                used_quads.add(tuple(sorted(path)))

        filling_by_score[scores].append(len(used_quads))

    if bits % 200 == 0 and bits > 0:
        print(f"  ... {bits}/{1 << m}")

print(f"\nFilling support sizes (# nonzero A₃ paths): {dict(sorted(filling_stats.items()))}")
print(f"\n# distinct 4-vertex subsets used in fillings, by score:")
for scores in sorted(filling_by_score.keys()):
    vals = filling_by_score[scores]
    print(f"  {scores}: {Counter(vals)}")

# KEY TEST: Is the filling always decomposable into 4-vertex sub-tournaments?
print("\n" + "=" * 70)
print("4-VERTEX DECOMPOSITION TEST")
print("=" * 70)

# For each 4-vertex sub-tournament S of T, compute Z₂(S) and B₂(S).
# If β₂(S) = 0 for all 4-vertex sub-tournaments (which we know),
# then can we use the sub-tournament fillings to fill the big cycle?

from itertools import combinations

decomp_works = 0
decomp_fails = 0

for bits in range(min(256, 1 << m)):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap0 = enumerate_allowed_paths(A, n, 0)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d2 = dim_om(om2)
    if d2 == 0:
        continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk2
    if z2_dim == 0:
        continue

    # Collect all ∂₃ images from each 4-vertex subset
    # These are "local fillings" — boundaries of paths within each 4-vertex sub-tournament
    ap2_list = [tuple(x) for x in ap2]
    ap2_idx = {p: i for i, p in enumerate(ap2_list)}

    local_boundaries = []  # each is a vector in A₂ coords

    for quad in combinations(range(n), 4):
        quad_set = set(quad)
        # Find all 3-paths within this quad
        for a in quad:
            for b in quad:
                if b == a: continue
                if not A[a][b]: continue
                for c in quad:
                    if c in (a,b): continue
                    if not A[b][c]: continue
                    for d in quad:
                        if d in (a,b,c): continue
                        if not A[c][d]: continue
                        # (a,b,c,d) is a 3-path in quad
                        # Compute its boundary in A₂
                        bd_vec = np.zeros(len(ap2_list))
                        faces = [
                            (1, (b,c,d)),
                            (-1, (a,c,d)),
                            (1, (a,b,d)),
                            (-1, (a,b,c))
                        ]
                        for sign, face in faces:
                            if face in ap2_idx:
                                bd_vec[ap2_idx[face]] += sign
                        local_boundaries.append(bd_vec)

    if not local_boundaries:
        continue

    # Project local boundaries into Ω₂ coords
    local_bd_mat = np.column_stack(local_boundaries)
    local_in_om2 = np.linalg.lstsq(om2, local_bd_mat, rcond=None)[0]

    # Check: do these span Z₂?
    # We need: the projection of local_in_om2 into Z₂ spans Z₂
    z2_basis = np.linalg.svd(coords2, full_matrices=True)[2][rk2:].T
    proj = z2_basis.T @ local_in_om2  # z2_dim × n_local
    rk_proj = np.linalg.matrix_rank(proj, tol=1e-8)

    if rk_proj == z2_dim:
        decomp_works += 1
    else:
        decomp_fails += 1
        if decomp_fails <= 3:
            scores = tuple(sorted(sum(A[i][j] for j in range(n) if j != i) for i in range(n)))
            print(f"  FAIL: bits={bits}, scores={scores}, z2_dim={z2_dim}, rk_proj={rk_proj}")

print(f"\n4-vertex decomposition: works={decomp_works}, fails={decomp_fails}")
if decomp_fails == 0:
    print("✓ All Z₂ cycles can be filled using 4-vertex sub-tournament paths!")
else:
    print("✗ 4-vertex decomposition doesn't always work")

print("\nDone.")
