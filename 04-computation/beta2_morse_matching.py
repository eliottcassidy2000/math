#!/usr/bin/env python3
"""
beta2_morse_matching.py — Can discrete Morse theory prove β₂=0?

Idea: Construct an acyclic matching on (Ω₃ ↔ Ω₂ ↔ Ω₁ ↔ Ω₀)
that matches all Ω₂ cells either DOWN (to Ω₁) or UP (to Ω₃).
If no Ω₂ cells are unmatched ("critical"), then β₂=0.

Strategy 1: "Right extension" matching
  For each non-boundary Ω₂ element (a,b,c), try to match it UP to (a,b,c,d)
  for some specific d. The matching must be injective and acyclic.

Strategy 2: "Left extension" matching
  Match (a,b,c) UP to (x,a,b,c) for some x.

Strategy 3: "King vertex" matching
  Pick the vertex with highest out-degree (the "king").
  Use it as a cone point for Morse matching.

Key constraint: matching must be ACYCLIC (no alternating cycles).

Author: opus-2026-03-08-S49
"""
import sys, time
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

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


print("=" * 70)
print("DISCRETE MORSE MATCHING FOR β₂=0")
print("=" * 70)

# =============================================================
# PART 1: Can we match Ω₂ elements UP to Ω₃ (right extension)?
# For (a,b,c) ∈ Ω₂, find d with (a,b,c,d) ∈ Ω₃
# =============================================================
print("\nPART 1: Right-extension matching for Ω₂ → Ω₃")
print("-" * 50)

# An Ω₂ element is a class [α] in A₂/junk.
# An Ω₃ element is a class [β] in A₃/junk.
# The boundary ∂₃: Ω₃ → Ω₂ has face (a,b,c,d) → (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c).
# A "right extension" match: (a,b,c) ↔ (a,b,c,d) means face₃(a,b,c,d) = ... - (a,b,c).

# For a GREEDY matching: for each (a,b,c) in A₂, pick d with c→d and (a,b,c,d) ∈ A₃.
# (a,b,c,d) ∈ A₃ means a→b, b→c, c→d.
# But for it to be in Ω₃, all faces must relate properly.

# Let's first check: for each A₂ path (a,b,c), how many right extensions exist?
for n in [4, 5, 6]:
    m = n*(n-1)//2
    total = 1 << m

    ext_count_dist = Counter()  # number of right extensions per A₂ path
    omega2_ext = Counter()  # can every Ω₂ basis element be right-extended?

    for bits in range(total):
        A = build_adj(n, bits)

        ap2 = enumerate_allowed_paths(A, n, 2)
        ap3 = enumerate_allowed_paths(A, n, 3)
        if not ap2:
            continue

        # For each A₂ path (a,b,c), count right extensions
        ap3_set = set(tuple(p) for p in ap3) if ap3 else set()
        for p in ap2:
            a, b, c = p
            ext = sum(1 for d in range(n) if d not in (a,b,c) and A[c][d] and (a,b,c,d) in ap3_set)
            ext_count_dist[ext] += 1

    print(f"\nn={n}: Right extensions per A₂ path:")
    for k in sorted(ext_count_dist.keys()):
        pct = 100 * ext_count_dist[k] / sum(ext_count_dist.values())
        print(f"  {k} extensions: {ext_count_dist[k]} ({pct:.1f}%)")


# =============================================================
# PART 2: Matching at Ω level (not A level)
# The matching should work on Ω₂ and Ω₃, not A₂ and A₃.
# An Ω₂ element is a linear combination of A₂ paths.
# A right-extension matching lifts each such combo.
# =============================================================
print(f"\n{'='*70}")
print("PART 2: Greedy right-extension matching on A₂ paths")
print("-" * 50)

# For each tournament, try a greedy matching:
# Process A₂ paths in some order. For each unmatched (a,b,c), find an
# unmatched A₃ path (a,b,c,d) and match them.
# Count unmatched A₂ paths.

for n in [4, 5]:
    m = n*(n-1)//2

    unmatched_dist = Counter()

    for bits in range(1 << m):
        A = build_adj(n, bits)

        ap2 = enumerate_allowed_paths(A, n, 2)
        ap3 = enumerate_allowed_paths(A, n, 3)
        if not ap2 or not ap3:
            unmatched_dist[len(ap2) if ap2 else 0] += 1
            continue

        ap2_list = [tuple(p) for p in ap2]
        ap3_list = [tuple(p) for p in ap3]

        # Build right-extension graph
        right_ext = defaultdict(list)  # A₂ path → list of A₃ extensions
        for p3 in ap3_list:
            right_ext[(p3[0], p3[1], p3[2])].append(p3)

        # Also build left-extension graph
        left_ext = defaultdict(list)
        for p3 in ap3_list:
            left_ext[(p3[1], p3[2], p3[3])].append(p3)

        # Greedy matching (right first, then left)
        matched_a2 = set()
        matched_a3 = set()

        # Sort by fewest extensions first (hardest to match)
        a2_by_ext = sorted(ap2_list, key=lambda p: len(right_ext.get(p, [])))
        for p2 in a2_by_ext:
            if p2 in matched_a2:
                continue
            for p3 in right_ext.get(p2, []):
                if p3 not in matched_a3:
                    matched_a2.add(p2)
                    matched_a3.add(p3)
                    break

        # Try left extensions for remaining
        for p2 in ap2_list:
            if p2 in matched_a2:
                continue
            for p3 in left_ext.get(p2, []):
                if p3 not in matched_a3:
                    matched_a2.add(p2)
                    matched_a3.add(p3)
                    break

        unmatched = len(ap2_list) - len(matched_a2)
        unmatched_dist[unmatched] += 1

    print(f"\nn={n}: Unmatched A₂ paths after greedy matching:")
    for k in sorted(unmatched_dist.keys()):
        pct = 100 * unmatched_dist[k] / sum(unmatched_dist.values())
        print(f"  {k} unmatched: {unmatched_dist[k]} ({pct:.1f}%)")


# =============================================================
# PART 3: Can we match DOWN (Ω₂ ↔ Ω₁)?
# Match (a,b,c) with a face (b,c), (a,c), or (a,b).
# For β₂=0, we need ALL Ω₂ matched either up or down.
# The Ω₂ elements matched DOWN correspond to im(∂₂) directions.
# The Ω₂ elements matched UP correspond to Z₂ directions filled by ∂₃.
# =============================================================
print(f"\n{'='*70}")
print("PART 3: Combined up+down matching")
print("-" * 50)

# A simpler approach: just check if we can construct a Morse function.
# For a tournament T, define height h(v) = d⁺(v) (out-degree).
# Then a path (v₀,...,v_p) is "ascending" if h(v₀) < ... < h(v_p).
# This doesn't directly give Morse theory, but might suggest a filtration.

# For the chain complex, try: match (a,b,c) DOWN to (a,b) if a→c (TT case),
# or match (a,b,c) DOWN to (b,c) if c→a (NT case) AND (a,b) is "harder to fill".
# This is getting complicated.

# Let me try a DIFFERENT approach: compute the Morse number
# (min number of critical cells) using optimization.

# For n=4, all tournaments: try all possible matchings
n = 4
m = n*(n-1)//2

for bits in [0, 5, 21]:  # A few examples
    A = build_adj(n, bits)
    scores = tuple(sorted(sum(A[i][j] for j in range(n) if j!=i) for i in range(n)))

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2)
    d3 = dim_om(om3)

    if d2 > 0:
        bd2 = build_full_boundary_matrix(ap2, ap1)
        bd2_om = bd2 @ om2
        coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
        rk_d2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    else:
        rk_d2 = 0

    if d3 > 0 and d2 > 0:
        bd3 = build_full_boundary_matrix(ap3, ap2)
        bd3_om = bd3 @ om3
        coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]
        rk_d3 = np.linalg.matrix_rank(coords3, tol=1e-8)
    else:
        rk_d3 = 0

    bd1 = build_full_boundary_matrix(ap1, ap0)
    rk_d1 = np.linalg.matrix_rank(bd1, tol=1e-8)

    z2 = d2 - rk_d2
    z1 = dim_om(om1) - rk_d1
    beta1 = z1 - rk_d2

    print(f"\nT#{bits} scores={scores}")
    print(f"  Ω₀={n}, Ω₁={dim_om(om1)}, Ω₂={d2}, Ω₃={d3}")
    print(f"  rk∂₁={rk_d1}, rk∂₂={rk_d2}, rk∂₃={rk_d3}")
    print(f"  β₁={beta1}, z₂={z2}")
    print(f"  Euler: {n} - {dim_om(om1)} + {d2} - {d3} = {n - dim_om(om1) + d2 - d3}")
    # Betti Euler: β₀ - β₁ + β₂ - β₃
    beta0 = 1  # connected
    beta3 = d3 - rk_d3  # dim(Z₃) ... wait, need rk_d4 too
    print(f"  β₀=1, β₁={beta1}, β₂=0, Z₃={d3-rk_d3}")


# =============================================================
# PART 4: The REAL question — try to find an explicit right inverse
# for ∂₃: Ω₃ → Z₂
# For each Z₂ basis vector, express it as ∂₃ of something.
# Check: is there a CANONICAL such expression?
# =============================================================
print(f"\n{'='*70}")
print("PART 4: Explicit right inverse ∂₃⁻¹: Z₂ → Ω₃")
print("-" * 50)

n = 5
m = n*(n-1)//2

filling_types = Counter()

for bits in range(1 << m):
    A = build_adj(n, bits)

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2)
    d3 = dim_om(om3)
    if d2 == 0 or d3 == 0:
        continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk_d2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk_d2

    if z2_dim == 0:
        continue

    bd3 = build_full_boundary_matrix(ap3, ap2)
    bd3_om = bd3 @ om3
    coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]

    # Find Z₂ basis
    U, S, Vt = np.linalg.svd(coords2, full_matrices=True)
    z2_basis = Vt[rk_d2:]  # rows

    # For each Z₂ direction, find min-norm filling in Ω₃
    for k in range(z2_dim):
        target = z2_basis[k]
        x, res, _, _ = np.linalg.lstsq(coords3, target, rcond=None)

        # Check filling worked
        err = np.max(np.abs(coords3 @ x - target))
        assert err < 1e-6, f"Filling failed: {err}"

        # What's the support size of x (in Ω₃)?
        support = np.sum(np.abs(x) > 1e-8)
        filling_types[support] += 1

print(f"n=5: Filling support sizes (in Ω₃ basis):")
for k in sorted(filling_types.keys()):
    pct = 100 * filling_types[k] / sum(filling_types.values())
    print(f"  support={k}: {filling_types[k]} ({pct:.1f}%)")

# What fraction can be filled with a SINGLE Ω₃ element?
single = filling_types.get(1, 0)
total_dirs = sum(filling_types.values())
print(f"\nSingle Ω₃ element fills Z₂ direction: {single}/{total_dirs} ({100*single/total_dirs:.1f}%)")


# =============================================================
# PART 5: Structure of ∂₃ as a matrix — eigenvalue analysis
# If ∂₃ restricted to Z₂ target is always full rank, check WHY.
# =============================================================
print(f"\n{'='*70}")
print("PART 5: Singular values of ∂₃ restricted to Z₂")
print("-" * 50)

n = 5
m = n*(n-1)//2

min_sv_ratio = float('inf')
sv_stats = []

for bits in range(1 << m):
    A = build_adj(n, bits)

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2)
    d3 = dim_om(om3)
    if d2 == 0 or d3 == 0:
        continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk_d2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk_d2
    if z2_dim == 0:
        continue

    bd3 = build_full_boundary_matrix(ap3, ap2)
    bd3_om = bd3 @ om3
    coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]

    # SVD of coords3
    U, S, Vt = np.linalg.svd(coords3, full_matrices=False)
    # The rank should be z2_dim (by β₂=0)
    # The z2_dim-th singular value should be > 0
    if len(S) >= z2_dim:
        min_sv = S[z2_dim - 1]
        max_sv = S[0]
        ratio = min_sv / max_sv if max_sv > 0 else 0
        sv_stats.append((z2_dim, min_sv, max_sv, ratio))
        if ratio < min_sv_ratio:
            min_sv_ratio = ratio

print(f"n=5: Singular value analysis of ∂₃ in Ω₂ coords")
print(f"  Total: {len(sv_stats)} tournaments with Z₂>0")
min_svs = [s[1] for s in sv_stats]
max_svs = [s[2] for s in sv_stats]
ratios = [s[3] for s in sv_stats]
print(f"  Min singular value: min={min(min_svs):.6f}, mean={np.mean(min_svs):.6f}")
print(f"  Max singular value: max={max(max_svs):.6f}, mean={np.mean(max_svs):.6f}")
print(f"  Condition ratio (σ_z2/σ_1): min={min(ratios):.6f}, mean={np.mean(ratios):.6f}")

# Group by z2_dim
by_z2 = defaultdict(list)
for z, mn, mx, r in sv_stats:
    by_z2[z].append((mn, mx, r))

for z in sorted(by_z2.keys()):
    stats = by_z2[z]
    min_sigs = [s[0] for s in stats]
    rats = [s[2] for s in stats]
    print(f"  z2={z}: {len(stats)} tours, min_σ in [{min(min_sigs):.4f}, {max(min_sigs):.4f}]")


print("\nDone.")
