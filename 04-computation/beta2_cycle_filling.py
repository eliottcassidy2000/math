#!/usr/bin/env python3
"""Analyze HOW each 2-cycle is filled by ∂_3 — explicit filling construction.

Key structural insight: each 2-path (a,b,c) has AT MOST one bad face (a,c),
occurring iff c→a. So the Ω_2 junk constraints are on DISJOINT variable sets.

Question: can we construct an EXPLICIT filling w ∈ Ω_3 for any z ∈ ker(∂_2) ∩ Ω_2?

Approach: for each tournament, find the actual filling map (right inverse of ∂_3).
Look for patterns in how the filling works.
"""
import numpy as np
from itertools import combinations
import sys
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
)

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

# ===== Part 1: Explicit cycle structure at n=5 =====
print("=" * 70)
print("EXPLICIT 2-CYCLE STRUCTURE AND FILLING")
print("=" * 70)

n = 5
# Pick a specific interesting tournament
# Use the regular tournament (2,2,2,2,2)
A_reg = [[0]*n for _ in range(n)]
# C_5: 0→1→2→3→4→0 and 0→2→4→1→3→0 (doubly regular)
for i in range(n):
    A_reg[i][(i+1)%n] = 1
    A_reg[i][(i+2)%n] = 1

print(f"\nRegular tournament C_5^{{1,2}}:")
for i in range(n):
    print(f"  {A_reg[i]}")

a = {}
for p in range(5):
    a[p] = [tuple(x) for x in enumerate_allowed_paths(A_reg, n, p)]

om = {}
for p in range(5):
    if p == 0:
        om[p] = np.eye(n)
    else:
        om[p] = compute_omega_basis(A_reg, n, p, a[p], a[p-1])

dim_om2 = om[2].shape[1]
dim_om3 = om[3].shape[1]
print(f"  dim(Ω_2)={dim_om2}, dim(Ω_3)={dim_om3}")

# Compute boundary maps in Ω coordinates
bd2 = build_full_boundary_matrix(a[2], a[1])
bd2_om = bd2 @ om[2]

# Express in Ω_1 coords (Ω_1 = A_1 for tournaments)
bd2_coords = bd2_om  # Already in A_1 = Ω_1 coords

U2, S2, Vt2 = np.linalg.svd(bd2_coords)
rank_d2 = sum(s > 1e-8 for s in S2)
z2_dim = dim_om2 - rank_d2
print(f"  rank(∂_2)={rank_d2}, dim(Z_2)={z2_dim}")

# Get Z_2 basis (kernel of ∂_2 in Ω_2 coords)
ker2_basis = Vt2[rank_d2:].T  # dim_om2 × z2_dim
print(f"  Z_2 basis shape: {ker2_basis.shape}")

# Express each Z_2 basis vector in A_2 coordinates
print(f"\n  Z_2 basis vectors (in A_2 coordinates):")
for j in range(z2_dim):
    cycle_om = ker2_basis[:, j]
    cycle_a2 = om[2] @ cycle_om
    nonzero = [(a[2][i], cycle_a2[i]) for i in range(len(cycle_a2)) if abs(cycle_a2[i]) > 1e-10]
    paths_str = ", ".join(f"{c:.3f}*{p}" for p, c in nonzero)
    print(f"    z_{j}: {paths_str}")

# Compute ∂_3 in Ω coords
bd3 = build_full_boundary_matrix(a[3], a[2])
bd3_om = bd3 @ om[3]
# Express in Ω_2 coords
im3_coords, _, _, _ = np.linalg.lstsq(om[2], bd3_om, rcond=None)
print(f"\n  ∂_3 map (Ω_3 → Ω_2): {im3_coords.shape}")

# Find filling: for each Z_2 basis vector, find w ∈ Ω_3 with ∂_3(w) = z
print(f"\n  Filling each 2-cycle:")
for j in range(z2_dim):
    cycle_om = ker2_basis[:, j]
    # Solve: im3_coords @ w = cycle_om (least norm solution)
    w, residual, rank, sv = np.linalg.lstsq(im3_coords, cycle_om, rcond=None)
    # Check: ∂_3(w) = z?
    check = im3_coords @ w
    error = np.linalg.norm(check - cycle_om)

    # Express w in A_3 coords
    w_a3 = om[3] @ w
    nonzero_w = [(a[3][i], w_a3[i]) for i in range(len(w_a3)) if abs(w_a3[i]) > 1e-10]

    print(f"    z_{j}: filled by {len(nonzero_w)} 3-paths (error={error:.2e})")
    for p, c in nonzero_w[:8]:
        # Check if DT
        is_dt = A_reg[p[0]][p[2]] and A_reg[p[1]][p[3]]
        print(f"      {c:.3f} * {p} (DT={is_dt})")
    if len(nonzero_w) > 8:
        print(f"      ... and {len(nonzero_w)-8} more")

# ===== Part 2: Bad face analysis for Ω_2 and Ω_3 =====
print(f"\n\n{'='*70}")
print("BAD FACE STRUCTURE")
print("=" * 70)

# For Ω_2: bad face (a,c) where c→a
# For Ω_3: bad faces (a,c,d) at pos 1 [c→a] and (a,b,d) at pos 2 [d→b]

print(f"\nAt n=5, analyzing bad face structure for ALL tournaments:")

for tidx, A in enumerate(all_tournaments_gen(n)):
    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    if scores != (2,2,2,2,2):
        continue  # Just look at regular tournament for now

    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]

    # Bad faces of 3-paths
    bad_face_groups = defaultdict(list)  # bad_face -> [(3-path, sign)]
    for path in a3:
        a_, b_, c_, d_ = path
        # Face at position 1: (a,c,d), bad if c→a
        if A[c_][a_]:
            bad_face_groups[(a_, c_, d_)].append((path, -1))
        # Face at position 2: (a,b,d), bad if d→b
        if A[d_][b_]:
            bad_face_groups[(a_, b_, d_)].append((path, +1))

    print(f"\n  Tournament with scores {scores}:")
    print(f"  |A_3|={len(a3)}, bad face groups: {len(bad_face_groups)}")

    # Classify bad face groups
    type_counter = Counter()
    shared_paths = 0

    # Check: does any 3-path appear in TWO bad face groups?
    path_to_groups = defaultdict(list)
    for face, members in bad_face_groups.items():
        for path, sign in members:
            path_to_groups[path].append((face, sign))

    multi_group = {p: g for p, g in path_to_groups.items() if len(g) > 1}
    print(f"  3-paths in multiple bad face groups: {len(multi_group)}/{len(a3)}")
    for path, groups in list(multi_group.items())[:5]:
        signs = [(f, s) for f, s in groups]
        print(f"    {path}: faces {signs}")

    # DT paths: those with NO bad faces
    dt_count = sum(1 for p in a3 if p not in path_to_groups)
    print(f"  DT paths (no bad faces): {dt_count}")

    # 1-bad-face paths
    one_bad = sum(1 for p in a3 if len(path_to_groups.get(p, [])) == 1)
    print(f"  1-bad-face paths: {one_bad}")

    # 2-bad-face paths
    two_bad = len(multi_group)
    print(f"  2-bad-face paths: {two_bad}")

    # Bad face group sizes
    group_sizes = Counter(len(members) for members in bad_face_groups.values())
    print(f"  Bad face group sizes: {dict(sorted(group_sizes.items()))}")

    break  # Just one tournament for now

# ===== Part 3: Verify the disjointness at Ω_2 level =====
print(f"\n\n{'='*70}")
print("VERIFYING Ω_2 DISJOINTNESS AND FORMULA")
print("=" * 70)

n = 5
all_correct = True
for A in all_tournaments_gen(n):
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]

    # Count transitive triples
    tt = sum(1 for (a_,b_,c_) in a2 if A[a_][c_])

    # Count bad pairs and their multiplicities
    bad_pair_mult = defaultdict(int)
    for (a_,b_,c_) in a2:
        if A[c_][a_]:  # c→a, so (a,c) is bad
            bad_pair_mult[(a_,c_)] += 1

    cancellation = sum(max(0, m-1) for m in bad_pair_mult.values())

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0

    predicted = tt + cancellation
    if predicted != dim_om2:
        print(f"  MISMATCH: tt={tt}, cancel={cancellation}, predicted={predicted}, actual={dim_om2}")
        all_correct = False

print(f"  dim(Ω_2) = TT + Σ(mult-1) for all n=5 tournaments: {'CONFIRMED' if all_correct else 'FAILED'}")

# ===== Part 4: Can we decompose ker(∂_2) into "local" pieces? =====
print(f"\n\n{'='*70}")
print("LOCAL DECOMPOSITION OF ker(∂_2)")
print("=" * 70)

# Question: can ker(∂_2) be decomposed into contributions from individual
# bad pairs? Each bad pair (a,c) with c→a and intermediaries b_1,...,b_k
# generates a (k-1)-dimensional subspace of Ω_2. The boundary ∂_2 of these
# elements involves edges incident to a and c. Are the ker(∂_2) constraints
# "local" to each bad pair?

n = 5
for tidx, A in enumerate(all_tournaments_gen(n)):
    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    t3 = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if A[i][j] and A[j][k] and A[k][i]: t3 += 1
                if A[j][i] and A[i][k] and A[k][j]: t3 += 1

    if t3 != 4 or scores != (1,2,2,2,3):
        continue

    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    # Decompose A_2 by bad pair
    tt_paths = []
    bad_pair_paths = defaultdict(list)

    for idx, (a_, b_, c_) in enumerate(a2):
        if A[a_][c_]:
            tt_paths.append(idx)
        else:
            bad_pair_paths[(a_, c_)].append(idx)

    print(f"\n  Tournament: scores={scores}, t3={t3}")
    print(f"  |A_2|={len(a2)}, TT={len(tt_paths)}, bad_pairs={len(bad_pair_paths)}")
    print(f"  dim(Ω_2)={dim_om2}, dim(Ω_3)={dim_om3}")

    # Ω_2 basis: TT paths are independent Ω_2 elements.
    # Bad pair (a,c) with mults b_1,...,b_k: contributes k-1 elements
    # of the form (a,b_i,c) - (a,b_j,c) to Ω_2.
    print(f"\n  Bad pairs and multiplicities:")
    for (a_, c_), indices in sorted(bad_pair_paths.items()):
        paths = [a2[i] for i in indices]
        mids = [p[1] for p in paths]
        print(f"    ({a_},{c_}) c→a: intermediaries {mids} (mult={len(mids)})")

    # Now: what does ∂_2 do to each TT and cancellation element?
    bd2 = build_full_boundary_matrix(a2, a1)

    print(f"\n  ∂_2 on TT paths:")
    for idx in tt_paths[:5]:
        path = a2[idx]
        bd = bd2[idx, :]
        nonzero = [(a1[j], bd[j]) for j in range(len(a1)) if abs(bd[j]) > 1e-10]
        print(f"    ∂({path}) = {' + '.join(f'{c:+.0f}*{p}' for p,c in nonzero)}")

    print(f"\n  ∂_2 on cancellation elements (a,b_i,c)-(a,b_j,c):")
    for (a_, c_), indices in sorted(bad_pair_paths.items()):
        if len(indices) < 2:
            continue
        # Take (a,b_1,c) - (a,b_2,c)
        i1, i2 = indices[0], indices[1]
        p1, p2 = a2[i1], a2[i2]
        bd_diff = bd2[i1, :] - bd2[i2, :]
        nonzero = [(a1[j], bd_diff[j]) for j in range(len(a1)) if abs(bd_diff[j]) > 1e-10]
        print(f"    ∂({p1}-{p2}) = {' + '.join(f'{c:+.0f}*{p}' for p,c in nonzero)}")

    # Key observation: ∂(a,b_1,c) - ∂(a,b_2,c) = (b_1,c)-(b_2,c) - (a,c)+(a,c) + (a,b_1)-(a,b_2)
    #   = (b_1,c) - (b_2,c) + (a,b_1) - (a,b_2)
    # The (a,c) terms CANCEL! So the boundary involves only edges incident to b_1 and b_2.
    print(f"\n  NOTE: ∂((a,b1,c)-(a,b2,c)) = (b1,c)-(b2,c)+(a,b1)-(a,b2)")
    print(f"  The (a,c) terms cancel. Boundary only involves edges at b1, b2.")

    break

# ===== Part 5: Universal pattern: ∂_2 on cancellation elements =====
print(f"\n\n{'='*70}")
print("UNIVERSAL BOUNDARY STRUCTURE")
print("=" * 70)

# ∂(a,b,c) = (b,c) - (a,c) + (a,b)
# For c→a: (a,c) is bad, but in Ω_2 the constraint is Σ_b z(a,b,c) = 0.
# Cancellation basis: (a,b_i,c) - (a,b_{i+1},c) for i=1,...,k-1
# ∂ of cancellation: (b_i,c) - (b_{i+1},c) + (a,b_i) - (a,b_{i+1})
#
# TT basis: (a,b,c) with a→c
# ∂ of TT: (b,c) - (a,c) + (a,b) where (a,c) ∈ A_1 since a→c
#
# The boundary ∂_2 has two types of terms:
# 1. TT contributions: each TT path gives 3 edge terms
# 2. Cancellation contributions: each gives 4 edge terms (2 pairs)
#
# For z = Σ α_i TT_i + Σ β_j cancel_j ∈ ker(∂_2):
# The boundary = 0 gives constraints on α_i and β_j.

print("Key structural observation for tournament β_2=0:")
print()
print("For tournament T on n vertices:")
print("  1. Ω_2 has a NATURAL basis:")
print("     - Transitive triples (a,b,c): a→b→c, a→c")
print("     - Cancellation pairs (a,b1,c)-(a,b2,c): c→a, a→b1→c, a→b2→c")
print()
print("  2. The boundary ∂_2 on cancellation pairs is SIMPLE:")
print("     ∂((a,b1,c)-(a,b2,c)) = (b1,c)-(b2,c)+(a,b1)-(a,b2)")
print("     No (a,c) terms! Only edges at b1, b2.")
print()
print("  3. A 2-cycle z ∈ ker(∂_2)∩Ω_2 must have:")
print("     For each edge (x,y): net coefficient = 0")
print("     TT and cancellation terms contribute independently")
print()

# ===== Part 6: Explicit cycle-filling construction attempt =====
print(f"\n{'='*70}")
print("CYCLE-FILLING CONSTRUCTION")
print("=" * 70)

# Try to find a universal filling rule:
# For a 2-cycle z = Σ z_{abc} (a,b,c), find w = Σ w_{abcd} (a,b,c,d) with ∂w = z.
#
# Idea: for each term z_{abc}(a,b,c), "extend" it to (v,a,b,c) or (a,b,c,v)
# by choosing an appropriate vertex v.
#
# Extension to the right: (a,b,c,v) where c→v
# Extension to the left: (v,a,b,c) where v→a
# The DT condition requires: v→b (for right extension) or v→b (for left extension)

n = 5
# Use the regular tournament
A = A_reg

a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]
a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]

om2 = compute_omega_basis(A, n, 2, a2, a1)
om3 = compute_omega_basis(A, n, 3, a3, a2)

bd2 = build_full_boundary_matrix(a2, a1)
bd2_om = bd2 @ om2
U2, S2, Vt2 = np.linalg.svd(bd2_om)
rank_d2 = sum(s > 1e-8 for s in S2)
ker2_basis = Vt2[rank_d2:].T

bd3 = build_full_boundary_matrix(a3, a2)
bd3_om = bd3 @ om3
im3_coords, _, _, _ = np.linalg.lstsq(om2, bd3_om, rcond=None)

print(f"\nRegular tournament C_5^{{1,2}}:")
print(f"  Z_2 dim: {ker2_basis.shape[1]}")

# For each 2-cycle, try the "right extension" construction:
# For z = Σ z_i (a_i, b_i, c_i), try w = Σ z_i * Σ_{v: c_i→v, (a,b,c,v)∈A_3} α_{v} (a_i,b_i,c_i,v)

for j in range(ker2_basis.shape[1]):
    cycle_om = ker2_basis[:, j]
    cycle_a2 = om2 @ cycle_om
    nonzero_paths = [(i, a2[i], cycle_a2[i]) for i in range(len(a2)) if abs(cycle_a2[i]) > 1e-10]

    print(f"\n  Cycle z_{j}:")
    for idx, path, coeff in nonzero_paths:
        a_, b_, c_ = path
        # What 3-paths extend this to the right?
        right_ext = [(a_,b_,c_,v) for v in range(n)
                     if v not in {a_,b_,c_} and A[c_][v] and (a_,b_,c_,v) in a3]
        # What 3-paths extend to the left?
        left_ext = [(v,a_,b_,c_) for v in range(n)
                    if v not in {a_,b_,c_} and A[v][a_] and (v,a_,b_,c_) in a3]

        # Are these DT?
        right_dt = [(a_,b_,c_,v) for v in range(n)
                    if v not in {a_,b_,c_} and A[c_][v] and A[a_][c_] and A[b_][v]
                    and (a_,b_,c_,v) in a3]
        left_dt = [(v,a_,b_,c_) for v in range(n)
                   if v not in {a_,b_,c_} and A[v][a_] and A[v][b_] and A[a_][c_]
                   and (v,a_,b_,c_) in a3]

        print(f"    {coeff:+.3f} * {path}: right_ext={len(right_ext)}, left_ext={len(left_ext)}, "
              f"right_DT={len(right_dt)}, left_DT={len(left_dt)}")

print("\nDone.")
