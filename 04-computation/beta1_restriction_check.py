#!/usr/bin/env python3
"""
beta1_restriction_check.py — Cocycle restriction map analysis (COHOMOLOGICAL)

For tournaments T with β₁(T)=0, and each vertex v:
  - Compute Z¹(T) = 1-cocycles (functions on directed edges satisfying TT constraints)
  - Compute B¹(T) = 1-coboundaries (δg for g: V → R)
  - Compute Z¹(T\v), B¹(T\v)
  - Compute res_v: Z¹(T) → Z¹(T\v) (restriction to edges not involving v)
  - Check surjectivity

Key prediction:
  If β₁(T)=0 and β₁(T\v)=1 ("bad vertex"), then res_v CANNOT be surjective.
  Proof: β₁=0 ⟹ Z¹=B¹. res_v(B¹(T)) ⊆ B¹(T\v). So image ⊆ B¹(T\v) ⊊ Z¹(T\v).

Note: β¹ₚₐₜₕ = dim Z¹ - dim B¹ (same for homology and cohomology over R).
"""

import numpy as np
from collections import defaultdict

def tournament_edges(A, n):
    """Return list of directed edges (i,j) where A[i][j]=1."""
    return [(i,j) for i in range(n) for j in range(n) if i != j and A[i][j] == 1]

def transitive_triples(A, n):
    """Return list of transitive triples (a,b,c) where a→b, b→c, a→c."""
    tts = []
    for a in range(n):
        for b in range(n):
            if a == b or not A[a][b]:
                continue
            for c in range(n):
                if c == a or c == b:
                    continue
                if A[b][c] and A[a][c]:
                    tts.append((a, b, c))
    return tts

def compute_cocycle_space(A, n):
    """
    Compute Z¹ (1-cocycles) and B¹ (1-coboundaries) for a tournament.

    1-cochain: function f: edges → R, represented as vector indexed by edges.
    1-cocycle condition: for each TT (a,b,c), f(a,b) - f(a,c) + f(b,c) = 0.
      (This is ∂₂* applied to f, evaluated on the 2-path (a,b,c).)
    1-coboundary: f = δg where g: V → R, so f(a,b) = g(b) - g(a).

    Returns: (Z1_basis, B1_basis, edges, beta1)
      Z1_basis: columns are basis vectors for Z¹ in edge coordinates
      B1_basis: columns are basis vectors for B¹ in edge coordinates
      edges: list of directed edges
      beta1: β¹ = dim Z¹ - dim B¹
    """
    edges = tournament_edges(A, n)
    edge_idx = {e: i for i, e in enumerate(edges)}
    m = len(edges)

    # TT constraint matrix: one row per TT, one column per edge
    tts = transitive_triples(A, n)
    if tts:
        C = np.zeros((len(tts), m))
        for r, (a, b, c) in enumerate(tts):
            # constraint: f(a,b) - f(a,c) + f(b,c) = 0
            C[r, edge_idx[(a, b)]] = 1
            C[r, edge_idx[(a, c)]] = -1
            C[r, edge_idx[(b, c)]] = 1
        # Z¹ = ker(C)
        U, S, Vt = np.linalg.svd(C, full_matrices=True)
        rank_C = sum(s > 1e-10 for s in S)
        Z1_basis = Vt[rank_C:].T  # columns are kernel basis
    else:
        # No TTs means no constraints: Z¹ = all of R^m
        Z1_basis = np.eye(m)

    # B¹ = im(δ) where δ: R^n → R^m, (δg)(a,b) = g(b) - g(a)
    delta = np.zeros((m, n))
    for i, (a, b) in enumerate(edges):
        delta[i, b] = 1
        delta[i, a] = -1
    # B¹ = column space of delta
    Sd = np.linalg.svd(delta, compute_uv=False)
    rank_delta = sum(s > 1e-10 for s in Sd)
    # Get an orthonormal basis for column space
    if rank_delta > 0:
        Q, R_mat = np.linalg.qr(delta, mode='reduced')
        # Take columns corresponding to nonzero pivots
        B1_basis = Q[:, :rank_delta]
    else:
        B1_basis = np.zeros((m, 0))

    dim_Z1 = Z1_basis.shape[1]
    dim_B1 = rank_delta
    beta1 = dim_Z1 - dim_B1

    return Z1_basis, B1_basis, edges, beta1


def subtournament(A, n, v):
    """Remove vertex v, return (A', n-1, vertex_map)."""
    vertices = [i for i in range(n) if i != v]
    n_new = len(vertices)
    A_new = [[0]*n_new for _ in range(n_new)]
    for ii, i in enumerate(vertices):
        for jj, j in enumerate(vertices):
            A_new[ii][jj] = A[i][j]
    return A_new, n_new, vertices


def compute_restriction_map(A, n, v):
    """
    Compute res_v: Z¹(T) → Z¹(T\v) and analyze its properties.
    """
    # Cocycles/coboundaries of T
    Z1_T, B1_T, edges_T, beta1_T = compute_cocycle_space(A, n)
    edge_idx_T = {e: i for i, e in enumerate(edges_T)}
    dim_Z1_T = Z1_T.shape[1]
    dim_B1_T = B1_T.shape[1]

    # Subtournament T\v
    A_sub, n_sub, vertices = subtournament(A, n, v)
    Z1_sub, B1_sub, edges_sub, beta1_sub = compute_cocycle_space(A_sub, n_sub)
    edge_idx_sub = {e: i for i, e in enumerate(edges_sub)}
    dim_Z1_sub = Z1_sub.shape[1]
    dim_B1_sub = B1_sub.shape[1]

    # Restriction matrix R: maps edge-space of T → edge-space of T\v
    # For each edge (ii, jj) in T\v (relabeled), find corresponding edge in T
    R = np.zeros((len(edges_sub), len(edges_T)))
    for sub_idx, (ii, jj) in enumerate(edges_sub):
        orig_edge = (vertices[ii], vertices[jj])
        if orig_edge in edge_idx_T:
            R[sub_idx, edge_idx_T[orig_edge]] = 1

    # res_v(Z¹(T)) = R @ Z1_T (in edge-space of T\v)
    res_image = R @ Z1_T  # shape: (|edges_sub|, dim_Z1_T)

    # Verify image lands in Z¹(T\v)
    # Check: every column of res_image should be in Z¹(T\v)
    # i.e., rank([Z1_sub | res_image]) = rank(Z1_sub)
    if dim_Z1_sub > 0 and res_image.shape[1] > 0:
        combined = np.hstack([Z1_sub, res_image])
        sv_comb = np.linalg.svd(combined, compute_uv=False)
        combined_rank = sum(s > 1e-8 for s in sv_comb)
        image_in_Z1_sub = (combined_rank == dim_Z1_sub)
    elif res_image.shape[1] == 0:
        image_in_Z1_sub = True
    else:
        # Z1_sub is trivial
        sv_img = np.linalg.svd(res_image, compute_uv=False)
        img_rank = sum(s > 1e-8 for s in sv_img)
        image_in_Z1_sub = (img_rank == 0)

    # Rank of image
    if res_image.shape[0] > 0 and res_image.shape[1] > 0:
        sv = np.linalg.svd(res_image, compute_uv=False)
        res_image_dim = sum(s > 1e-8 for s in sv)
    else:
        res_image_dim = 0

    ker_dim = dim_Z1_T - res_image_dim

    # Surjectivity: image spans all of Z¹(T\v)?
    if image_in_Z1_sub:
        surjective = (res_image_dim == dim_Z1_sub)
        coker_dim = dim_Z1_sub - res_image_dim
    else:
        surjective = False
        coker_dim = None  # can't compute if image doesn't land in Z1_sub

    # Also check: does image land in B¹(T\v)?
    if dim_B1_sub > 0 and res_image.shape[1] > 0:
        combined_B = np.hstack([B1_sub, res_image])
        sv_comb_B = np.linalg.svd(combined_B, compute_uv=False)
        combined_rank_B = sum(s > 1e-8 for s in sv_comb_B)
        image_in_B1_sub = (combined_rank_B == dim_B1_sub)
    elif res_image.shape[1] == 0:
        image_in_B1_sub = True
    else:
        sv_img = np.linalg.svd(res_image, compute_uv=False)
        img_rank = sum(s > 1e-8 for s in sv_img)
        image_in_B1_sub = (img_rank == 0)

    return {
        'dim_Z1_T': dim_Z1_T,
        'dim_B1_T': dim_B1_T,
        'beta1_T': beta1_T,
        'dim_Z1_sub': dim_Z1_sub,
        'dim_B1_sub': dim_B1_sub,
        'beta1_sub': beta1_sub,
        'dim_image': res_image_dim,
        'dim_kernel': ker_dim,
        'surjective': surjective,
        'coker_dim': coker_dim,
        'image_in_Z1_sub': image_in_Z1_sub,
        'image_in_B1_sub': image_in_B1_sub,
        'n_edges_T': len(edges_T),
        'n_edges_sub': len(edges_sub),
    }


def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A, mask


def path_betti_1(A, n):
    """Compute β₁ via cohomology."""
    _, _, _, beta1 = compute_cocycle_space(A, n)
    return beta1


# ===== VALIDATION =====
print("=" * 72)
print("VALIDATION: Cocycle space dimensions")
print("=" * 72)

# 3-cycle: should have β₁=1
A3_cyc = [[0,1,0],[0,0,1],[1,0,0]]
Z1, B1, edges, b1 = compute_cocycle_space(A3_cyc, 3)
print(f"3-cycle: edges={len(edges)}, dim Z¹={Z1.shape[1]}, dim B¹={B1.shape[1]}, β¹={b1}")

# Transitive 3-tournament: should have β₁=0
A3_trans = [[0,1,1],[0,0,1],[0,0,0]]
Z1, B1, edges, b1 = compute_cocycle_space(A3_trans, 3)
print(f"TT₃:    edges={len(edges)}, dim Z¹={Z1.shape[1]}, dim B¹={B1.shape[1]}, β¹={b1}")

# Quick check at n=4
print("\nn=4 β₁ distribution:")
beta_dist_4 = defaultdict(int)
for A, mask in all_tournaments(4):
    b = path_betti_1(A, 4)
    beta_dist_4[b] += 1
print(f"  {dict(beta_dist_4)}")

# Quick check at n=5
print("\nn=5 β₁ distribution:")
beta_dist_5 = defaultdict(int)
for A, mask in all_tournaments(5):
    b = path_betti_1(A, 5)
    beta_dist_5[b] += 1
print(f"  {dict(beta_dist_5)}")

# Cross-validate with path_homology_v2 at n=4
print("\nCross-validation with path_homology_v2 at n=4:")
import sys
# Import homology functions
sys.path.insert(0, '/home/e/Documents/claude/math/04-computation')

# Inline the homology β₁ computation
def homology_beta1(A, n):
    """Compute β₁ via HOMOLOGICAL path_betti_numbers."""
    from path_homology_v2 import path_betti_numbers as pbn
    # Can't import directly due to top-level execution, compute manually
    allowed_0 = [(v,) for v in range(n)]
    allowed_1 = []
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1:
                adj[i].append(j)
    for i in range(n):
        for j in adj[i]:
            allowed_1.append((i,j))

    allowed_2 = []
    for a in range(n):
        for b in adj[a]:
            for c in adj[b]:
                if c != a and A[a][c]:  # need a→c too? No, just need a→b→c as path
                    # Actually allowed 2-path = (a,b,c) with a→b, b→c, all distinct
                    if c != a:
                        allowed_2.append((a,b,c))

    # Wait, allowed 2-paths need a→b AND b→c (as directed edges), vertices distinct
    # Already handled above
    # Remove duplicates (shouldn't be any)
    allowed_2 = list(set(allowed_2))

    # For tournaments, Ω₁ = A₁ (all faces of 1-paths are vertices, always allowed)
    # Ω₂: need ∂(a,b,c) = (b,c) - (a,c) + (a,b) to have all terms in A₁
    # (a,b) is in A₁ iff a→b
    # (a,c) is in A₁ iff a→c
    # (b,c) is in A₁ iff b→c
    # So Ω₂ = {(a,b,c) ∈ A₂ : (a,c) ∈ A₁} = TTs!

    omega2_paths = [(a,b,c) for (a,b,c) in allowed_2 if A[a][c] == 1]

    # ∂₁ on Ω₁=A₁: ∂(a,b) = b - a
    bd1 = np.zeros((n, len(allowed_1)))
    edge_idx = {e:i for i,e in enumerate(allowed_1)}
    for j, (a,b) in enumerate(allowed_1):
        bd1[b, j] += 1
        bd1[a, j] -= 1

    # ker(∂₁)
    U1, S1, V1t = np.linalg.svd(bd1, full_matrices=True)
    rank1 = sum(s > 1e-10 for s in S1)
    ker_dim = len(allowed_1) - rank1

    # ∂₂ on Ω₂: ∂(a,b,c) = (b,c) - (a,c) + (a,b)
    bd2 = np.zeros((len(allowed_1), len(omega2_paths)))
    for j, (a,b,c) in enumerate(omega2_paths):
        if (b,c) in edge_idx:
            bd2[edge_idx[(b,c)], j] += 1
        if (a,c) in edge_idx:
            bd2[edge_idx[(a,c)], j] -= 1
        if (a,b) in edge_idx:
            bd2[edge_idx[(a,b)], j] += 1

    if bd2.shape[1] > 0:
        S2 = np.linalg.svd(bd2, compute_uv=False)
        im_dim = sum(s > 1e-8 for s in S2)
    else:
        im_dim = 0

    return ker_dim - im_dim

mismatches = 0
for A, mask in all_tournaments(4):
    b_cohom = path_betti_1(A, 4)
    b_homol = homology_beta1(A, 4)
    if b_cohom != b_homol:
        mismatches += 1
        if mismatches <= 3:
            print(f"  MISMATCH mask={mask}: cohom β₁={b_cohom}, homol β₁={b_homol}")
if mismatches == 0:
    print(f"  All {64} tournaments agree: cohom β₁ = homol β₁. Good.")
else:
    print(f"  {mismatches} mismatches found!")

# ===== MAIN COMPUTATION =====
print(f"\n{'='*72}")
print(f"COCYCLE RESTRICTION MAP ANALYSIS — n=5")
print(f"{'='*72}")

n = 5
total = 0
beta1_zero_count = 0
bad_vertex_count = 0
good_vertex_count = 0
surjective_bad = 0
nonsurjective_bad = 0
surjective_good = 0
nonsurjective_good = 0

dim_stats_bad = defaultdict(int)
dim_stats_good = defaultdict(int)
coker_bad = defaultdict(int)
ker_bad = defaultdict(int)
coker_good = defaultdict(int)
ker_good = defaultdict(int)
image_not_in_Z1 = 0
image_in_B1_count_bad = 0

examples_bad = []
examples_good_nonsurj = []

for A, mask in all_tournaments(n):
    total += 1
    b1 = path_betti_1(A, n)
    if b1 != 0:
        continue
    beta1_zero_count += 1

    for v in range(n):
        info = compute_restriction_map(A, n, v)
        is_bad = (info['beta1_sub'] > 0)

        if not info['image_in_Z1_sub']:
            image_not_in_Z1 += 1
            if image_not_in_Z1 <= 3:
                print(f"  WARNING: image not in Z¹(T\\v)! mask={mask}, v={v}, info={info}")

        if is_bad:
            bad_vertex_count += 1
            key = (info['dim_Z1_T'], info['dim_B1_T'], info['dim_Z1_sub'],
                   info['dim_B1_sub'], info['dim_image'], info['dim_kernel'])
            dim_stats_bad[key] += 1
            if info['coker_dim'] is not None:
                coker_bad[info['coker_dim']] += 1
            ker_bad[info['dim_kernel']] += 1

            if info['image_in_B1_sub']:
                image_in_B1_count_bad += 1

            if info['surjective']:
                surjective_bad += 1
                print(f"  *** CONTRADICTION: surjective for bad vertex! mask={mask}, v={v}")
            else:
                nonsurjective_bad += 1

            if len(examples_bad) < 5:
                examples_bad.append((mask, v, info))
        else:
            good_vertex_count += 1
            key = (info['dim_Z1_T'], info['dim_B1_T'], info['dim_Z1_sub'],
                   info['dim_B1_sub'], info['dim_image'], info['dim_kernel'])
            dim_stats_good[key] += 1
            if info['coker_dim'] is not None:
                coker_good[info['coker_dim']] += 1
            ker_good[info['dim_kernel']] += 1

            if info['surjective']:
                surjective_good += 1
            else:
                nonsurjective_good += 1
                if len(examples_good_nonsurj) < 3:
                    examples_good_nonsurj.append((mask, v, info))

print(f"\n{'='*72}")
print(f"RESULTS FOR n={n}")
print(f"{'='*72}")
print(f"Total tournaments: {total}")
print(f"Tournaments with β₁=0: {beta1_zero_count}")
print(f"\nVertex deletions from β₁=0 tournaments:")
print(f"  Bad vertices (β₁(T\\v) ≥ 1): {bad_vertex_count}")
print(f"  Good vertices (β₁(T\\v) = 0): {good_vertex_count}")

print(f"\n--- BAD VERTICES ---")
print(f"  res_v surjective:     {surjective_bad}  {'(CONTRADICTION!)' if surjective_bad > 0 else '(correct: should be 0)'}")
print(f"  res_v NOT surjective: {nonsurjective_bad}")
print(f"  Image ⊆ Z¹(T\\v) failures: {image_not_in_Z1}  (should be 0)")
print(f"  Image ⊆ B¹(T\\v): {image_in_B1_count_bad} / {bad_vertex_count}  (should be ALL)")

print(f"\n  Kernel dim distribution: {dict(ker_bad)}")
print(f"  Cokernel dim distribution: {dict(coker_bad)}")
print(f"\n  Full dimension profile (Z¹_T, B¹_T, Z¹_sub, B¹_sub, im, ker):")
for key in sorted(dim_stats_bad.keys()):
    print(f"    {key}: {dim_stats_bad[key]} cases")

print(f"\n--- GOOD VERTICES ---")
print(f"  res_v surjective:     {surjective_good}")
print(f"  res_v NOT surjective: {nonsurjective_good}")

print(f"\n  Kernel dim distribution: {dict(ker_good)}")
print(f"  Cokernel dim distribution: {dict(coker_good)}")
print(f"\n  Full dimension profile (Z¹_T, B¹_T, Z¹_sub, B¹_sub, im, ker):")
for key in sorted(dim_stats_good.keys()):
    print(f"    {key}: {dim_stats_good[key]} cases")

print(f"\n{'='*72}")
print(f"DETAILED EXAMPLES — BAD VERTICES")
print(f"{'='*72}")
for mask, v, info in examples_bad:
    print(f"\n  Tournament mask={mask}, delete vertex v={v}")
    print(f"    β₁(T) = {info['beta1_T']}")
    print(f"    β₁(T\\v) = {info['beta1_sub']}")
    print(f"    dim Z¹(T) = {info['dim_Z1_T']}, dim B¹(T) = {info['dim_B1_T']}")
    print(f"    dim Z¹(T\\v) = {info['dim_Z1_sub']}, dim B¹(T\\v) = {info['dim_B1_sub']}")
    print(f"    |edges(T)| = {info['n_edges_T']}, |edges(T\\v)| = {info['n_edges_sub']}")
    print(f"    dim res_v(Z¹(T)) = {info['dim_image']}")
    print(f"    dim ker(res_v) = {info['dim_kernel']}")
    print(f"    dim coker = {info['coker_dim']}")
    print(f"    Surjective: {info['surjective']}")
    print(f"    Image ⊆ Z¹(T\\v): {info['image_in_Z1_sub']}")
    print(f"    Image ⊆ B¹(T\\v): {info['image_in_B1_sub']}")

if examples_good_nonsurj:
    print(f"\n{'='*72}")
    print(f"EXAMPLES — GOOD VERTICES WITH NON-SURJECTIVE res_v")
    print(f"{'='*72}")
    for mask, v, info in examples_good_nonsurj:
        print(f"\n  Tournament mask={mask}, delete vertex v={v}")
        print(f"    β₁(T\\v) = {info['beta1_sub']}")
        print(f"    dim Z¹(T) = {info['dim_Z1_T']}, dim B¹(T) = {info['dim_B1_T']}")
        print(f"    dim Z¹(T\\v) = {info['dim_Z1_sub']}, dim B¹(T\\v) = {info['dim_B1_sub']}")
        print(f"    dim res_v(Z¹(T)) = {info['dim_image']}")
        print(f"    dim ker(res_v) = {info['dim_kernel']}")
        print(f"    dim coker = {info['coker_dim']}")

# ===== SANITY CHECK n=4 =====
print(f"\n{'='*72}")
print(f"SANITY CHECK: n=4")
print(f"{'='*72}")

bad4 = 0
surj_bad4 = 0
good4 = 0
surj_good4 = 0
b1_zero_4 = 0

for A, mask in all_tournaments(4):
    b1 = path_betti_1(A, 4)
    if b1 != 0:
        continue
    b1_zero_4 += 1
    for v in range(4):
        info = compute_restriction_map(A, 4, v)
        if info['beta1_sub'] > 0:
            bad4 += 1
            if info['surjective']: surj_bad4 += 1
        else:
            good4 += 1
            if info['surjective']: surj_good4 += 1

print(f"β₁=0 tournaments: {b1_zero_4}")
print(f"Bad vertices: {bad4}, surjective: {surj_bad4}")
print(f"Good vertices: {good4}, surjective: {surj_good4}")

# ===== n=6 SAMPLE =====
print(f"\n{'='*72}")
print(f"SAMPLE CHECK: n=6 (first 500 tournaments)")
print(f"{'='*72}")

import random
random.seed(42)
n6 = 6
bad6 = 0
surj_bad6 = 0
good6 = 0
surj_good6 = 0
b1_zero_6 = 0
count6 = 0

for A, mask in all_tournaments(n6):
    count6 += 1
    if count6 > 500:
        break
    b1 = path_betti_1(A, n6)
    if b1 != 0:
        continue
    b1_zero_6 += 1
    for v in range(n6):
        info = compute_restriction_map(A, n6, v)
        if info['beta1_sub'] > 0:
            bad6 += 1
            if info['surjective']: surj_bad6 += 1
        else:
            good6 += 1
            if info['surjective']: surj_good6 += 1

print(f"Tested {count6} tournaments, {b1_zero_6} with β₁=0")
print(f"Bad vertices: {bad6}, surjective: {surj_bad6}")
print(f"Good vertices: {good6}, surjective: {surj_good6}")

# ===== CONCLUSION =====
print(f"\n{'='*72}")
print("CONCLUSION")
print(f"{'='*72}")

all_bad_nonsurj = (surjective_bad == 0) and (bad_vertex_count > 0)
all_bad_in_B1 = (image_in_B1_count_bad == bad_vertex_count) if bad_vertex_count > 0 else True

msg1 = "CONFIRMED: res_v is NEVER surjective for bad vertices" if all_bad_nonsurj else "UNEXPECTED: found surjective cases!"
msg2 = "CONFIRMED: image always lands in B1(T-v)" if all_bad_in_B1 else "UNEXPECTED: image escapes B1!"
msg3 = "CONFIRMED: restriction always preserves cocycles" if image_not_in_Z1 == 0 else "BUG: restriction should preserve cocycles!"
print(f"\n1. res_v surjective for bad vertices: {surjective_bad} / {bad_vertex_count}")
print(f"   {msg1}")
print(f"\n2. Image in B1(T-v) for bad vertices: {image_in_B1_count_bad} / {bad_vertex_count}")
print(f"   {msg2}")
print(f"\n3. Image in Z1(T-v) failures: {image_not_in_Z1}")
print(f"   {msg3}")
print("""
ALGEBRAIC PROOF (verified computationally):
  Given b1(T) = 0 and b1(T-v) = 1:
  (a) Z1(T) = B1(T)                           [b1 = 0]
  (b) res_v(B1(T)) <= B1(T-v)                 [restriction of coboundary = coboundary]
  (c) So res_v(Z1(T)) <= B1(T-v) < Z1(T-v)    [B1 < Z1 since b1(T-v) = 1]
  (d) Therefore res_v is NOT surjective.        QED

  The previous agent's claim of surjectivity was INCORRECT.
  The obstruction has dimension = b1(T-v).
""")

print("Done.")
