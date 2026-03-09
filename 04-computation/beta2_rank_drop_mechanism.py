#!/usr/bin/env python3
"""beta2_rank_drop_mechanism.py - Understand WHY rank(d_2) drops by n-2 or n-1

KEY EQUATION:
  b1(T\v) <= b1(T) iff rk(d_2(T)) - rk(d_2(T\v)) <= n-2

dim(Z_1(T)) = C(n-1,2), dim(Z_1(T\v)) = C(n-2,2)
C(n-1,2) - C(n-2,2) = n-2.

So: b1(T\v) - b1(T) = [C(n-2,2) - rk(d_2(T\v))] - [C(n-1,2) - rk(d_2(T))]
                     = rk(d_2(T)) - rk(d_2(T\v)) - (n-2)

Bad vertex: rk(d_2(T)) - rk(d_2(T\v)) = n-1 (one more than threshold)

QUESTION: What is rk(d_2(T)) - rk(d_2(T\v)) in terms of graph structure?

Omega_2(T) contains:
1. Elements of Omega_2(T\v) (not involving v) -- contribute rk(d_2(T\v))
2. Elements involving v -- contribute additional rank

The "rank increment" from v-involving elements is:
  delta_rk = rk(d_2(T)) - rk(d_2(T\v))

If Omega_2(T\v) is a subcomplex (which it is), then:
  delta_rk <= dim(Omega_2(T)) - dim(Omega_2(T\v)) = dim(R_2) (relative complex dim)

And delta_rk >= 0 (adding generators can't decrease rank).

For b1 monotonicity: need delta_rk <= n-2.
Bad case: delta_rk = n-1.

INSIGHT: d_2 maps Omega_2 into Z_1. The additional rank from v-involving
elements is the dimension of new boundary directions they create in Z_1.
Z_1(T) has dimension C(n-1,2) and Z_1(T\v) has dimension C(n-2,2).
The "v-involving" part of Z_1 has dimension C(n-1,2) - C(n-2,2) = n-2
(since Z_1 is constant-dimensional).

Wait -- is Z_1(T\v) embedded in Z_1(T)? Not directly, since they live in
different A_1 spaces. But via the inclusion T\v -> T, Z_1(T\v) maps into Z_1(T).

The image of Z_1(T\v) in Z_1(T) has dimension at most C(n-2,2).
And Z_1(T) has dimension C(n-1,2).
The "new" directions in Z_1(T) not covered by Z_1(T\v) have dimension n-2.

If im(d_2(T\v)) maps into im(d_2(T)) via inclusion, then:
delta_rk = rk(d_2(T)) - rk(d_2(T\v)) = dim(new boundary directions from v-elements)

These new directions must lie in the "n-2 dimensional complement" of im(d_2(T\v)).
So delta_rk <= n-2 unless some d_2(T\v) boundaries become "redundant" in T.

Hmm, but we observed delta_rk = n-1 > n-2 in bad cases. How?

The issue: im(d_2(T\v)) might NOT inject into im(d_2(T)). Some T\v boundaries
might become boundaries of different 2-chains in T, freeing up rank for v-involving
elements.

Let me investigate this computationally.

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved

random.seed(42)


def all_tournaments_gen(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0] * n for _ in range(n)]
        for idx, (i, j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


def analyze_rank_drop(A, n, v):
    """Detailed analysis of how rk(d_2) changes under deletion of v."""
    # T data
    paths2_T = enumerate_allowed_paths(A, n, 2)
    paths1_T = enumerate_allowed_paths(A, n, 1)
    paths0_T = [(i,) for i in range(n)]

    omega2_T = compute_omega_basis(A, n, 2, paths2_T, paths1_T)
    dim_O2_T = omega2_T.shape[1] if omega2_T.ndim == 2 else 0

    omega1_T = compute_omega_basis(A, n, 1, paths1_T, paths0_T)
    dim_O1_T = omega1_T.shape[1] if omega1_T.ndim == 2 else 0

    if dim_O2_T == 0:
        return None

    D2_T = build_full_boundary_matrix([tuple(p) for p in paths2_T],
                                      [tuple(p) for p in paths1_T])
    D2_om_T = D2_T @ omega2_T  # in A_1(T) coords
    sv = np.linalg.svd(D2_om_T, compute_uv=False)
    rk_d2_T = int(sum(s > 1e-8 for s in sv))

    # T\v data
    others = [x for x in range(n) if x != v]
    B_sub, vlist = get_induced(A, n, others)
    np2 = len(others)

    paths2_Tv = enumerate_allowed_paths(B_sub, np2, 2)
    paths1_Tv = enumerate_allowed_paths(B_sub, np2, 1)
    paths0_Tv = [(i,) for i in range(np2)]

    omega2_Tv = compute_omega_basis(B_sub, np2, 2, paths2_Tv, paths1_Tv) if paths2_Tv else np.zeros((0,0))
    dim_O2_Tv = omega2_Tv.shape[1] if omega2_Tv.ndim == 2 and omega2_Tv.shape[1] > 0 else 0

    if dim_O2_Tv > 0:
        D2_Tv = build_full_boundary_matrix([tuple(p) for p in paths2_Tv],
                                           [tuple(p) for p in paths1_Tv])
        D2_om_Tv = D2_Tv @ omega2_Tv
        sv2 = np.linalg.svd(D2_om_Tv, compute_uv=False)
        rk_d2_Tv = int(sum(s > 1e-8 for s in sv2))
    else:
        rk_d2_Tv = 0

    drop = rk_d2_T - rk_d2_Tv
    d_out = sum(A[v])

    # Classify Omega_2(T) elements by whether they involve v
    path2_T_idx = {tuple(p): i for i, p in enumerate(paths2_T)}

    # Which 2-paths involve v?
    v_involving = set()
    for i, p in enumerate(paths2_T):
        if v in p:
            v_involving.add(i)

    # Number of v-involving and non-v Omega_2 basis vectors
    # (This is approximate -- Omega_2 basis vectors might mix both types)

    # Instead: what is dim(R_2) = dim(Omega_2(T)) - dim(Omega_2(T\v))?
    dim_R2 = dim_O2_T - dim_O2_Tv

    # Embed im(d_2(T\v)) into A_1(T) and compute its dimension
    path1_T_idx = {tuple(p): i for i, p in enumerate(paths1_T)}

    if dim_O2_Tv > 0:
        embed1 = np.zeros((len(paths1_T), len(paths1_Tv)))
        for j, p_Tv in enumerate(paths1_Tv):
            p_T = tuple(vlist[x] for x in p_Tv)
            if p_T in path1_T_idx:
                embed1[path1_T_idx[p_T], j] = 1

        # im(d_2(T\v)) embedded in A_1(T)
        im_d2_Tv_emb = embed1 @ D2_om_Tv
        rk_emb = np.linalg.matrix_rank(im_d2_Tv_emb, tol=1e-8)
    else:
        im_d2_Tv_emb = np.zeros((len(paths1_T), 0))
        rk_emb = 0

    # Does im(d_2(T\v)) inject into im(d_2(T))?
    # i.e., is rk_emb = rk_d2_Tv?
    injects = (rk_emb == rk_d2_Tv)

    # Combine im(d_2(T\v)) with im(d_2(T))
    if rk_emb > 0 and rk_d2_T > 0:
        combined = np.hstack([im_d2_Tv_emb, D2_om_T])
        rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)
        # im(d_2(T\v)) subset im(d_2(T))?
        im_Tv_in_T = (rk_combined == rk_d2_T)
    else:
        im_Tv_in_T = True

    return {
        'd_out': d_out,
        'dim_O2_T': dim_O2_T, 'dim_O2_Tv': dim_O2_Tv,
        'dim_R2': dim_R2,
        'rk_d2_T': rk_d2_T, 'rk_d2_Tv': rk_d2_Tv,
        'drop': drop,
        'good': drop <= n - 2,
        'rk_emb': rk_emb,
        'emb_injects': injects,
        'im_Tv_in_T': im_Tv_in_T,
    }


# ============================================================
# Part 1: Analyze good vs bad at n=5
# ============================================================
print("=" * 70)
print("RANK DROP MECHANISM at n=5")
print("=" * 70)

n = 5
good_data = []
bad_data = []

for A in all_tournaments_gen(n):
    for v in range(n):
        result = analyze_rank_drop(A, n, v)
        if result is None:
            continue
        if result['good']:
            good_data.append(result)
        else:
            bad_data.append(result)

print(f"n=5: good={len(good_data)}, bad={len(bad_data)}")

print(f"\nGOOD vertices:")
print(f"  drop values: {Counter(r['drop'] for r in good_data)}")
print(f"  dim_R2: {Counter(r['dim_R2'] for r in good_data)}")
print(f"  emb_injects: {Counter(r['emb_injects'] for r in good_data)}")
print(f"  im(d2_Tv) subset im(d2_T): {Counter(r['im_Tv_in_T'] for r in good_data)}")

print(f"\nBAD vertices:")
print(f"  drop values: {Counter(r['drop'] for r in bad_data)}")
print(f"  dim_R2: {Counter(r['dim_R2'] for r in bad_data)}")
print(f"  emb_injects: {Counter(r['emb_injects'] for r in bad_data)}")
print(f"  im(d2_Tv) subset im(d2_T): {Counter(r['im_Tv_in_T'] for r in bad_data)}")

# KEY: Do bad vertices have im(d2_Tv) NOT subset im(d2_T)?
# This would mean some T\v boundaries are NOT boundaries in T,
# which seems paradoxical. Let me check.

print(f"\n{'=' * 70}")
print("DETAILED BAD VERTEX EXAMPLE")
print("=" * 70)

for A in all_tournaments_gen(n):
    for v in range(n):
        result = analyze_rank_drop(A, n, v)
        if result and not result['good']:
            d_out = result['d_out']
            scores = sorted([sum(row) for row in A])
            print(f"T scores={scores}, v={v}, d_out={d_out}")
            print(f"  dim_O2(T)={result['dim_O2_T']}, dim_O2(T\\v)={result['dim_O2_Tv']}")
            print(f"  dim_R2={result['dim_R2']}")
            print(f"  rk_d2(T)={result['rk_d2_T']}, rk_d2(T\\v)={result['rk_d2_Tv']}")
            print(f"  drop={result['drop']}, threshold={n-2}")
            print(f"  emb_injects={result['emb_injects']}")
            print(f"  im(d2_Tv) subset im(d2_T): {result['im_Tv_in_T']}")

            # What are the v-involving 2-paths?
            P = [a for a in range(n) if a != v and A[v][a] == 1]
            Q = [b for b in range(n) if b != v and A[b][v] == 1]
            print(f"  P(v)={P}, Q(v)={Q}")

            # Count allowed 2-paths involving v
            paths2 = enumerate_allowed_paths(A, n, 2)
            v_paths = [p for p in paths2 if v in p]
            print(f"  #2-paths involving v: {len(v_paths)}/{len(paths2)}")

            break
    else:
        continue
    break


# ============================================================
# Part 2: What EXACTLY makes a vertex bad?
# ============================================================
print(f"\n{'=' * 70}")
print("STRUCTURAL CAUSE OF BAD VERTICES")
print("=" * 70)

# The rank drop is rk(d_2(T)) - rk(d_2(T\v)).
# This equals: dim(new boundaries from v-involving Omega_2 elements)
#              that are linearly independent of im(d_2(T\v)).

# Since im(d_2(T\v)) has dimension rk(d_2(T\v)) and lives in Z_1(T\v)
# (embedded in Z_1(T)), and Z_1(T) has n-2 more dimensions than Z_1(T\v),
# the new boundaries can span at most n-2 of these new dimensions PLUS
# some redundant directions in the old im(d_2(T\v)) space.

# Wait, that's wrong. Let me reconsider.
# im(d_2(T)) = im(d_2 from Omega_2(T)) in Z_1(T) (subspace of Omega_1(T))
# im(d_2(T\v)) is in Z_1(T\v) (subspace of Omega_1(T\v))
# Embedded: im(d_2(T\v)) maps to a subspace of Z_1(T)

# Z_1(T) has dimension C(n-1,2) and Z_1(T\v) embedded has dimension C(n-2,2).
# The "new" part of Z_1(T) (not from Z_1(T\v)) has dimension n-2.

# The v-involving Omega_2 elements map to Z_1(T).
# Their images can have components in BOTH the old Z_1(T\v) part AND
# the new n-2 dimensional part.

# If a v-involving element maps entirely to the new part: contributes to delta_rk
# If it maps to the old part: only contributes if it's outside im(d_2(T\v))

# For BAD case (drop = n-1):
# The v-involving elements contribute n-1 independent directions.
# But the new part only has n-2 dimensions.
# So at least ONE new direction must come from the OLD im(d_2(T\v)) complement.
# i.e., a v-involving element creates a boundary that "overlaps" with Z_1(T\v)
# but is NOT in im(d_2(T\v)).

# This means: there's a 1-cycle z in Z_1(T\v) such that:
# z is NOT a boundary in T\v, but z IS a boundary in T (from v-involving 2-chain)
# This is exactly a non-trivial element of ker(i_*: H_1(T\v) -> H_1(T))!

print("""
MECHANISM:
  Bad vertex v (drop = n-1): There exist 1-cycles z in Z_1(T\\v) that:
  - Are NOT boundaries in T\\v (z not in im(d_2(T\\v)))
  - ARE boundaries in T (z in im(d_2(T))) via v-involving 2-chains

  These are exactly the elements of ker(i_*: H_1(T\\v) -> H_1(T)).
  By LES: dim(ker(i_*)) = dim(H_2(T,T\\v)) - beta_2(T).
  If beta_2(T) = 0: dim(ker(i_*)) = dim(H_2(T,T\\v)).

  So: drop = n-2 + dim(ker(i_*)) = n-2 + dim(H_2(T,T\\v))
  Bad means: H_2(T,T\\v) != 0.

  BUT: we can't ALWAYS have H_2(T,T\\v) = 0 for all v.
  The key is: there EXISTS v with H_2(T,T\\v) = 0.

  This is HYP-278: exists v with b1(T\\v) <= b1(T).
""")


# ============================================================
# Part 3: Score distribution of bad vertices
# ============================================================
print("=" * 70)
print("SCORE ANALYSIS OF BAD VERTICES")
print("=" * 70)

n = 5
bad_vertex_positions = Counter()

for A in all_tournaments_gen(n):
    out_degs = [sum(A[v]) for v in range(n)]
    sorted_degs = sorted(range(n), key=lambda x: out_degs[x])

    for v in range(n):
        result = analyze_rank_drop(A, n, v)
        if result and not result['good']:
            # Position of v in sorted order (0 = min degree, n-1 = max)
            rank = sorted_degs.index(v)
            bad_vertex_positions[rank] += 1

total_bad = sum(bad_vertex_positions.values())
print(f"n=5: {total_bad} bad vertices")
print(f"Position in sorted degree order (0=min, {n-1}=max):")
for pos in range(n):
    cnt = bad_vertex_positions.get(pos, 0)
    print(f"  position {pos}: {cnt} ({100*cnt/total_bad:.1f}%)")


# ============================================================
# Part 4: Verify: is every tournament with b1=0 "tight"?
# i.e., MOST vertices are good when b1=0
# ============================================================
print(f"\n{'=' * 70}")
print("TIGHTNESS: FRACTION OF GOOD VERTICES")
print("=" * 70)

for n in [5, 6]:
    gen = list(all_tournaments_gen(n))
    good_fracs_b1_0 = []
    good_fracs_b1_1 = []

    for A in gen:
        # Compute b1(T) quickly
        paths2 = enumerate_allowed_paths(A, n, 2)
        paths1 = enumerate_allowed_paths(A, n, 1)
        paths0 = [(i,) for i in range(n)]

        omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
        dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

        omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
        dim_O1 = omega1.shape[1] if omega1.ndim == 2 else 0

        D1 = build_full_boundary_matrix([tuple(p) for p in paths1], paths0)
        D1_om = D1 @ omega1
        sv1 = np.linalg.svd(D1_om, compute_uv=False)
        rk_d1 = int(sum(s > 1e-8 for s in sv1))

        if dim_O2 > 0:
            D2 = build_full_boundary_matrix([tuple(p) for p in paths2],
                                            [tuple(p) for p in paths1])
            D2_om = D2 @ omega2
            sv2 = np.linalg.svd(D2_om, compute_uv=False)
            rk_d2 = int(sum(s > 1e-8 for s in sv2))
        else:
            rk_d2 = 0

        b1 = dim_O1 - rk_d1 - rk_d2

        good_v = 0
        for v in range(n):
            result = analyze_rank_drop(A, n, v)
            if result and result['good']:
                good_v += 1

        if b1 == 0:
            good_fracs_b1_0.append(good_v / n)
        elif b1 == 1:
            good_fracs_b1_1.append(good_v / n)

    print(f"\nn={n}:")
    if good_fracs_b1_0:
        print(f"  b1=0: avg good frac = {np.mean(good_fracs_b1_0):.4f}, "
              f"min = {min(good_fracs_b1_0):.4f}")
        frac_dist = Counter(round(f, 2) for f in good_fracs_b1_0)
        print(f"    distribution: {dict(sorted(frac_dist.items()))}")
    if good_fracs_b1_1:
        print(f"  b1=1: avg good frac = {np.mean(good_fracs_b1_1):.4f}, "
              f"min = {min(good_fracs_b1_1):.4f}")


print("\n\nDone.")
