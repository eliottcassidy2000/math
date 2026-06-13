#!/usr/bin/env python3
"""
beta2_algebraic_mechanism.py - Identify the algebraic mechanism for beta_2 = 0

Key insight: at p=2, Omega constraints come from non-allowed 1-paths.
For a tournament T on [n]:
- Non-allowed (a,c) means c->a (the "reverse" edge)
- Constraint: sum_{b: a->b, b->c} f(a,b,c) = 0

Combined with d_2(f) = 0:
- For each allowed (a,c) (a->c): sum_b f(a,b,c) - f(a,b,c') ... (complex)

Actually, let's be precise about what d_2 looks like.
For f = sum alpha_{a,b,c} * e_{a,b,c}, d_2(f) = sum alpha_{a,b,c} * (e_{b,c} - e_{a,c} + e_{a,b})

d_2(f) = 0 means: for each allowed 1-path (x,y):
    sum_{a,b,c: (b,c)=(x,y)} alpha_{a,b,c} - sum_{a,b,c: (a,c)=(x,y)} alpha_{a,b,c} + sum_{a,b,c: (a,b)=(x,y)} alpha_{a,b,c} = 0

For allowed (x,y) (i.e., x->y in T):
    sum_{a: a->x, x->y in T} alpha_{a,x,y} - sum_{b: x->b->y} alpha_{x,b,y} + sum_{c: x->y->c} alpha_{x,y,c} = 0

This is a "flow through (x,y)" equation.

Omega constraint for non-allowed (a,c) (i.e., c->a):
    sum_{b: a->b->c} alpha_{a,b,c} = 0

This says: for each NON-EDGE (a,c) with a-/->c, the "flow through the
non-allowed face (a,c)" vanishes.

KEY QUESTION: Do these two constraints (d_2=0 and Omega) together force
f to be in im(d_3|Omega_3)?

Let's analyze the constraint structure dimension-by-dimension.

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time
import numpy as np
from collections import defaultdict
from itertools import combinations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis, path_betti_numbers
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


# ============================================================
# PART 1: Explicit constraint analysis at n=5
# ============================================================
print("=" * 70)
print("CONSTRAINT STRUCTURE FOR beta_2 = 0")
print("=" * 70)

n = 5

# Pick a specific tournament with beta_1 = 1, ker(d_2) > 0
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=3)
    if len(betti) > 1 and betti[1] == 1:
        scores = tuple(sorted([sum(row) for row in A]))
        paths2 = enumerate_allowed_paths(A, n, 2)
        paths1 = enumerate_allowed_paths(A, n, 1)
        paths3 = enumerate_allowed_paths(A, n, 3)
        omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
        dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

        D2 = build_full_boundary_matrix(paths2, paths1)
        D2_omega = D2 @ omega2
        S2 = np.linalg.svd(D2_omega, compute_uv=False)
        rank_d2 = sum(s > 1e-8 for s in S2)
        ker_d2 = dim_O2 - rank_d2

        if ker_d2 > 0:
            print(f"\nExample: T#{bits} scores={scores} beta={list(betti)}")
            print(f"  dim(A_2)={len(paths2)}, dim(Omega_2)={dim_O2}, ker(d_2)={ker_d2}")
            print(f"  Allowed 2-paths ({len(paths2)}):")
            for p in paths2:
                a,b,c = p
                face_type = "TT" if A[a][c] else "3C"
                print(f"    {p} [{face_type}]")

            # Non-allowed 1-paths and their mediators
            print(f"\n  Non-allowed 1-paths and mediators:")
            for a in range(n):
                for c in range(n):
                    if a != c and A[a][c] == 0:  # c->a, so (a,c) is non-allowed
                        meds = [b for b in range(n) if b != a and b != c
                                and A[a][b] == 1 and A[b][c] == 1]
                        if meds:
                            print(f"    ({a},{c}): c->a edge, mediators={meds}")
                            for b in meds:
                                idx = paths2.index((a,b,c)) if (a,b,c) in paths2 else -1
                                print(f"      2-path ({a},{b},{c}) at index {idx}")

            # Get explicit ker(d_2) generator
            U2, S2f, V2t = np.linalg.svd(D2_omega, full_matrices=True)
            rank = sum(s > 1e-8 for s in S2f)
            ker_basis = V2t[rank:]  # rows = ker vectors in Omega_2 coords
            ker_A2 = omega2 @ ker_basis.T  # columns in A_2 coords

            print(f"\n  ker(d_2|Omega_2) generators ({ker_d2} total):")
            for g in range(ker_A2.shape[1]):
                v = ker_A2[:, g]
                print(f"    Generator #{g+1}:")
                for j, path in enumerate(paths2):
                    if abs(v[j]) > 1e-10:
                        a,b,c = path
                        face_type = "TT" if A[a][c] else "3C"
                        triple = tuple(sorted([a,b,c]))
                        print(f"      ({a},{b},{c}): {v[j]:+.6f} [{face_type}] triple={triple}")

            # Now check: does the generator satisfy the Omega constraint?
            # (It must, since it's in Omega_2, but let's verify explicitly)
            for g in range(ker_A2.shape[1]):
                v = ker_A2[:, g]
                print(f"\n  Verify Omega constraints for generator #{g+1}:")
                for a in range(n):
                    for c in range(n):
                        if a != c and A[a][c] == 0:
                            meds = [b for b in range(n) if b != a and b != c
                                    and A[a][b] == 1 and A[b][c] == 1]
                            if meds:
                                s = 0
                                for b in meds:
                                    if (a,b,c) in paths2:
                                        idx = paths2.index((a,b,c))
                                        s += v[idx]
                                print(f"    ({a},{c}): sum = {s:.10f} [{'OK' if abs(s) < 1e-8 else 'FAIL'}]")

            # And d_2 = 0:
            D2_full = build_full_boundary_matrix(paths2, paths1)
            for g in range(ker_A2.shape[1]):
                v = ker_A2[:, g]
                dv = D2_full @ v
                print(f"\n  d_2(generator #{g+1}) max component: {np.max(np.abs(dv)):.2e}")

            # Now: explicitly build the filling from Omega_3
            omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
            dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0
            print(f"\n  dim(Omega_3)={dim_O3}, dim(A_3)={len(paths3)}")

            if dim_O3 > 0:
                D3 = build_full_boundary_matrix(paths3, paths2)
                D3_omega = D3 @ omega3
                # Does the generator lie in im(D3_omega)?
                for g in range(ker_A2.shape[1]):
                    v = ker_A2[:, g]
                    # Solve D3_omega @ w = v
                    w, res, _, _ = np.linalg.lstsq(D3_omega, v, rcond=None)
                    err = np.max(np.abs(D3_omega @ w - v))
                    print(f"    Generator #{g+1} in im(d_3|Omega_3)? err={err:.2e}")

                    if err < 1e-8:
                        w_A3 = omega3 @ w
                        print(f"    Filling w in Omega_3 (A_3 coords):")
                        for j, path in enumerate(paths3):
                            if abs(w_A3[j]) > 1e-10:
                                print(f"      w[{path}] = {w_A3[j]:+.6f}")

            break  # Just analyze first example


# ============================================================
# PART 2: Dimension counting identity
# ============================================================
print(f"\n\n{'='*70}")
print("DIMENSION COUNTING: WHY ker(d_2) = im(d_3)")
print("=" * 70)

# For beta_2 = 0, we need: dim(ker d_2) = dim(im d_3)
# Let's track how these dimensions relate to tournament invariants.

n = 5
total = 1 << (n*(n-1)//2)

print(f"\nn=5 exhaustive: ker(d_2) vs im(d_3) + invariants")
data_5 = []
for bits in range(total):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    c3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
             if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j]))

    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths3 = enumerate_allowed_paths(A, n, 3)
    paths0 = enumerate_allowed_paths(A, n, 0)

    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

    # rank(d_2)
    D2 = build_full_boundary_matrix(paths2, paths1)
    D2_omega = D2 @ omega2
    S2 = np.linalg.svd(D2_omega, compute_uv=False)
    rank_d2 = sum(s > 1e-8 for s in S2) if len(S2) > 0 else 0
    ker_d2 = dim_O2 - rank_d2

    # rank(d_3|Omega_3) = im(d_3)
    if dim_O3 > 0 and len(paths2) > 0:
        D3 = build_full_boundary_matrix(paths3, paths2)
        D3_omega = D3 @ omega3
        S3 = np.linalg.svd(D3_omega, compute_uv=False)
        im_d3 = sum(s > 1e-8 for s in S3)
    else:
        im_d3 = 0

    # Non-allowed pairs with mediators (= Omega constraints)
    n_constraints = 0
    for a in range(n):
        for c in range(n):
            if a != c and A[a][c] == 0:
                meds = [b for b in range(n) if b != a and b != c
                        and A[a][b] == 1 and A[b][c] == 1]
                if len(meds) > 0:
                    n_constraints += 1

    data_5.append((scores, c3, len(paths2), dim_O2, rank_d2, ker_d2, dim_O3, im_d3, n_constraints))

# Group by score sequence
score_data = defaultdict(list)
for d in data_5:
    score_data[d[0]].append(d)

print(f"\nGrouped by score sequence:")
for scores in sorted(score_data.keys()):
    entries = score_data[scores]
    # Check if all same
    ker_vals = set(e[5] for e in entries)
    im_vals = set(e[7] for e in entries)
    dim_O2_vals = set(e[3] for e in entries)
    dim_O3_vals = set(e[6] for e in entries)
    c3_vals = set(e[1] for e in entries)

    e = entries[0]
    print(f"\n  scores={scores} ({len(entries)} tournaments)")
    print(f"    c3={c3_vals}")
    print(f"    dim(A_2)={set(e2[2] for e2 in entries)}")
    print(f"    dim(Omega_2)={dim_O2_vals}")
    print(f"    rank(d_2)={set(e2[4] for e2 in entries)}")
    print(f"    ker(d_2)={ker_vals}")
    print(f"    dim(Omega_3)={dim_O3_vals}")
    print(f"    im(d_3)={im_vals}")
    print(f"    #constraints={set(e2[8] for e2 in entries)}")
    print(f"    beta_2 = ker-im = {set(k-i for k,i in zip([e2[5] for e2 in entries], [e2[7] for e2 in entries]))}")


# ============================================================
# PART 3: What formula gives dim(Omega_2)?
# ============================================================
print(f"\n\n{'='*70}")
print("FORMULA FOR dim(Omega_2)")
print("=" * 70)

# We know dim(A_2) = C(n,3) + 2*c3
# What is dim(Omega_2) = dim(A_2) - #constraints?

print(f"\nn=5:")
for scores in sorted(score_data.keys()):
    entries = score_data[scores]
    for e in entries[:1]:
        c3 = e[1]
        dim_A2 = e[2]
        dim_O2 = e[3]
        n_constr = e[8]
        n_edges = n*(n-1)//2  # = 10 for n=5
        n_reverse = n_edges  # also 10 (each pair has exactly one direction)
        # Each pair (a,c) with c->a: mediators = #{b: a->b->c}
        print(f"  scores={e[0]}: dim(A_2)={dim_A2}=C(5,3)+2*{c3}={10+2*c3}, "
              f"dim(O_2)={dim_O2}, #constr={n_constr}, "
              f"dim_O2 = dim_A2 - #constr = {dim_A2 - n_constr}")


# ============================================================
# PART 4: Key formula test: rank(d_2) = n - 1 - beta_1?
# ============================================================
print(f"\n\n{'='*70}")
print("rank(d_2) FORMULA")
print("=" * 70)

# rank(d_2|Omega_2) = dim(Omega_2) - ker(d_2)
# If beta_2=0: rank(d_2) = dim(O_2) - im(d_3)
# Also: rank(d_2) = im(d_2|O_2) in O_1 = dim(im d_2)
# And beta_1 = dim(O_1) - rank(d_1) - im(d_2)
# Since dim(O_1) = n(n-1)/2, rank(d_1) = n-1:
# beta_1 = n(n-1)/2 - (n-1) - im_d2 = (n-1)(n-2)/2 - im_d2
# => im_d2 = (n-1)(n-2)/2 - beta_1

print(f"\nn=5: rank(d_2) = (n-1)(n-2)/2 - beta_1 = 6 - beta_1")
for bits in range(min(100, total)):
    A = build_adj(5, bits)
    betti = path_betti_numbers(A, 5, max_dim=2)
    b1 = int(betti[1]) if len(betti) > 1 else 0

    paths2 = enumerate_allowed_paths(A, 5, 2)
    paths1 = enumerate_allowed_paths(A, 5, 1)
    omega2 = compute_omega_basis(A, 5, 2, paths2, paths1)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

    D2 = build_full_boundary_matrix(paths2, paths1)
    D2_omega = D2 @ omega2
    S2 = np.linalg.svd(D2_omega, compute_uv=False)
    rank_d2 = sum(s > 1e-8 for s in S2)

    predicted_rank = 6 - b1
    if rank_d2 != predicted_rank:
        print(f"  MISMATCH T#{bits}: rank(d_2)={rank_d2}, predicted={predicted_rank}, beta_1={b1}")
        break
else:
    print(f"  First 100: rank(d_2) = 6 - beta_1 CONFIRMED")

# Check all
all_match = True
for bits in range(total):
    A = build_adj(5, bits)
    betti = path_betti_numbers(A, 5, max_dim=2)
    b1 = int(betti[1]) if len(betti) > 1 else 0

    paths2 = enumerate_allowed_paths(A, 5, 2)
    paths1 = enumerate_allowed_paths(A, 5, 1)
    omega2 = compute_omega_basis(A, 5, 2, paths2, paths1)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

    D2 = build_full_boundary_matrix(paths2, paths1)
    D2_omega = D2 @ omega2
    S2 = np.linalg.svd(D2_omega, compute_uv=False)
    rank_d2 = sum(s > 1e-8 for s in S2)

    predicted = (5-1)*(5-2)//2 - b1
    if rank_d2 != predicted:
        all_match = False
        break

print(f"  ALL n=5: rank(d_2) = (n-1)(n-2)/2 - beta_1? {all_match}")


# ============================================================
# PART 5: So beta_2=0 iff ker(d_2) = im(d_3)
# iff dim(Omega_2) - rank(d_2) = rank(d_3|Omega_3)
# iff dim(Omega_2) - ((n-1)(n-2)/2 - beta_1) = rank(d_3|Omega_3)
# ============================================================
print(f"\n\n{'='*70}")
print("BETA_2 = 0 EQUIVALENT CONDITION")
print("=" * 70)

# beta_2 = 0 iff:
# dim(Omega_2) - (n-1)(n-2)/2 + beta_1 = rank(d_3|Omega_3)
#
# LHS = dim(Omega_2) - C(n-1,2) + beta_1
# Using dim(Omega_2) = dim(A_2) - #constraints = C(n,3) + 2*c3 - #constraints:
# LHS = C(n,3) + 2*c3 - #constraints - C(n-1,2) + beta_1

# For n=5: C(5,3)=10, C(4,2)=6
# LHS = 10 + 2*c3 - #constraints - 6 + beta_1 = 4 + 2*c3 - #constraints + beta_1

print(f"\nn=5: beta_2=0 iff rank(d_3|O_3) = 4 + 2*c3 - #constraints + beta_1")
for scores in sorted(score_data.keys()):
    entries = score_data[scores]
    e = entries[0]
    c3 = e[1]
    ker_d2 = e[5]
    im_d3 = e[7]
    n_constr = e[8]

    betti = path_betti_numbers(build_adj(5, 0), 5, max_dim=2)
    # Need beta_1 for this score class
    A_ex = None
    for bits in range(total):
        A = build_adj(5, bits)
        if tuple(sorted([sum(row) for row in A])) == scores:
            A_ex = A
            break
    if A_ex:
        betti = path_betti_numbers(A_ex, 5, max_dim=2)
        b1 = int(betti[1]) if len(betti) > 1 else 0
        lhs = 4 + 2*c3 - n_constr + b1
        print(f"  scores={scores}: c3={c3}, #constr={n_constr}, beta_1={b1}, "
              f"predicted_im(d3)={lhs}, actual_im(d3)={im_d3}, match={lhs==im_d3}")


# ============================================================
# PART 6: Direct approach - does d_2 kill ALL of Omega_2 at p=2
# in a way that differs from p=3?
# ============================================================
print(f"\n\n{'='*70}")
print("RATIO: rank(d_p) / dim(Omega_p) FOR DIFFERENT p")
print("=" * 70)

# If rank(d_p)/dim(Omega_p) is consistently high at p=2,
# it means d_2 is "almost injective" on Omega_2

n = 5
ratios_by_p = defaultdict(list)

for bits in range(total):
    A = build_adj(n, bits)
    allowed = {}
    for p in range(-1, n):
        allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

    for p in range(1, n-1):
        omega_p = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
        dim_Op = omega_p.shape[1] if omega_p.ndim == 2 else 0
        if dim_Op == 0:
            continue
        D = build_full_boundary_matrix(allowed[p], allowed[p-1])
        D_omega = D @ omega_p
        S = np.linalg.svd(D_omega, compute_uv=False)
        rank = sum(s > 1e-8 for s in S)
        ratios_by_p[p].append(rank / dim_Op)

print(f"\nn=5: rank(d_p)/dim(Omega_p) statistics:")
for p in sorted(ratios_by_p.keys()):
    vals = ratios_by_p[p]
    print(f"  p={p}: min={min(vals):.4f}, max={max(vals):.4f}, "
          f"mean={np.mean(vals):.4f}, count={len(vals)}")
    # Distribution of ratio
    from collections import Counter
    ratio_counts = Counter(round(v, 3) for v in vals)
    for r, c in sorted(ratio_counts.items()):
        print(f"    ratio={r:.3f}: {c}")


print("\n\nDone.")
