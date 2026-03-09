#!/usr/bin/env python3
"""beta2_omega2_arcflip.py - Understand dim(Omega_2) change under arc flip

KEY: dim(Om_2) != #allowed_2_paths for tournaments.
Om_2 = ker(d_1) restricted to allowed 2-paths.
d_1(a,b,c) = (a,c). So (a,b,c) in Om_2 iff a->b->c allowed AND
the formal combination cancels under d_1.

Actually, let me re-read the definition. From Grigoryan et al:
A_p = allowed p-paths (those not hitting a vertex twice).
Omega_p = {alpha in R^{A_p} : d_i(alpha) = 0 for 0 < i < p}

For p=2: Omega_2 = {alpha in R^{A_2} : d_1(alpha) = 0}
d_1(a,b,c) = (a,c). So d_1(alpha) = 0 means:
sum_{b: (a,b,c) in A_2} alpha(a,b,c) = 0 for every (a,c) in A_1.

Wait, that's not right either. d_1 is the "inner" face map.
For a 2-path (v_0, v_1, v_2):
  d_0(v0,v1,v2) = (v1,v2)
  d_1(v0,v1,v2) = (v0,v2)
  d_2(v0,v1,v2) = (v0,v1)

Omega_p = intersection of ker(d_i) for 0 < i < p.
For p=2, this means d_1(alpha) = 0.
For the coefficient of each 1-path (a,c), d_1(alpha)(a,c) = sum of alpha(a,b,c) over all b s.t. (a,b,c) is allowed.
So the condition is: for each pair (a,c) with a->c, sum_{b: a->b->c} alpha(a,b,c) = 0.

This means: if (a,c) has exactly one mediator b, then alpha(a,b,c) = 0.
If (a,c) has k mediators b1,...,bk, then sum alpha(a,bi,c) = 0 (k-1 free parameters).

So dim(Om_2) = |A_2| - |{(a,c): a->c with at least one mediator b}|
             = |A_2| - |{(a,c): a->c AND exists b with a->b, b->c, b!=a, b!=c}|

For tournaments: a->c always. So (a,c) has mediators iff exists b with a->b and b->c.
This is equivalent to (a,c) being a TT pair for some 2-path (a,b,c).

Actually dim(Om_2) = |A_2| - #{TT pairs}. Wait, let me reconsider:
Each TT pair (a,c) with k mediators gives k-1 constraints.
So dim(Om_2) = |A_2| - sum_{(a,c) TT} 1 (one constraint per pair).
Wait no, the constraint is one LINEAR equation per (a,c) pair that has mediators.
Each such pair removes 1 dimension. So:
dim(Om_2) = |A_2| - #{pairs (a,c) that have at least one mediator}

For a tournament: every edge a->c is present. The pair (a,c) has a mediator iff
there exists b with a->b and b->c (i.e., a 2-path (a,b,c)).
If no such b exists, then there's no 2-path with first/last vertex (a,c), so no constraint.

Actually, this double-counts. Let me think again...

Each 2-path (a,b,c) contributes to the d_1 constraint for pair (a,c).
The number of constraints is the number of pairs (a,c) such that the number of
2-paths (a,*,c) is >= 1. Each such pair contributes one constraint.

So dim(Om_2) = |A_2| - |{(a,c) : exists b with (a,b,c) in A_2}|

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


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def compute_omega2_dim_formula(A, n):
    """Compute dim(Om_2) using the formula: |A_2| - #{constrained pairs}."""
    # Enumerate all allowed 2-paths
    paths2 = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]:
                continue
            for c in range(n):
                if c == a or c == b or not A[b][c]:
                    continue
                paths2.append((a, b, c))

    # For each pair (a,c), count mediators
    pair_mediators = defaultdict(int)
    for a, b, c in paths2:
        pair_mediators[(a, c)] += 1

    # Number of constrained pairs = those with >= 1 mediator
    num_constrained = len(pair_mediators)

    return len(paths2), num_constrained, len(paths2) - num_constrained


def compute_omega2_dim_library(A, n):
    """Compute dim(Om_2) using the library."""
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths2 = enumerate_allowed_paths(A, n, 2)
    if not paths2 or not paths1:
        return 0
    paths0 = [(i,) for i in range(n)]
    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    return omega2.shape[1] if omega2.ndim == 2 else 0


def flip_arc(A, n, u, v):
    B = [row[:] for row in A]
    B[u][v] = 0
    B[v][u] = 1
    return B


# ============================================================
# Part 1: Verify the formula dim(Om_2) = |A_2| - #constrained_pairs
# ============================================================
print("=" * 70)
print("VERIFY: dim(Om_2) = |A_2| - #constrained_pairs")
print("=" * 70)

for n in [4, 5, 6]:
    ok = 0
    fail = 0
    random.seed(42)
    trials = {4: 64, 5: 200, 6: 100}[n]
    for trial in range(trials):
        A = random_tournament(n)
        p2, nc, formula_dim = compute_omega2_dim_formula(A, n)
        lib_dim = compute_omega2_dim_library(A, n)
        if formula_dim == lib_dim:
            ok += 1
        else:
            fail += 1
            if fail <= 3:
                print(f"  n={n}: formula={formula_dim} ({p2}-{nc}), lib={lib_dim}")
    print(f"\nn={n}: formula matches: {ok}/{ok+fail}")


# ============================================================
# Part 2: Change in dim(Om_2) under arc flip using formula
# ============================================================
print(f"\n{'=' * 70}")
print("D(dim Om_2) UNDER ARC FLIP — FORMULA")
print("=" * 70)

# When we flip u->v to v->u:
# A_2 changes: paths using u->v disappear, paths using v->u appear.
# Constrained pairs change: some (a,c) gain/lose mediators.
#
# Gained paths: (w, v, u) for w->v, and (v, u, w) for u->w
#   Also paths not using the flipped edge but newly valid? No — only edge u->v changes.
#   Other edges unchanged, so only paths THROUGH u->v or v->u change.
# Lost paths: (w, u, v) for w->u, and (u, v, w) for v->w
#
# But wait: there are ALSO 2-paths that use u and v but not the edge between them:
#   (u, w, v) or (v, w, u) for some w.
#   These DON'T change because w->u/v arcs are unchanged.
# And (a, u, c) where c != v: unchanged (u's edges to non-v unchanged).
# And (a, v, c) where a != u: unchanged.
#
# So the ONLY changes are paths with the flipped edge as one of its two edges.

for n in [5, 6]:
    random.seed(42)
    trials = {5: 300, 6: 100}[n]
    formula_data = Counter()

    for trial in range(trials):
        A = random_tournament(n)
        edges = [(i, j) for i in range(n) for j in range(n) if A[i][j]]
        u, v = edges[random.randint(0, len(edges) - 1)]
        B = flip_arc(A, n, u, v)

        # ABCD
        a_set = [w for w in range(n) if w != u and w != v and A[u][w] and A[w][v]]
        b_set = [w for w in range(n) if w != u and w != v and A[w][u] and A[v][w]]
        c_set = [w for w in range(n) if w != u and w != v and A[u][w] and A[v][w]]
        d_set = [w for w in range(n) if w != u and w != v and A[w][u] and A[w][v]]
        a, b, c, d = len(a_set), len(b_set), len(c_set), len(d_set)

        p2_A, nc_A, dim_A = compute_omega2_dim_formula(A, n)
        p2_B, nc_B, dim_B = compute_omega2_dim_formula(B, n)

        d_p2 = p2_B - p2_A
        d_nc = nc_B - nc_A
        d_dim = dim_B - dim_A

        lib_dim_A = compute_omega2_dim_library(A, n)
        lib_dim_B = compute_omega2_dim_library(B, n)
        d_lib = lib_dim_B - lib_dim_A

        formula_data[(a, b, c, d, d_p2, d_nc, d_dim, d_lib)] += 1

    print(f"\nn={n}: ABCD + path/constraint changes under arc flip")
    print(f"{'(a,b,c,d)':>12} {'Dp2':>5} {'Dnc':>5} {'Ddim':>5} {'Dlib':>5} {'count':>6}")
    for (a, b, c, d, dp2, dnc, ddim, dlib), cnt in sorted(formula_data.items()):
        marker = " <-- MISMATCH" if ddim != dlib else ""
        print(f"  ({a},{b},{c},{d})    {dp2:>5} {dnc:>5} {ddim:>5} {dlib:>5}   {cnt:>5}{marker}")


# ============================================================
# Part 3: Verify D(p2) = 2(a-b) AND find D(nc) formula
# ============================================================
print(f"\n{'=' * 70}")
print("D(nc) AS FUNCTION OF a,b,c,d")
print("=" * 70)

# We know D(p2) = 2(a-b). What about D(nc)?
# nc = number of pairs (x,y) with x->y that have at least one 2-path mediator.
#
# Under arc flip u->v to v->u:
# LOST mediators: pair (w,v) had mediator u (w->u->v). After flip, u no longer beats v.
#   But w->u still holds. v->u now. So (w,v) path (w,u,v) is gone.
#   Is there another mediator for (w,v)? Only if exists z!=u with w->z and z->v.
#   For w in B or D: w->u was true.
# GAINED mediators: pair (w,u) now has mediator v (w->v->u). Before flip, v->u didn't hold.
#   For w in A or D: w->v and v->u. So (w,u) gains mediator.
# Also: pair (v,w): (v,u,w) with u->w. Before: u->v, so (v,u,w) not a path (v doesn't beat u).
#   After: v->u is true. (v,u,w) with u->w. So pair (v,w) gains mediator u.
#   For w in A or C: u->w.
# Also: pair (u,w): (u,v,w) with v->w. Before: u->v, so (u,v,w) was valid.
#   After: u no longer beats v. So (u,v,w) disappears.
#   For w in B or C: v->w. Pair (u,w) loses mediator v.

# This is getting complex. Let me just compute.

for n in [5]:
    random.seed(42)
    dnc_by_abcd = defaultdict(list)
    for trial in range(1000):
        A = random_tournament(n)
        edges = [(i, j) for i in range(n) for j in range(n) if A[i][j]]
        u, v = edges[random.randint(0, len(edges) - 1)]
        B_mat = flip_arc(A, n, u, v)

        a_set = [w for w in range(n) if w != u and w != v and A[u][w] and A[w][v]]
        b_set = [w for w in range(n) if w != u and w != v and A[w][u] and A[v][w]]
        c_set = [w for w in range(n) if w != u and w != v and A[u][w] and A[v][w]]
        d_set = [w for w in range(n) if w != u and w != v and A[w][u] and A[w][v]]
        a, b, c, d = len(a_set), len(b_set), len(c_set), len(d_set)

        _, nc_A, _ = compute_omega2_dim_formula(A, n)
        _, nc_B, _ = compute_omega2_dim_formula(B_mat, n)
        dnc = nc_B - nc_A

        dnc_by_abcd[(a, b, c, d)].append(dnc)

    print(f"\nn={n}: D(nc) by ABCD partition")
    for key in sorted(dnc_by_abcd.keys()):
        vals = dnc_by_abcd[key]
        unique = sorted(set(vals))
        dist = Counter(vals)
        print(f"  {key}: values={unique}, counts={dict(sorted(dist.items()))}")

    print(f"\n  D(nc) is NOT constant per ABCD partition!")
    print(f"  But D(dim Om_2) = D(p2) - D(nc) = 2(a-b) - D(nc)")
    print(f"  So D(dim Om_2) depends on more than just (a,b,c,d).")


# ============================================================
# Part 4: What DOES D(dim Om_2) depend on?
# ============================================================
print(f"\n{'=' * 70}")
print("D(dim Om_2) DECOMPOSITION")
print("=" * 70)

# Since D(nc) is not determined by (a,b,c,d), maybe D(dim Om_2) isn't either.
# But beta_2 IS invariant. So D(rk d3) must track D(dim Z_2) perfectly.

# Let me compute D(dim Om_2) by ABCD to see if IT is also variable.
for n in [5]:
    random.seed(42)
    ddim_by_abcd = defaultdict(list)
    for trial in range(1000):
        A = random_tournament(n)
        edges = [(i, j) for i in range(n) for j in range(n) if A[i][j]]
        u, v = edges[random.randint(0, len(edges) - 1)]
        B_mat = flip_arc(A, n, u, v)

        a_set = [w for w in range(n) if w != u and w != v and A[u][w] and A[w][v]]
        b_set = [w for w in range(n) if w != u and w != v and A[w][u] and A[v][w]]
        c_set = [w for w in range(n) if w != u and w != v and A[u][w] and A[v][w]]
        d_set = [w for w in range(n) if w != u and w != v and A[w][u] and A[w][v]]
        a, b, c, d = len(a_set), len(b_set), len(c_set), len(d_set)

        dim_A = compute_omega2_dim_library(A, n)
        dim_B = compute_omega2_dim_library(B_mat, n)
        ddim = dim_B - dim_A

        ddim_by_abcd[(a, b, c, d)].append(ddim)

    print(f"\nn={n}: D(dim Om_2) by ABCD partition")
    for key in sorted(ddim_by_abcd.keys()):
        vals = ddim_by_abcd[key]
        unique = sorted(set(vals))
        dist = Counter(vals)
        print(f"  {key}: values={unique}, counts={dict(sorted(dist.items()))}")


print("\n\nDone.")
