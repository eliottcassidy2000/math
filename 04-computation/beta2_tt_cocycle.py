#!/usr/bin/env python3
r"""beta2_tt_cocycle.py - TT cocycle characterization of b1

b1(T) = dim of the space of 1-cycles z such that:
  z(a,b) + z(b,c) - z(a,c) = 0 for all TT triples (a,b,c)

This is z in ker(M_TT^T) intersect Z_1.

The condition z(a,b) + z(b,c) = z(a,c) whenever a->b->c and a->c
means z is "additive on transitive paths" — it's a quasimorphism-like condition.

If we think of z as a function on directed edges, the condition says:
for any transitive triple, the edge-value is determined by the
"shortcut" edge value and the "detour" edge values.

KEY QUESTION: For which tournaments does a nonzero TT-cocycle exist?
And can we prove the cocycle space has dimension at most 1?

APPROACH: Define the "TT graph" where vertices are edges of T, and
the constraint z(a,b) + z(b,c) = z(a,c) defines a linear system.
The solution space dimension is b1.

We need: the constraint system has corank at most 1 in Z_1.
Since Z_1 has dimension (n-1)(n-2)/2, we need at least
(n-1)(n-2)/2 - 1 independent constraints.

Each TT triple gives one constraint. We have C(n,3) - c_3 TT triples.
At n=5: C(5,3) - c_3 = 10 - c_3. Need at least 5 independent.
The 10 - c_3 constraints are not all independent (they satisfy relations).

What if we think of this as: the TT constraints define a "partial
additive structure" on the tournament edges?

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random, time
import numpy as np
from collections import Counter
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


def is_sc(A, n):
    if n <= 1:
        return True
    visited = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in range(n):
            if A[u][v] and v not in visited:
                visited.add(v)
                stack.append(v)
    if len(visited) < n:
        return False
    visited = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in range(n):
            if A[v][u] and v not in visited:
                visited.add(v)
                stack.append(v)
    return len(visited) == n


def count_c3(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    c3 += 1
                if A[i][k] and A[k][j] and A[j][i]:
                    c3 += 1
    return c3


def compute_tt_cocycle(A, n):
    """Compute the TT-cocycle space (= cokernel of TT boundary in Z_1).
    Returns: (dim_cocycle, cocycle_vectors_in_edge_coords, dim_Z1)
    """
    edges = []
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                edges.append((i, j))
    edge_idx = {e: i for i, e in enumerate(edges)}
    num_edges = len(edges)

    # d_1 matrix
    D1 = np.zeros((n, num_edges))
    for idx, (a, b) in enumerate(edges):
        D1[b, idx] += 1
        D1[a, idx] -= 1

    # Z_1 basis
    U, S, Vt = np.linalg.svd(D1, full_matrices=True)
    rk_d1 = int(sum(s > 1e-8 for s in S))
    Z1_basis = Vt[rk_d1:]  # shape (dim_Z1, num_edges)
    dim_Z1 = Z1_basis.shape[0]

    # TT triples
    tt_triples = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]:
                continue
            for c in range(n):
                if c == a or c == b or not A[b][c] or not A[a][c]:
                    continue
                tt_triples.append((a, b, c))

    if not tt_triples:
        return dim_Z1, Z1_basis, dim_Z1, edges

    # TT boundary matrix in edge coordinates
    M_TT = np.zeros((num_edges, len(tt_triples)))
    for j, (a, b, c) in enumerate(tt_triples):
        M_TT[edge_idx[(b, c)], j] += 1
        M_TT[edge_idx[(a, c)], j] -= 1
        M_TT[edge_idx[(a, b)], j] += 1

    # Project into Z_1
    M_TT_Z1 = Z1_basis @ M_TT  # shape (dim_Z1, num_tt)

    # Cocycle space = left null space of M_TT_Z1^T = cokernel of M_TT_Z1
    U2, S2, Vt2 = np.linalg.svd(M_TT_Z1, full_matrices=True)
    rk = int(sum(s > 1e-8 for s in S2))
    dim_cocycle = dim_Z1 - rk

    if dim_cocycle > 0:
        cocycles_Z1 = U2[:, rk:]  # shape (dim_Z1, dim_cocycle)
        cocycles_edge = Z1_basis.T @ cocycles_Z1  # shape (num_edges, dim_cocycle)
    else:
        cocycles_edge = np.zeros((num_edges, 0))

    return dim_cocycle, cocycles_edge, dim_Z1, edges


# ============================================================
# Part 1: Characterize the TT-cocycle at n=4,5
# ============================================================
print("=" * 70)
print("TT-COCYCLE CHARACTERIZATION")
print("=" * 70)

# At n=4, b1=1 tournaments are exactly the score (1,1,2,2) ones.
# These are the 4-cycles with one "diagonal".
# The cocycle should have a specific structure.

n = 4
print(f"\nn={n}:")
for bits in range(2 ** (n*(n-1)//2)):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    dim_coc, coc_edge, dim_Z1, edges = compute_tt_cocycle(A, n)
    if dim_coc == 0:
        continue

    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    c3 = count_c3(A, n)

    # Analyze cocycle
    edge_vals = {}
    for i in range(len(edges)):
        if abs(coc_edge[i, 0]) > 1e-8:
            edge_vals[edges[i]] = round(float(coc_edge[i, 0]), 4)

    print(f"  bits={bits}, scores={scores}, c3={c3}, dim_coc={dim_coc}")
    # Normalize so max coefficient is 1
    max_val = max(abs(v) for v in edge_vals.values())
    norm_vals = {e: round(v/max_val, 4) for e, v in edge_vals.items()}
    print(f"    Cocycle (normalized): {norm_vals}")

    # Check: which edges are in the 3-cycles?
    cycle_edges = set()
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    cycle_edges |= {(i,j), (j,k), (k,i)}
                if A[i][k] and A[k][j] and A[j][i]:
                    cycle_edges |= {(i,k), (k,j), (j,i)}
    non_cycle_edges = set(edges) - cycle_edges
    print(f"    3-cycle edges: {cycle_edges}")
    print(f"    Non-cycle edges: {non_cycle_edges}")

    if bits > 30:  # Just show a few
        break


# ============================================================
# Part 2: At n=5, cocycle structure
# ============================================================
print(f"\n{'=' * 70}")
print(f"n=5: COCYCLE STRUCTURE FOR b1=1")
print("=" * 70)

n = 5
count = 0
cocycle_support_dist = Counter()
cocycle_signs = Counter()

for bits in range(2 ** (n*(n-1)//2)):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    dim_coc, coc_edge, dim_Z1, edges = compute_tt_cocycle(A, n)
    if dim_coc == 0:
        continue

    # Support size
    support = sum(abs(coc_edge[i, 0]) > 1e-8 for i in range(len(edges)))
    cocycle_support_dist[support] += 1

    # Sign pattern (how many positive vs negative)
    pos = sum(coc_edge[i, 0] > 1e-8 for i in range(len(edges)))
    neg = sum(coc_edge[i, 0] < -1e-8 for i in range(len(edges)))
    cocycle_signs[(pos, neg)] += 1

    count += 1

    if count <= 5:
        scores = tuple(sorted(sum(A[i]) for i in range(n)))
        c3 = count_c3(A, n)
        edge_vals = {}
        for i in range(len(edges)):
            if abs(coc_edge[i, 0]) > 1e-8:
                edge_vals[edges[i]] = round(float(coc_edge[i, 0]), 4)
        max_val = max(abs(v) for v in edge_vals.values())
        norm_vals = {e: round(v/max_val, 4) for e, v in edge_vals.items()}
        print(f"\n  bits={bits}, scores={scores}, c3={c3}")
        print(f"    Cocycle (normalized, {len(norm_vals)} edges): {norm_vals}")

print(f"\nn=5: {count} tournaments with b1=1")
print(f"  Support size distribution: {dict(sorted(cocycle_support_dist.items()))}")
print(f"  Sign pattern (pos, neg): {dict(sorted(cocycle_signs.items()))}")


# ============================================================
# Part 3: Key insight — vertex potentials?
# ============================================================
print(f"\n{'=' * 70}")
print("VERTEX POTENTIAL STRUCTURE")
print("=" * 70)

# The cocycle z satisfies z(a,b) + z(b,c) = z(a,c) for TT triples.
# This is ALMOST saying z(a,b) = phi(b) - phi(a) for some potential phi.
# But that would make z a coboundary, hence in im(d_1).
# Since z is in Z_1, it has zero boundary already.
# A coboundary z(a,b) = phi(b) - phi(a) IS in Z_1 and IS in im(d_1).
# But z is in ker(d_1) and NOT in im(d_2), so z is NOT a coboundary.
#
# Wait: Z_1 = ker(d_1). The 1-cocycles (= elements of ker(d_1^T))
# are NOT the same as 1-cycles (= elements of ker(d_1)).
#
# Let me reconsider. The TT constraint z(a,b) + z(b,c) = z(a,c) is:
# z(a,c) - z(a,b) - z(b,c) = 0, which is the same as:
# <z, d_2(a,b,c)> = 0 where d_2(a,b,c) = (a,b) + (b,c) - (a,c).
#
# So z is in the ORTHOGONAL COMPLEMENT of im(d_2) within Z_1.
# Since Z_1/im(d_2) = H_1, and b1 = dim H_1, the cocycle space = H_1.
#
# The interesting question is: what does the cocycle LOOK like?
# If z(a,b) + z(b,c) = z(a,c) for all TT triples, does this force
# some kind of "potential-like" structure on the 3-cycle edges?

# Let me check: for a 3-cycle a->b->c->a, is z(a,b) + z(b,c) + z(c,a) = 0?
# From the cycle condition: d_1(z) = 0, so sum of z over any directed cycle = 0.
# Wait: d_1(z) = 0 means sum_{edges leaving v} z(v,w) - sum_{edges entering v} z(w,v) = 0.
# This is NOT the same as cycle-sum = 0.
#
# Actually for z in Z_1:
# z is a 1-chain alpha in Om_1 = R^{edges}.
# d_1(alpha) is the 0-chain. d_1(alpha)(v) = sum_w alpha(w,v) - sum_w alpha(v,w)?
# No: d_1 sends edge (a,b) to vertex b minus vertex a.
# d_1(alpha) = sum_e alpha(e) * d_1(e) = sum_{a->b} alpha(a,b) * [(b) - (a)].
# So d_1(alpha)(v) = sum_{a: a->v} alpha(a,v) - sum_{b: v->b} alpha(v,b).
# d_1(alpha) = 0 means this is zero for all v.
# This is the "flow conservation" / "balanced flow" condition.
# It does NOT directly say cycle sums are zero.
#
# But for a 3-cycle a->b->c->a:
# The formal sum (a,b) + (b,c) + (c,a) in R^{edges} has:
# d_1[(a,b)+(b,c)+(c,a)] = (b-a) + (c-b) + (a-c) = 0.
# So the cycle sum IS a 1-cycle. And z is a 1-cycle too.
# But <z, cycle> is NOT necessarily 0.

# Let me compute: for b1=1 tournaments, what are the cycle-sums of the cocycle?
n = 5
print(f"\nn=5: cycle sums of TT-cocycle")
count_printed = 0
for bits in range(2 ** (n*(n-1)//2)):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    dim_coc, coc_edge, dim_Z1, edges = compute_tt_cocycle(A, n)
    if dim_coc == 0:
        continue

    edge_vals = {}
    for i in range(len(edges)):
        edge_vals[edges[i]] = float(coc_edge[i, 0])

    # Compute cycle sums for all 3-cycles
    cycle_sums = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    s = edge_vals.get((i,j), 0) + edge_vals.get((j,k), 0) + edge_vals.get((k,i), 0)
                    cycle_sums.append(((i,j,k), round(s, 6)))
                if A[i][k] and A[k][j] and A[j][i]:
                    s = edge_vals.get((i,k), 0) + edge_vals.get((k,j), 0) + edge_vals.get((j,i), 0)
                    cycle_sums.append(((i,k,j), round(s, 6)))

    count_printed += 1
    if count_printed <= 5:
        scores = tuple(sorted(sum(A[i]) for i in range(n)))
        c3 = count_c3(A, n)
        print(f"\n  bits={bits}, scores={scores}, c3={c3}")
        print(f"    3-cycle sums: {cycle_sums}")
        # Normalize
        max_sum = max(abs(s) for _, s in cycle_sums) if cycle_sums else 1
        if max_sum > 1e-8:
            norm_sums = [(c, round(s/max_sum, 4)) for c, s in cycle_sums]
            print(f"    Normalized: {norm_sums}")

    if count_printed > 5:
        break


# ============================================================
# Part 4: Is the cocycle always supported on cycle edges?
# ============================================================
print(f"\n{'=' * 70}")
print("COCYCLE SUPPORT vs CYCLE EDGES")
print("=" * 70)

n = 5
always_on_cycles = 0
not_always = 0

for bits in range(2 ** (n*(n-1)//2)):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    dim_coc, coc_edge, dim_Z1, edges = compute_tt_cocycle(A, n)
    if dim_coc == 0:
        continue

    # Find cycle edges
    cycle_edges = set()
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    cycle_edges |= {(i,j), (j,k), (k,i)}
                if A[i][k] and A[k][j] and A[j][i]:
                    cycle_edges |= {(i,k), (k,j), (j,i)}

    # Check support
    support = set()
    for i in range(len(edges)):
        if abs(coc_edge[i, 0]) > 1e-8:
            support.add(edges[i])

    if support <= cycle_edges:
        always_on_cycles += 1
    else:
        not_always += 1

print(f"\nn=5: Cocycle supported on 3-cycle edges: {always_on_cycles}/{always_on_cycles+not_always}")
print(f"  NOT fully on cycle edges: {not_always}")


# ============================================================
# Part 5: Check b1 <= 1 for larger n (already done but reconfirm)
# ============================================================
print(f"\n{'=' * 70}")
print("b1 <= 1 VERIFICATION (EXTENDED)")
print("=" * 70)

for n in [7, 8, 9, 10, 12]:
    random.seed(42)
    trials = max(10, min(200, 2000 // n))
    max_b1 = 0

    for trial in range(trials):
        A = random_tournament(n)
        dim_coc, _, _, _ = compute_tt_cocycle(A, n)
        max_b1 = max(max_b1, dim_coc)

    print(f"  n={n}: max b1 = {max_b1} ({trials} trials)")


print("\n\nDone.")
