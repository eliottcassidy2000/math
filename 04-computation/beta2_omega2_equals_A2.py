#!/usr/bin/env python3
"""
beta2_omega2_equals_A2.py - Check if Omega_2 = A_2 for tournaments

For a tournament T, every pair of vertices has an arc, so A_1 = all arcs
(complete set of directed edges). A 2-path (a,b,c) with a->b, b->c is in
Omega_2 iff all faces are in A_1:
  face (b,c) in A_1: always yes (b->c)
  face (a,c) in A_1: yes iff there's an arc a->c or c->a. In a tournament, YES.
  face (a,b) in A_1: always yes (a->b)

So the faces of (a,b,c) are (b,c), (a,c), (a,b) — all single arcs.
In a tournament, every arc exists in one direction, so (a,c) IS in A_1
(as either a->c or c->a). But wait: A_1 consists of ALLOWED 1-paths,
i.e., arcs that actually exist. So (a,c) is in A_1 iff either a->c or c->a
is an arc. In a tournament, exactly one of these holds. But A_1 is the
set of 1-paths (x,y) with x->y. So (a,c) is in A_1 iff a->c.

Wait: the face map for a 2-path (a,b,c) is:
  d_2(a,b,c) = (b,c) - (a,c) + (a,b)

For (a,b,c) to be in Omega_2, we need d_2(a,b,c) to be in span(A_1).
The term (a,c) appears with coefficient -1. If (a,c) is NOT in A_1
(i.e., c->a instead of a->c), then (a,c) is not a basis element of
the A_1 space, and the boundary would have a component outside A_1.

But actually: Omega_p = {u in A_p : d(u) is in span of A_{p-1}}.
The boundary d_2(a,b,c) = (b,c) - (a,c) + (a,b).
This is always a vector in the span of all possible 1-paths.
But for the GLMY complex, Omega_2 requires that d_2(a,b,c) is in
the span of ALLOWED 1-paths A_1.

If a->c, then (a,c) in A_1, and d_2(a,b,c) = (b,c) - (a,c) + (a,b) in span(A_1). OK.
If c->a, then (a,c) NOT in A_1, but (c,a) IS in A_1. The term -(a,c) is
NOT in span(A_1). So (a,b,c) is NOT in Omega_2.

So Omega_2 = {(a,b,c) : a->b, b->c, AND a->c} = TT triples!
(Transitive triples, where a->b->c and a->c.)

Wait, but from the data above, dim(Omega_2) for the transitive tournament on n=5
is 10 = C(5,3), which equals |A_2| = 10 (all 3-subsets are TT in transitive tournament).

But for non-transitive tournaments, Omega_2 < A_2 because some triples are 3-cycles
(a->b->c->a), which are NOT TT.

Hmm, but the path_homology_v2 code computes Omega_2 differently. Let me verify.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis
)
sys.stdout = _saved

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_TT(A, n):
    """Count transitive triples: (a,b,c) with a->b, b->c, a->c."""
    count = 0
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b or not A[b][c]: continue
                if A[a][c]:  # Transitive: a->c
                    count += 1
    return count

n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

print("=" * 70)
print(f"OMEGA_2 vs A_2 vs TT at n={n}")
print("=" * 70)

matches_TT = 0
matches_A2 = 0
total_checked = 0

for bits in range(total):
    A = build_adj(n, bits)
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    allowed_1 = enumerate_allowed_paths(A, n, 1)
    omega2_basis = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
    dim_O2 = omega2_basis.shape[1] if omega2_basis.ndim == 2 else 0
    n_A2 = len(allowed_2)
    n_TT = count_TT(A, n)

    total_checked += 1

    if dim_O2 == n_TT:
        matches_TT += 1
    if dim_O2 == n_A2:
        matches_A2 += 1

    if total_checked <= 5:
        print(f"  bits={bits}: |A_2|={n_A2}, dim(O_2)={dim_O2}, |TT|={n_TT}")

print(f"\n  dim(O_2) == |TT|: {matches_TT}/{total} ({100*matches_TT/total:.1f}%)")
print(f"  dim(O_2) == |A_2|: {matches_A2}/{total} ({100*matches_A2/total:.1f}%)")

# Also check: is Omega_2 ALWAYS the TT subspace?
# That is, are the Omega_2 basis vectors exactly the TT 2-paths?
print(f"\n  Checking Omega_2 basis structure for first 20 tournaments:")
for bits in range(20):
    A = build_adj(n, bits)
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    allowed_1 = enumerate_allowed_paths(A, n, 1)
    omega2_basis = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
    dim_O2 = omega2_basis.shape[1] if omega2_basis.ndim == 2 else 0

    # Which allowed 2-paths are TT?
    tt_indices = []
    ntt_indices = []
    for i, p in enumerate(allowed_2):
        if A[p[0]][p[2]]:  # a->c (transitive)
            tt_indices.append(i)
        else:
            ntt_indices.append(i)

    # Is Omega_2 basis exactly the TT paths?
    # Check if omega2 basis has nonzero entries only at TT indices
    if dim_O2 > 0:
        ntt_contributions = 0
        for j in range(dim_O2):
            col = omega2_basis[:, j]
            for idx in ntt_indices:
                if abs(col[idx]) > 1e-8:
                    ntt_contributions += 1
                    break

        if ntt_contributions > 0:
            print(f"  bits={bits}: dim(O2)={dim_O2}, |TT|={len(tt_indices)}, |NTT|={len(ntt_indices)}, O2 has NTT components: {ntt_contributions}")

# So the question is: does Omega_2 contain ONLY TT paths, or also
# linear combinations involving NTT paths?
# From THM-101 notes: "dim(Omega_2) = |TT| + dim(NT cancellation space)"
# This means Omega_2 contains BOTH TT paths AND certain NTT combinations.

print(f"\n{'='*70}")
print("OMEGA_2 DIMENSION FORMULA")
print(f"{'='*70}")

# Check: dim(Omega_2) = |TT| always? Or dim(Omega_2) > |TT|?
for n in [5, 6]:
    n_arcs = n*(n-1)//2
    total = 1 << n_arcs
    O2_eq_TT = 0
    O2_gt_TT = 0
    max_diff = 0

    for bits in range(total):
        A = build_adj(n, bits)
        allowed_2 = enumerate_allowed_paths(A, n, 2)
        allowed_1 = enumerate_allowed_paths(A, n, 1)
        omega2_basis = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
        dim_O2 = omega2_basis.shape[1] if omega2_basis.ndim == 2 else 0
        n_TT = count_TT(A, n)

        if dim_O2 == n_TT:
            O2_eq_TT += 1
        elif dim_O2 > n_TT:
            O2_gt_TT += 1
            max_diff = max(max_diff, dim_O2 - n_TT)

    print(f"  n={n}: O2=TT: {O2_eq_TT}/{total}, O2>TT: {O2_gt_TT}/{total}, max_diff={max_diff}")

print("\nDone.")
