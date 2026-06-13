#!/usr/bin/env python3
"""
beta2_simp_hole_filling.py — How do DT_c + cancel fill the flag complex holes?

KEY DISCOVERY: H₂^simp(flag) ≠ 0 for 40/1024 n=5 tournaments (all c₃=4).
Yet path β₂ = 0 always. The difference comes from DT_c paths and cancel pairs.

This script identifies:
1. The simplicial 2-cycle that is NOT a simplicial boundary
2. The DT_c path or cancel pair that fills it in path homology
3. The algebraic mechanism by which DT_c "couples" TT and NT

Author: opus-2026-03-08-S49
"""
import sys
import numpy as np
from collections import defaultdict
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
print("SIMPLICIAL HOLE FILLING BY DT_c + CANCEL")
print("=" * 70)

n = 5
m = n*(n-1)//2

count = 0
for bits in range(1 << m):
    A = build_adj(n, bits)

    scores = tuple(sorted(sum(A[i][j] for j in range(n) if j!=i) for i in range(n)))
    c3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
             if A[i][j]+A[j][k]+A[k][i] in (0,3))
    if c3 != 4:
        continue

    # Check flag H₂
    # TT triples
    tt_triples = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c in (a,b) or not A[b][c]: continue
                if A[a][c]:
                    tt_triples.append((a, b, c))

    arcs = [(i, j) for i in range(n) for j in range(n) if i != j and A[i][j]]
    arc_idx = {a: i for i, a in enumerate(arcs)}

    d2_simp = np.zeros((len(arcs), len(tt_triples)))
    for j, (a, b, c) in enumerate(tt_triples):
        if (b, c) in arc_idx: d2_simp[arc_idx[(b, c)], j] += 1
        if (a, c) in arc_idx: d2_simp[arc_idx[(a, c)], j] -= 1
        if (a, b) in arc_idx: d2_simp[arc_idx[(a, b)], j] += 1

    rk_d2_simp = np.linalg.matrix_rank(d2_simp, tol=1e-8)
    z2_simp_dim = len(tt_triples) - rk_d2_simp

    # Transitive 4-subsets
    trans4 = []
    for quad in combinations(range(n), 4):
        # Check if transitive
        has_cycle = False
        for i, a in enumerate(quad):
            for j, b in enumerate(quad):
                if j <= i: continue
                for k, c in enumerate(quad):
                    if k <= j: continue
                    total = A[a][b] + A[b][c] + A[c][a]
                    if total == 0 or total == 3:
                        has_cycle = True
        if not has_cycle:
            # Find total order
            sc = [(sum(A[v][w] for w in quad if w != v), v) for v in quad]
            sc.sort(reverse=True)
            trans4.append(tuple(v for _, v in sc))

    tt_idx = {t: i for i, t in enumerate(tt_triples)}
    d3_simp = np.zeros((len(tt_triples), len(trans4)))
    for j, (a, b, c, d) in enumerate(trans4):
        faces = [(1, (b,c,d)), (-1, (a,c,d)), (1, (a,b,d)), (-1, (a,b,c))]
        for sign, face in faces:
            if face in tt_idx:
                d3_simp[tt_idx[face], j] += sign

    rk_d3_simp = np.linalg.matrix_rank(d3_simp, tol=1e-8)
    h2_simp = z2_simp_dim - rk_d3_simp

    if h2_simp == 0:
        continue

    count += 1
    if count > 3:
        continue

    print(f"\n{'='*60}")
    print(f"T#{bits} scores={scores}, c₃={c3}, H₂^simp={h2_simp}")
    print(f"{'='*60}")
    print(f"  TT triples: {len(tt_triples)}, trans4: {len(trans4)}")
    print(f"  rk(∂₂^simp)={rk_d2_simp}, Z₂^simp={z2_simp_dim}, rk(∂₃^simp)={rk_d3_simp}")

    # Find the simplicial cycle not in im(∂₃^simp)
    U, S, Vt = np.linalg.svd(d2_simp, full_matrices=True)
    z2_simp_basis = Vt[rk_d2_simp:].T  # columns are Z₂^simp basis

    # Find which Z₂^simp directions are NOT in im(∂₃^simp)
    if d3_simp.shape[1] > 0:
        # Project Z₂^simp basis into the range of d3_simp
        proj = z2_simp_basis.T @ d3_simp  # z2_simp_dim × n_trans4
        U2, S2, Vt2 = np.linalg.svd(proj, full_matrices=True)
        uncovered_idx = sum(s > 1e-8 for s in S2)
        uncovered = z2_simp_basis @ U2[:, uncovered_idx:]
    else:
        uncovered = z2_simp_basis

    # Show the uncovered simplicial cycle
    cycle_tt = uncovered[:, 0]
    print(f"\n  Simplicial cycle NOT in im(∂₃^simp):")
    for j, (a, b, c) in enumerate(tt_triples):
        if abs(cycle_tt[j]) > 1e-8:
            print(f"    {cycle_tt[j]:+.4f} * ({a},{b},{c})")

    # Now compute the PATH homology filling
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2)
    d3 = dim_om(om3)

    bd3 = build_full_boundary_matrix(ap3, ap2)
    bd3_om = bd3 @ om3
    coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]

    # Express the simplicial cycle in A₂ coordinates
    ap2_list = [tuple(x) for x in ap2]
    cycle_a2 = np.zeros(len(ap2_list))
    for j, (a, b, c) in enumerate(tt_triples):
        if abs(cycle_tt[j]) > 1e-8:
            if (a, b, c) in ap2_list:
                idx = ap2_list.index((a, b, c))
                cycle_a2[idx] = cycle_tt[j]

    # Express in Ω₂ coordinates
    cycle_om2 = np.linalg.lstsq(om2, cycle_a2, rcond=None)[0]

    # Find the Ω₃ element that fills it
    x, _, _, _ = np.linalg.lstsq(coords3, cycle_om2, rcond=None)
    err = np.max(np.abs(coords3 @ x - cycle_om2))
    print(f"\n  Path filling error: {err:.2e}")

    # What is the filling in A₃ coordinates?
    filling = om3 @ x
    ap3_list = [tuple(p) for p in ap3]

    print(f"  Filling element:")
    for j in range(len(ap3_list)):
        if abs(filling[j]) > 1e-8:
            a, b, c, d = ap3_list[j]
            is_dt = A[a][c] and A[b][d]
            ad_type = "a→d" if A[a][d] else "d→a"
            dt_label = f"DT({ad_type})" if is_dt else "non-DT"
            print(f"    {filling[j]:+.4f} * ({a},{b},{c},{d}) [{dt_label}]")

    # Classify the filling
    n_dtf = sum(1 for j in range(len(ap3_list))
                if abs(filling[j]) > 1e-8 and A[ap3_list[j][0]][ap3_list[j][2]]
                and A[ap3_list[j][1]][ap3_list[j][3]] and A[ap3_list[j][0]][ap3_list[j][3]])
    n_dtc = sum(1 for j in range(len(ap3_list))
                if abs(filling[j]) > 1e-8 and A[ap3_list[j][0]][ap3_list[j][2]]
                and A[ap3_list[j][1]][ap3_list[j][3]] and not A[ap3_list[j][0]][ap3_list[j][3]])
    n_nondt = sum(1 for j in range(len(ap3_list))
                  if abs(filling[j]) > 1e-8 and not (A[ap3_list[j][0]][ap3_list[j][2]]
                  and A[ap3_list[j][1]][ap3_list[j][3]]))
    print(f"\n  Filling composition: DT_f={n_dtf}, DT_c={n_dtc}, non-DT={n_nondt}")

    # KEY: The simplicial hole is filled by DT_c paths!
    # Show which DT_c paths contribute
    if n_dtc > 0:
        print(f"\n  DT_c paths in the filling:")
        for j in range(len(ap3_list)):
            if abs(filling[j]) < 1e-8: continue
            a, b, c, d = ap3_list[j]
            if A[a][c] and A[b][d] and not A[a][d]:
                # DT_c: d→a creates 3-cycle + 4-cycle structure
                # (a,c,d) is NT: a→c→d→a
                # (a,b,d) is NT: a→b→d→a
                print(f"    {filling[j]:+.4f} * ({a},{b},{c},{d}) [DT_c: d→a]")
                print(f"      Cycles: a→c→d→a and a→b→d→a")
                print(f"      NT boundary faces: ({a},{c},{d}) and ({a},{b},{d})")

    # Also check: does this cycle have NT terms when expressed as a path Z₂ cycle?
    # The simplicial cycle is in TT only. But in Ω₂, it maps to TT + possibly NT.
    # Wait — the cycle is already expressed in TT only. Let me check if it's in Z₂.
    bd2 = build_full_boundary_matrix(ap2, ap1)
    boundary = bd2 @ cycle_a2
    bd_norm = np.linalg.norm(boundary)
    print(f"\n  Is the simplicial cycle in Z₂? ||∂₂(cycle)|| = {bd_norm:.6f}")

    if bd_norm < 1e-8:
        print("  ✓ Yes! The simplicial cycle IS in path Z₂.")
        print("  So the simplicial hole in flag complex IS also a hole in path complex.")
        print("  But path homology fills it via DT_c paths!")
    else:
        print("  ✗ No! The simplicial cycle is NOT in path Z₂.")
        print("  The simplicial cycle isn't even a path cycle — it needs NT corrections.")
        # Express the boundary
        for j in range(len(ap1)):
            if abs(boundary[j]) > 1e-8:
                arc = tuple(ap1[j])
                print(f"    boundary[{arc}] = {boundary[j]:+.4f}")

print(f"\n\nTotal tournaments with flag H₂ > 0 at n=5, c₃=4: {count}")
print("\nDone.")
