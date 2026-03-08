#!/usr/bin/env python3
"""
β_2 = 0 FOR TOURNAMENTS: DO DT 4-PATHS SUFFICE?

KEY QUESTION: Is im(∂_3 restricted to DT 4-paths) = ker(∂_2|Ω_2)?

If YES: proof reduces to showing DT-path boundaries span all 2-cycles.
If NO: need the full Ω_3 cancellation structure.

DT = "doubly transitive": (a,b,c,d) with a→b→c→d AND a→c AND b→d.
These are the 3-paths whose ALL 4 faces lie in A_2.

∂_3(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
For DT: all 4 faces are transitive triples, so ∂_3(DT) ⊂ Ω_2 automatically.
"""
import numpy as np
from itertools import permutations
from collections import Counter
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import enumerate_allowed_paths, compute_omega_basis

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def analyze_tournament(A, n):
    """Returns (ker_dim, im_dt_rank, im_full_rank, dim_om3, beta2)."""
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    # Ω_2 = transitive triples
    tt = [tuple(p) for p in a2 if A[p[0]][p[2]] == 1]
    if len(tt) == 0:
        return (0, 0, 0, 0, 0)

    # ∂_2: Ω_2 → A_1
    a1_idx = {tuple(p): i for i, p in enumerate(a1)}
    bd2 = np.zeros((len(a1), len(tt)))
    for j, (a, b, c) in enumerate(tt):
        bd2[a1_idx[(b,c)], j] += 1
        bd2[a1_idx[(a,c)], j] -= 1
        bd2[a1_idx[(a,b)], j] += 1

    rank2 = np.linalg.matrix_rank(bd2, tol=1e-8)
    ker_dim = len(tt) - rank2

    if ker_dim == 0:
        return (0, 0, 0, 0, 0)

    # DT 4-paths
    dt_paths = [tuple(p) for p in a3 if A[p[0]][p[2]]==1 and A[p[1]][p[3]]==1]
    tt_idx = {t: i for i, t in enumerate(tt)}

    # ∂_3 restricted to DT → Ω_2
    if len(dt_paths) > 0:
        bd3_dt = np.zeros((len(tt), len(dt_paths)))
        for j, path in enumerate(dt_paths):
            v0, v1, v2, v3 = path
            faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
            signs = [1, -1, 1, -1]
            for face, sign in zip(faces, signs):
                if face in tt_idx:
                    bd3_dt[tt_idx[face], j] += sign
        im_dt_rank = np.linalg.matrix_rank(bd3_dt, tol=1e-8)
    else:
        im_dt_rank = 0

    # Full Ω_3 and ∂_3
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    if dim_om3 > 0:
        a2_idx = {tuple(p): i for i, p in enumerate(a2)}
        bd3_full = np.zeros((len(a2), len(a3)))
        a3_tuples = [tuple(p) for p in a3]
        for j, path in enumerate(a3_tuples):
            v0, v1, v2, v3 = path
            faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
            signs = [1, -1, 1, -1]
            for face, sign in zip(faces, signs):
                if face in a2_idx:
                    bd3_full[a2_idx[face], j] += sign

        bd3_omega = bd3_full @ om3
        bd3_tt = np.zeros((len(tt), dim_om3))
        for i, t in enumerate(tt):
            bd3_tt[i, :] = bd3_omega[a2_idx[t], :]
        im_full_rank = np.linalg.matrix_rank(bd3_tt, tol=1e-8)
    else:
        im_full_rank = 0

    beta2 = ker_dim - im_full_rank
    return (ker_dim, im_dt_rank, im_full_rank, dim_om3, beta2)

# ===== n=4 =====
print("=" * 70)
print("β_2 = 0 ANALYSIS: DO DT 4-PATHS SUFFICE?")
print("=" * 70)

for n in [4, 5]:
    print(f"\n--- n={n} ---")
    stats = Counter()
    dt_sufficient = 0
    dt_insufficient = 0
    total_nontrivial = 0

    for A in all_tournaments_gen(n):
        ker, im_dt, im_full, om3, b2 = analyze_tournament(A, n)
        if ker > 0:
            total_nontrivial += 1
            if im_dt >= ker:
                dt_sufficient += 1
            else:
                dt_insufficient += 1
            stats[(ker, im_dt, im_full, om3, b2)] += 1

    print(f"  Tournaments with nontrivial ker(∂_2|Ω_2): {total_nontrivial}")
    print(f"  DT-only sufficient: {dt_sufficient}")
    print(f"  DT-only insufficient: {dt_insufficient}")
    print(f"\n  (ker, im_dt, im_full, dim_Ω_3, β_2): count")
    for key in sorted(stats.keys()):
        ker, im_dt, im_full, om3, b2 = key
        print(f"    ({ker}, {im_dt}, {im_full}, {om3}, {b2}): {stats[key]}")

# ===== n=6 (random sample) =====
print(f"\n--- n=6 (random sample of 1000) ---")
import random
random.seed(42)
n = 6
dt_sufficient = 0
dt_insufficient = 0
total_nontrivial = 0

for _ in range(1000):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    A = [[0]*n for _ in range(n)]
    for (i,j) in edges:
        if random.random() < 0.5:
            A[i][j] = 1
        else:
            A[j][i] = 1

    ker, im_dt, im_full, om3, b2 = analyze_tournament(A, n)
    if ker > 0:
        total_nontrivial += 1
        if im_dt >= ker:
            dt_sufficient += 1
        else:
            dt_insufficient += 1

print(f"  Tournaments with nontrivial ker: {total_nontrivial}")
print(f"  DT-only sufficient: {dt_sufficient}")
print(f"  DT-only insufficient: {dt_insufficient}")

# ===== Detailed analysis of DT-insufficient cases =====
print(f"\n\n{'='*70}")
print("DETAILED ANALYSIS: WHEN DT PATHS ARE INSUFFICIENT")
print("="*70)

n = 5
example_count = 0
for t_idx, A in enumerate(all_tournaments_gen(n)):
    ker, im_dt, im_full, om3, b2 = analyze_tournament(A, n)
    if ker > 0 and im_dt < ker and example_count < 3:
        example_count += 1
        print(f"\nTournament #{t_idx}: ker={ker}, im_dt={im_dt}, im_full={im_full}")

        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)
        tt = [tuple(p) for p in a2 if A[p[0]][p[2]] == 1]

        # Show the tournament
        print("  Edges:", end=" ")
        for i in range(n):
            for j in range(n):
                if A[i][j]: print(f"{i}→{j}", end=" ")
        print()

        # DT paths
        dt = [tuple(p) for p in a3 if A[p[0]][p[2]]==1 and A[p[1]][p[3]]==1]
        print(f"  DT 4-paths: {len(dt)}")
        for p in dt:
            print(f"    {p}")

        # Non-DT 3-paths and their face patterns
        non_dt = [tuple(p) for p in a3 if not (A[p[0]][p[2]]==1 and A[p[1]][p[3]]==1)]
        print(f"  Non-DT 3-paths: {len(non_dt)}")
        for p in non_dt[:5]:
            v0,v1,v2,v3 = p
            f1 = "a→c" if A[v0][v2] else "c→a"
            f2 = "b→d" if A[v1][v3] else "d→b"
            print(f"    {p}: {f1}, {f2}")

        # Full Ω_3 basis
        om3_basis = compute_omega_basis(A, n, 3, a3, a2)
        if om3_basis.ndim == 2 and om3_basis.shape[1] > 0:
            a3_tuples = [tuple(p) for p in a3]
            print(f"  Ω_3 basis ({om3_basis.shape[1]} vectors):")
            for col in range(min(om3_basis.shape[1], 5)):
                vec = om3_basis[:, col]
                terms = [(a3_tuples[j], vec[j]) for j in range(len(a3_tuples)) if abs(vec[j]) > 1e-8]
                dt_terms = sum(1 for p,_ in terms if A[p[0]][p[2]]==1 and A[p[1]][p[3]]==1)
                nondt_terms = len(terms) - dt_terms
                print(f"    v{col}: {len(terms)} terms ({dt_terms} DT, {nondt_terms} non-DT)")

if example_count == 0:
    print("\n*** NO DT-INSUFFICIENT CASES AT n=5! ***")
    print("DT 4-paths ALWAYS suffice to fill ker(∂_2|Ω_2)!")

print("\nDone.")
