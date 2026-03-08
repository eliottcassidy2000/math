#!/usr/bin/env python3
"""
DT 4-PATHS vs TRANSITIVE 4-CLIQUES

DT 4-path (a,b,c,d): a→b→c→d, a→c, b→d (6 conditions on 4 vertices + 2 non-adj)
Transitive 4-clique: total order a→b→c→d with ALL 6 edges (a→b,a→c,a→d,b→c,b→d,c→d)

DT does NOT require a→d. In tournament, either a→d or d→a.
- If a→d: it's a transitive 4-clique
- If d→a: NOT a 4-clique, but still DT

So DT paths ⊃ 4-clique paths. The "extra" DT paths are those with d→a.

The simplicial complex uses only 4-clique boundaries → H_2^{simp} can be nonzero.
GLMY uses DT boundaries (and more) → H_2^{GLMY} = 0.

The gap: DT paths with d→a provide the EXTRA boundaries that kill simplicial H_2.

Let's verify this explicitly.
"""
import numpy as np
from itertools import combinations
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import enumerate_allowed_paths

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

n = 5
print("=" * 70)
print("DT 4-PATHS vs TRANSITIVE 4-CLIQUES at n=5")
print("=" * 70)

dt_clique_count = 0
dt_nonclique_count = 0
h2_simp_nonzero = []

for t_idx, A in enumerate(all_tournaments_gen(n)):
    a3 = enumerate_allowed_paths(A, n, 3)
    a2 = enumerate_allowed_paths(A, n, 2)
    tt = [tuple(p) for p in a2 if A[p[0]][p[2]] == 1]
    dt_paths = [tuple(p) for p in a3 if A[p[0]][p[2]]==1 and A[p[1]][p[3]]==1]

    # Classify DT paths
    for p in dt_paths:
        v0,v1,v2,v3 = p
        if A[v0][v3]:
            dt_clique_count += 1
        else:
            dt_nonclique_count += 1

    # Check simplicial H_2 (4-clique boundaries only)
    if len(tt) == 0:
        continue

    a1 = enumerate_allowed_paths(A, n, 1)
    a1_idx = {tuple(p): i for i, p in enumerate(a1)}
    bd2 = np.zeros((len(a1), len(tt)))
    for j, (a, b, c) in enumerate(tt):
        bd2[a1_idx[(b,c)], j] += 1
        bd2[a1_idx[(a,c)], j] -= 1
        bd2[a1_idx[(a,b)], j] += 1

    rank2 = np.linalg.matrix_rank(bd2, tol=1e-8)
    ker_dim = len(tt) - rank2
    if ker_dim == 0:
        continue

    # 4-clique paths only (DT with a→d)
    clique_paths = [p for p in dt_paths if A[p[0]][p[3]]]
    nonclique_dt = [p for p in dt_paths if not A[p[0]][p[3]]]

    tt_idx = {t: i for i, t in enumerate(tt)}

    # im(∂_3 | clique only)
    if len(clique_paths) > 0:
        bd3_clique = np.zeros((len(tt), len(clique_paths)))
        for j, path in enumerate(clique_paths):
            v0,v1,v2,v3 = path
            faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
            signs = [1, -1, 1, -1]
            for face, sign in zip(faces, signs):
                if face in tt_idx:
                    bd3_clique[tt_idx[face], j] += sign
        im_clique = np.linalg.matrix_rank(bd3_clique, tol=1e-8)
    else:
        im_clique = 0

    # im(∂_3 | all DT)
    if len(dt_paths) > 0:
        bd3_dt = np.zeros((len(tt), len(dt_paths)))
        for j, path in enumerate(dt_paths):
            v0,v1,v2,v3 = path
            faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
            signs = [1, -1, 1, -1]
            for face, sign in zip(faces, signs):
                if face in tt_idx:
                    bd3_dt[tt_idx[face], j] += sign
        im_dt = np.linalg.matrix_rank(bd3_dt, tol=1e-8)
    else:
        im_dt = 0

    h2_clique = ker_dim - im_clique
    h2_dt = ker_dim - im_dt

    if h2_clique > 0:
        h2_simp_nonzero.append((t_idx, ker_dim, im_clique, im_dt,
                                 len(clique_paths), len(nonclique_dt)))

print(f"\nDT 4-paths that are 4-cliques: {dt_clique_count}")
print(f"DT 4-paths that are NOT 4-cliques (d→a): {dt_nonclique_count}")
print(f"Ratio non-clique: {dt_nonclique_count/(dt_clique_count+dt_nonclique_count):.4f}")

print(f"\nTournaments with simplicial H_2 > 0: {len(h2_simp_nonzero)}")
print(f"(ker_dim, im_clique, im_dt, #clique_paths, #nonclique_DT)")
from collections import Counter
pattern_counter = Counter()
for t_idx, ker, im_c, im_d, nc, nnd in h2_simp_nonzero:
    pattern_counter[(ker, im_c, im_d, nc, nnd)] += 1

for key in sorted(pattern_counter.keys()):
    ker, im_c, im_d, nc, nnd = key
    h2_c = ker - im_c
    h2_d = ker - im_d
    print(f"  ker={ker}, im_clique={im_c} (H_2^simp={h2_c}), im_DT={im_d} (H_2^DT={h2_d}), "
          f"clique={nc}, non-clique-DT={nnd}: {pattern_counter[key]}")

# ===== Detailed example =====
print(f"\n\n{'='*70}")
print("DETAILED EXAMPLE: TOURNAMENT WITH H_2^{simp} = 1 but H_2^{DT} = 0")
print("="*70)

for t_idx, ker, im_c, im_d, nc, nnd in h2_simp_nonzero[:3]:
    A = list(all_tournaments_gen(n))[t_idx]
    a3 = enumerate_allowed_paths(A, n, 3)
    a2 = enumerate_allowed_paths(A, n, 2)
    tt = [tuple(p) for p in a2 if A[p[0]][p[2]] == 1]
    dt_paths = [tuple(p) for p in a3 if A[p[0]][p[2]]==1 and A[p[1]][p[3]]==1]

    print(f"\nT#{t_idx}:")
    print(f"  Edges: ", end="")
    for i in range(n):
        for j in range(n):
            if A[i][j]: print(f"{i}→{j}", end=" ")
    print()

    # 3-cycle count
    t3 = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[i][k] and A[k][j] and A[j][i]):
                    t3 += 1
    print(f"  t3 = {t3}")

    print(f"  Transitive triples ({len(tt)}):")
    for t in tt:
        print(f"    {t}")

    print(f"\n  DT 4-paths ({len(dt_paths)}):")
    for p in dt_paths:
        v0,v1,v2,v3 = p
        clique = "4-clique" if A[v0][v3] else "d→a"
        print(f"    {p}  [{clique}]")

    # The non-clique DT path provides the extra boundary
    nonclique = [p for p in dt_paths if not A[p[0]][p[3]]]
    if nonclique:
        p = nonclique[0]
        v0,v1,v2,v3 = p
        print(f"\n  Non-clique DT {p}: edges {v0}→{v1}→{v2}→{v3}, {v0}→{v2}, {v1}→{v3}, {v3}→{v0}")
        faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
        signs = [1, -1, 1, -1]
        print(f"  ∂({p}) =")
        for face, sign in zip(faces, signs):
            in_tt = face in set(tt)
            status = '∈ Ω_2' if in_tt else '∈ A_2 minus Ω_2'
            print(f"    {'+' if sign>0 else '-'}({face})  {status}")

        # IMPORTANT: which faces are in Ω_2 vs A_2\Ω_2?
        a2_set = set(tuple(p) for p in a2)
        for face, sign in zip(faces, signs):
            if face not in a2_set:
                print(f"    *** {face} NOT IN A_2! ***")
            elif face not in set(tt):
                a,b,c = face
                print(f"    *** {face} in A_2 but NOT Ω_2: {a}→{b}={A[a][b]}, {b}→{c}={A[b][c]}, {a}→{c}={A[a][c]} ***")

print("\nDone.")
