#!/usr/bin/env python3
"""
β_2 = 0: CAN EVERY 2-CYCLE BE DECOMPOSED INTO 4-VERTEX CYCLES?

At n=4: every cycle = ∂(DT 4-path) on a transitive 4-clique.
At n=5: cycles have support 4-10. Can larger cycles be written as
        sums of 4-vertex cycles?

If YES: β_2 = 0 follows from β_2 = 0 at n=4 (trivially true).
If NO: need a different proof strategy.

THEOREM ATTEMPT: Every 2-cycle in Ω_2(T) for a tournament T is a
Z-linear combination of boundaries of DT 4-paths.

PROOF STRATEGY: The DT 4-paths come from 4-vertex transitive
subtournaments. Each such gives a simplex boundary. The cycles on
larger vertex sets should decompose into these.
"""
import numpy as np
from collections import Counter
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

# ===== At n=5: understand cycle generators =====
print("=" * 70)
print("CYCLE DECOMPOSITION INTO 4-VERTEX PIECES")
print("=" * 70)

n = 5
detailed_count = 0
for t_idx, A in enumerate(all_tournaments_gen(n)):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    tt = [tuple(p) for p in a2 if A[p[0]][p[2]] == 1]
    if len(tt) == 0:
        continue

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

    # Get all DT 4-paths
    a3 = enumerate_allowed_paths(A, n, 3)
    dt_paths = [tuple(p) for p in a3 if A[p[0]][p[2]]==1 and A[p[1]][p[3]]==1]

    # Each DT path lives on 4 vertices.
    # Group DT paths by their vertex set:
    from itertools import combinations
    dt_by_vset = {}
    for p in dt_paths:
        vset = frozenset(p)
        if vset not in dt_by_vset:
            dt_by_vset[vset] = []
        dt_by_vset[vset].append(p)

    # How many 4-cliques (transitive 4-subtournaments)?
    clique4_count = 0
    for quad in combinations(range(n), 4):
        # Check if the subtournament on quad is transitive
        is_trans = True
        for i_idx in range(4):
            for j_idx in range(i_idx+2, 4):
                if not A[quad[i_idx]][quad[j_idx]] and not A[quad[j_idx]][quad[i_idx]]:
                    is_trans = False
                    break
            if not is_trans:
                break
        # Actually: in a tournament, all edges exist. "Transitive" = no 3-cycle.
        # A tournament on 4 vertices is transitive iff it has 0 3-cycles.
        # Equivalently: it has a total order a→b→c→d with a→c, a→d, b→d.
        t3 = 0
        q = list(quad)
        for i_idx in range(4):
            for j_idx in range(i_idx+1, 4):
                for k_idx in range(j_idx+1, 4):
                    a,b,c = q[i_idx], q[j_idx], q[k_idx]
                    if (A[a][b] and A[b][c] and A[c][a]) or \
                       (A[a][c] and A[c][b] and A[b][a]):
                        t3 += 1
        if t3 == 0:
            clique4_count += 1

    if detailed_count < 5:
        detailed_count += 1
        print(f"\nT#{t_idx}: |Ω_2|={len(tt)}, ker_dim={ker_dim}, |DT|={len(dt_paths)}")
        print(f"  4-vertex transitive subtournaments: {clique4_count}")
        print(f"  DT paths by vertex set:")
        for vset, paths in sorted(dt_by_vset.items()):
            print(f"    {sorted(vset)}: {len(paths)} paths")
            for p in paths:
                print(f"      {p}")

# ===== Key: count DT paths per 4-vertex subset =====
print(f"\n\n{'='*70}")
print("DT PATHS PER TRANSITIVE 4-SUBSET")
print("="*70)

n = 5
dt_per_clique = Counter()
for A in all_tournaments_gen(n):
    a3 = enumerate_allowed_paths(A, n, 3)
    dt_paths = [tuple(p) for p in a3 if A[p[0]][p[2]]==1 and A[p[1]][p[3]]==1]

    from itertools import combinations
    for quad in combinations(range(n), 4):
        q = list(quad)
        t3 = 0
        for i_idx in range(4):
            for j_idx in range(i_idx+1, 4):
                for k_idx in range(j_idx+1, 4):
                    a,b,c = q[i_idx], q[j_idx], q[k_idx]
                    if (A[a][b] and A[b][c] and A[c][a]) or \
                       (A[a][c] and A[c][b] and A[b][a]):
                        t3 += 1
        if t3 == 0:
            # Count DT paths on this quad
            quad_set = frozenset(quad)
            dt_on_quad = [p for p in dt_paths if frozenset(p) == quad_set]
            dt_per_clique[len(dt_on_quad)] += 1

print("DT paths per transitive 4-subset: count")
for k in sorted(dt_per_clique.keys()):
    print(f"  {k}: {dt_per_clique[k]}")

# ===== How many DT paths does a transitive 4-tournament have? =====
print(f"\n\nA transitive tournament on 4 vertices (0→1→2→3, 0→2, 0→3, 1→3):")
print(f"DT 4-paths = sequences (a,b,c,d) with a→b→c→d, a→c, b→d")
A4 = [[0,1,1,1],[0,0,1,1],[0,0,0,1],[0,0,0,0]]  # 0→1→2→3
a3_4 = enumerate_allowed_paths(A4, 4, 3)
dt4 = [tuple(p) for p in a3_4 if A4[p[0]][p[2]]==1 and A4[p[1]][p[3]]==1]
print(f"DT paths: {dt4}")
for p in dt4:
    v0,v1,v2,v3 = p
    print(f"  {p}: a→c={A4[v0][v2]}, b→d={A4[v1][v3]}")

# Now check: how many 3-paths total?
a3_4_all = enumerate_allowed_paths(A4, 4, 3)
print(f"\nAll 3-paths on transitive T_4: {len(a3_4_all)}")
for p in a3_4_all:
    v0,v1,v2,v3 = p
    ac = A4[v0][v2]
    bd = A4[v1][v3]
    is_dt = ac and bd
    print(f"  {tuple(p)}: a→c={ac}, b→d={bd} {'DT' if is_dt else ''}")

# ===== CRITICAL: For transitive T_4, is the UNIQUE cycle = ∂(unique DT path)? =====
print(f"\n\n{'='*70}")
print("TRANSITIVE T_4: CYCLE = ∂(DT PATH)?")
print("="*70)

# On transitive T_4 = 0→1→2→3 with 0→2, 0→3, 1→3:
# Ω_2 = all transitive triples:
a2_4 = enumerate_allowed_paths(A4, 4, 2)
tt4 = [tuple(p) for p in a2_4 if A4[p[0]][p[2]] == 1]
print(f"Ω_2: {tt4}")

# Kernel
a1_4 = enumerate_allowed_paths(A4, 4, 1)
a1_idx = {tuple(p): i for i, p in enumerate(a1_4)}
bd2 = np.zeros((len(a1_4), len(tt4)))
for j, (a, b, c) in enumerate(tt4):
    bd2[a1_idx[(b,c)], j] += 1
    bd2[a1_idx[(a,c)], j] -= 1
    bd2[a1_idx[(a,b)], j] += 1

print(f"\n∂_2 matrix ({len(a1_4)} x {len(tt4)}):")
for i in range(len(a1_4)):
    print(f"  {a1_4[i]}: [{', '.join(f'{bd2[i,j]:+.0f}' for j in range(len(tt4)))}]")

rank2 = np.linalg.matrix_rank(bd2, tol=1e-8)
print(f"\nrank(∂_2) = {rank2}, ker_dim = {len(tt4) - rank2}")

if len(tt4) - rank2 > 0:
    U, S, Vt = np.linalg.svd(bd2, full_matrices=True)
    ker = Vt[rank2:].T[:, 0]
    print(f"Kernel vector: {ker}")
    # Normalize
    ker = ker / ker[np.argmax(np.abs(ker))]
    print(f"Normalized: {ker}")

    # DT boundary
    if len(dt4) > 0:
        tt_idx = {t: i for i, t in enumerate(tt4)}
        bd3 = np.zeros(len(tt4))
        v0,v1,v2,v3 = dt4[0]
        faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
        signs = [1, -1, 1, -1]
        for face, sign in zip(faces, signs):
            if face in tt_idx:
                bd3[tt_idx[face]] += sign

        bd3_norm = bd3 / bd3[np.argmax(np.abs(bd3))]
        print(f"∂(DT) normalized: {bd3_norm}")
        print(f"ker ∝ ∂(DT)? {np.allclose(ker, bd3_norm) or np.allclose(ker, -bd3_norm)}")

print("\nDone.")
