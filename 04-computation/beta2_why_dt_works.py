#!/usr/bin/env python3
"""
WHY DT PATHS FILL ALL 2-CYCLES: ALGEBRAIC ANALYSIS

DT 4-path (a,b,c,d): a→b→c→d AND a→c AND b→d
Boundary: ∂(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)

KEY OBSERVATION: A DT 4-path is exactly a directed 4-clique (4 vertices
with all 6 edges going in a consistent direction forming a total order).

No wait — not quite. A DT 4-path (a,b,c,d) requires:
  a→b, b→c, c→d (path edges)
  a→c, b→d (DT condition)
But does NOT require a→d. In a tournament, either a→d or d→a.

Actually: a→b→c→d with a→c and b→d. What about a→d?
a→b→d (transitive) but we don't know a→d.
a→c→d (transitive) but we don't know a→d.
So a→d is NOT guaranteed.

HOWEVER: the 4 FACES are:
  (b,c,d): b→c→d and b→d ✓ (transitive)
  (a,c,d): a→c→d and a→d? We need a→d for this to be transitive!
  (a,b,d): a→b→d and a→d? Same issue!
  (a,b,c): a→b→c and a→c ✓

Wait, Face 1 = (a,c,d) requires a→c (yes) and c→d (yes) and a→d.
If d→a, then (a,c,d) has edges a→c, c→d, d→a — that's NOT in A_2!

HOLD ON. Let me re-examine what "DT 4-path" actually means for faces.

For path (v0,v1,v2,v3):
  Faces: F0=(v1,v2,v3), F1=(v0,v2,v3), F2=(v0,v1,v3), F3=(v0,v1,v2)

  F0=(v1,v2,v3): edges v1→v2, v2→v3 always exist. Is a 2-path in A_2.
  F3=(v0,v1,v2): edges v0→v1, v1→v2 always exist. Is a 2-path in A_2.
  F1=(v0,v2,v3): IS in A_2 iff v0→v2 edge exists AND v2→v3 exists (yes).
    So F1 ∈ A_2 iff v0→v2.
  F2=(v0,v1,v3): IS in A_2 iff v0→v1 (yes) AND v1→v3 edge exists.
    So F2 ∈ A_2 iff v1→v3.

"DT" = v0→v2 AND v1→v3 (both conditions).
So F1 ∈ A_2 and F2 ∈ A_2.

But wait: F1=(v0,v2,v3) is in A_2 (is an allowed 2-path) iff v0→v2 AND v2→v3.
Being in A_2 does NOT mean it's in Ω_2!
Ω_2 = transitive triples = those (a,b,c) with a→b→c AND a→c.

F1=(v0,v2,v3): v0→v2, v2→v3, and v0→v3? Not guaranteed!
F2=(v0,v1,v3): v0→v1, v1→v3, and v0→v3? Not guaranteed!

So DT condition ensures F1,F2 ∈ A_2, but NOT necessarily ∈ Ω_2!

However, the boundary ∂_3(v0,v1,v2,v3) is in A_2 (all faces are allowed 2-paths).
But is it in Ω_2? For ∂_3 to map to Ω_2, we need the faces themselves to be
∂-invariant (in Ω_2). That's the Ω_3 condition!

Wait — ∂_3 maps Ω_3 → Ω_2. So for a DT path, ∂_3(DT) is in A_2.
But is ∂_3(DT) in Ω_2?

For any x ∈ Ω_3, ∂x ∈ A_2, and then automatically ∂x ∈ Ω_2 because
∂²x = 0 which means ∂(∂x) ∈ A_1, so ∂x ∈ Ω_2.

So yes! DT ⊂ Ω_3 and ∂(DT) ⊂ Ω_2.

But the COEFFICIENTS of ∂(DT) in the Ω_2 basis include ALL transitive triple
faces, even those that aren't individually DT faces. The point is that
∂(DT) is a chain in A_2 that lies in Ω_2.

OK so the question is really:
∂_3: span(DT) → Ω_2 surjects onto ker(∂_2|Ω_2).

Let me compute this more carefully.
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

# ===== Understand what DT means for 4-vertex subsets =====
print("=" * 70)
print("DT 4-PATHS vs TRANSITIVE 4-CLIQUES")
print("=" * 70)

n = 5
dt_is_clique = Counter()  # does DT => a→d?

for A in all_tournaments_gen(n):
    a3 = enumerate_allowed_paths(A, n, 3)
    for p in a3:
        v0,v1,v2,v3 = p
        if A[v0][v2]==1 and A[v1][v3]==1:  # DT
            has_ad = A[v0][v3] == 1
            dt_is_clique[has_ad] += 1

print(f"DT 4-paths with a→d (full 4-clique): {dt_is_clique[True]}")
print(f"DT 4-paths with d→a (NOT 4-clique): {dt_is_clique[False]}")
print(f"Fraction that are 4-cliques: {dt_is_clique[True]/(dt_is_clique[True]+dt_is_clique[False]):.4f}")

# ===== What ARE the 2-cycles? =====
print(f"\n\n{'='*70}")
print("STRUCTURE OF 2-CYCLES IN ker(∂_2|Ω_2)")
print("="*70)

# A 2-cycle z = Σ c_i (a_i,b_i,c_i) with ∂z = 0.
# ∂(a,b,c) = (b,c) - (a,c) + (a,b)
# So ∂z = 0 means: for each edge (u,v), the sum of coefficients where
# (u,v) appears as a face equals zero.
#
# Edge (u,v) appears:
#   as face (b,c) with +1 when (a,u,v) is a transitive triple
#   as face (a,c) with -1 when (u,a,v) is... wait, face=(a,c) means (a,c)=(u,v),
#   so a=u, c=v. That means we need triple (u,b,v) with u→b→v and u→v (transitive).
#
# More precisely: (u,v) appears as:
#   Face 0 = (b,c) when b=u,c=v: from any triple (a,u,v) with coefficient +1
#   Face 1 = (a,c) when a=u,c=v: from any triple (u,b,v) with coefficient -1
#   Face 2 = (a,b) when a=u,b=v: from any triple (u,v,c) with coefficient +1
#
# So ∂z = 0 requires: for each edge u→v:
#   Σ_{a: (a,u,v)∈Ω_2} c_{(a,u,v)} - Σ_{b: (u,b,v)∈Ω_2} c_{(u,b,v)}
#   + Σ_{c: (u,v,c)∈Ω_2} c_{(u,v,c)} = 0

n = 5
cycle_patterns = Counter()

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

    # Get kernel basis
    U, S, Vt = np.linalg.svd(bd2, full_matrices=True)
    ker_basis = Vt[rank2:].T

    for k in range(ker_dim):
        vec = ker_basis[:, k]
        support = sum(1 for v in vec if abs(v) > 1e-8)
        # What vertices are involved?
        verts = set()
        for j in range(len(tt)):
            if abs(vec[j]) > 1e-8:
                verts.update(tt[j])
        cycle_patterns[(support, len(verts))] += 1

    if t_idx < 20 and ker_dim > 0:
        print(f"\nT#{t_idx}: ker_dim={ker_dim}")
        for k in range(min(ker_dim, 2)):
            vec = ker_basis[:, k]
            terms = [(tt[j], vec[j]) for j in range(len(tt)) if abs(vec[j]) > 1e-8]
            print(f"  cycle {k} ({len(terms)} terms):")
            for triple, coeff in terms:
                print(f"    {coeff:+.4f} · {triple}")

print(f"\n\nCycle patterns (support_size, n_vertices): count")
for key in sorted(cycle_patterns.keys()):
    print(f"  {key}: {cycle_patterns[key]}")

# ===== THE 4-VERTEX CYCLE =====
print(f"\n\n{'='*70}")
print("4-VERTEX CYCLES: THE BASIC BUILDING BLOCK")
print("="*70)

# At n=4, there are 24 tournaments with ker_dim=1.
# These all come from 4-cliques (transitive tournaments on 4 vertices).
# What does the single cycle look like?

n = 4
for t_idx, A in enumerate(all_tournaments_gen(n)):
    a2 = enumerate_allowed_paths(A, n, 2)
    tt = [tuple(p) for p in a2 if A[p[0]][p[2]] == 1]
    if len(tt) < 4:
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

    U, S, Vt = np.linalg.svd(bd2, full_matrices=True)
    ker_basis = Vt[rank2:].T

    print(f"\nT#{t_idx}: |Ω_2|={len(tt)}, ker_dim={ker_dim}")
    print(f"  Edges: ", end="")
    for i in range(n):
        for j in range(n):
            if A[i][j]: print(f"{i}→{j}", end=" ")
    print()
    print(f"  Transitive triples:")
    for t in tt:
        print(f"    {t}")

    vec = ker_basis[:, 0]
    terms = [(tt[j], vec[j]) for j in range(len(tt)) if abs(vec[j]) > 1e-8]
    print(f"  2-cycle: {len(terms)} terms")
    for triple, coeff in terms:
        print(f"    {coeff:+.4f} · {triple}")

    # DT paths
    a3 = enumerate_allowed_paths(A, n, 3)
    dt = [tuple(p) for p in a3 if A[p[0]][p[2]]==1 and A[p[1]][p[3]]==1]
    print(f"  DT 4-paths: {len(dt)}")
    for p in dt:
        v0,v1,v2,v3 = p
        print(f"    {p}:  ∂ = {(v1,v2,v3)} - {(v0,v2,v3)} + {(v0,v1,v3)} - {(v0,v1,v2)}")

    # Show that ∂(DT) = cycle
    if len(dt) > 0:
        tt_idx = {t: i for i, t in enumerate(tt)}
        bd3_dt = np.zeros((len(tt), len(dt)))
        for j, path in enumerate(dt):
            v0,v1,v2,v3 = path
            faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
            signs = [1, -1, 1, -1]
            for face, sign in zip(faces, signs):
                if face in tt_idx:
                    bd3_dt[tt_idx[face], j] += sign

        print(f"  ∂_3(DT) image vectors:")
        for j in range(len(dt)):
            col = bd3_dt[:, j]
            terms = [(tt[i], col[i]) for i in range(len(tt)) if abs(col[i]) > 0]
            print(f"    ∂{dt[j]} = ", end="")
            for triple, coeff in terms:
                print(f" {coeff:+.0f}·{triple}", end="")
            print()

    if t_idx > 5:
        break

print("\nDone.")
