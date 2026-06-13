#!/usr/bin/env python3
r"""beta2_proof_final.py - Complete proof strategy for beta_2 = 0

THEOREM: For every tournament T on n >= 3 vertices, beta_2(T) = 0.

PROOF STRATEGY:
Base case: n <= 4 (verified exhaustively, beta_2 = 0 trivially since dim(Z_2)=rk(d_3)=0).
Actually: at n=3, Om_3 is empty, Z_2 is trivially 0. Same at n=4.
At n=5: first non-trivial case.

Induction step: Assume beta_2(T') = 0 for all tournaments T' on < n vertices.
For an n-vertex tournament T, we use the LES:
  0 = H_2(T\v) -> H_2(T) -> H_2(T,T\v) -> H_1(T\v) -> H_1(T)

By induction, H_2(T\v) = 0.
So H_2(T) injects into H_2(T,T\v).
H_2(T,T\v) = 0 iff i_*: H_1(T\v) -> H_1(T) is injective.
i_* injective iff b_1(T\v) <= b_1(T).

We need: EXISTS v with b_1(T\v) <= b_1(T) ("good vertex").

PROOF OF GOOD VERTEX EXISTENCE:
Claim 1 (HYP-279): b_1(T) in {0, 1} for all tournaments T.
Claim 2 (HYP-282): Sum_v b_1(T\v) <= 3 for all tournaments T.

Together: If b_1(T) = 0, then Sum_v b_1(T\v) <= 3, so at most 3 vertices
have b_1(T\v) = 1. For n >= 4, at least n-3 >= 1 vertex has b_1(T\v) = 0 = b_1(T). GOOD.
If b_1(T) = 1, then ALL vertices are good (b_1(T\v) in {0,1} <= 1 = b_1(T)). GOOD.

NOW: Can we PROVE Claims 1 and 2?

CLAIM 1 APPROACH: b_1(T) = dim(Z_1) - rk(d_2).
dim(Z_1) = (n-1)(n-2)/2 (constant for all tournaments on n vertices).
We need: rk(d_2) >= (n-1)(n-2)/2 - 1.

d_2: Om_2 -> Om_1. The image of d_2 lives in Z_1 (since d_1 d_2 = 0).
We need: the image of d_2 has codimension at most 1 in Z_1.

Key structural fact: Z_1 is the "cycle space" of the tournament viewed as a graph.
It has dimension C(n,2) - (n-1) = (n-1)(n-2)/2.
A basis consists of fundamental cycles of any spanning tree.

For a tournament T, consider any 3-cycle i->j->k->i.
The boundary d_2(i,j,k) = (j,k) - (i,k) + (i,j).
This is a 1-cycle (its boundary is k - i + j - j + i - k = 0).
So im(d_2) contains all such "triangle cycles".

QUESTION: Do triangle cycles span Z_1 (up to codimension 1)?

CLAIM 2 APPROACH: Double-counting.
Sum_v b_1(T\v) = Sum_v [dim Z_1(T\v) - rk d_2(T\v)]
= Sum_v [(n-2)(n-3)/2 - rk d_2(T\v)]
= n*(n-2)(n-3)/2 - Sum_v rk d_2(T\v)

We need: Sum_v rk d_2(T\v) >= n*(n-2)(n-3)/2 - 3.

Let me verify this and look for patterns.

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


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


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


def compute_b1(A, n):
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths0 = [(i,) for i in range(n)]
    paths2 = enumerate_allowed_paths(A, n, 2)
    if not paths1:
        return 0
    omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
    dim_O1 = omega1.shape[1] if omega1.ndim == 2 else 0
    if dim_O1 == 0:
        return 0
    D1 = build_full_boundary_matrix([tuple(p) for p in paths1], paths0)
    D1_om = D1 @ omega1
    sv = np.linalg.svd(D1_om, compute_uv=False)
    rk_d1 = int(sum(s > 1e-8 for s in sv))
    if paths2:
        omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
        dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
        if dim_O2 > 0:
            D2 = build_full_boundary_matrix([tuple(p) for p in paths2],
                                            [tuple(p) for p in paths1])
            D2_om = D2 @ omega2
            sv2 = np.linalg.svd(D2_om, compute_uv=False)
            rk_d2 = int(sum(s > 1e-8 for s in sv2))
        else:
            rk_d2 = 0
    else:
        rk_d2 = 0
    return dim_O1 - rk_d1 - rk_d2


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


# ============================================================
# Part 1: Triangle cycles vs Z_1 — do they span?
# ============================================================
print("=" * 70)
print("TRIANGLE CYCLES vs Z_1 — CODIMENSION")
print("=" * 70)

# For each tournament, compute:
# 1. dim(Z_1)
# 2. rank of the matrix of all triangle-cycle vectors in Z_1
# 3. Codimension = dim(Z_1) - rank

for n in [4, 5, 6]:
    total = 2 ** (n * (n - 1) // 2)
    codim_dist = Counter()
    t0 = time.time()

    for bits in range(total):
        A = [[0] * n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        # Enumerate all 3-cycles
        triangles = []
        for i in range(n):
            for j in range(n):
                if j == i or not A[i][j]:
                    continue
                for k in range(n):
                    if k == i or k == j or not A[j][k] or not A[k][i]:
                        continue
                    triangles.append((i, j, k))

        # Build triangle cycle vectors in edge space
        # Each triangle (i,j,k) with i->j->k->i gives cycle:
        # (i,j) + (j,k) + (k,i) = (i,j) + (j,k) - (i,k)
        # But we need signs: d_2(i,j,k) = (j,k) - (i,k) + (i,j)
        # which is in the edge space R^{C(n,2)} indexed by ordered pairs.

        num_edges = n * (n - 1) // 2
        edge_to_idx = {}
        idx_e = 0
        for i in range(n):
            for j in range(i + 1, n):
                edge_to_idx[(i, j)] = idx_e
                idx_e += 1

        if not triangles:
            codim_dist[('no_triangles', n)] += 1
            continue

        # Triangle cycle matrix
        # Each triangle (a,b,c) with a->b->c->a gives cycle vector in Z_1.
        # In the ORIENTED edge space (using edges as a->b for a<b or b->a for a>b):
        # We work in R^{edges} where edge (a,b) with a<b carries +1 if a->b, -1 if b->a.
        # The triangle cycle is: e_{a,b} + e_{b,c} + e_{c,a}
        # where e_{x,y} means +1 at edge {min(x,y),max(x,y)} if x->y, -1 if y->x.

        # Actually, for path homology the edge space is R^{A_1} = R^{directed edges}.
        # Each directed edge i->j is a basis vector.
        # d_2(i,j,k) = (j,k) - (i,k) + (i,j) in R^{directed edges}.
        # But (i,k) means the directed edge i->k. If k->i in the tournament, (i,k) is NOT an edge.
        # In that case, d_2(i,j,k) = (j,k) + (i,j) - (i,k) where (i,k) may or may not be an edge.

        # Wait, in the path homology formulation:
        # A_1 = set of allowed 1-paths = edges of the tournament (directed).
        # d_2(a,b,c) = d_0(a,b,c) - d_1(a,b,c) + d_2(a,b,c)
        #            = (b,c) - (a,c) + (a,b)
        # where (a,c) is a directed edge a->c. For a TOURNAMENT, either a->c or c->a.
        # If a->c: (a,c) is in A_1. So the boundary is (b,c) - (a,c) + (a,b).
        # If c->a: (a,c) is NOT in A_1. But then the 2-path (a,b,c) is an "NT" path.
        #
        # For a 3-cycle a->b->c->a, we have c->a (NOT a->c). So (a,b,c) is NOT a TT path.
        # d_2(a,b,c) = (b,c) - (a,c) + (a,b). But (a,c) is not in A_1 since c->a.
        # Hmm, but d_2 maps to R^{A_1}. If (a,c) is not in A_1, then d_2(a,b,c) puts -1 on
        # a path that doesn't exist?

        # I think I need to go back to the definition more carefully.
        # In Grigoryan's path homology:
        # d_i : R^{A_p} -> R^{A_{p-1}} is defined as:
        # d_i(v_0,...,v_p) = (v_0,...,v_{i-1},v_{i+1},...,v_p) IF this is an allowed path
        #                  = 0 otherwise.
        # The boundary map is partial = sum (-1)^i d_i.
        #
        # For (a,b,c) with a->b->c:
        # d_0(a,b,c) = (b,c) [always allowed since b->c]
        # d_1(a,b,c) = (a,c) [allowed iff a->c, i.e., TT]
        # d_2(a,b,c) = (a,b) [always allowed since a->b]
        # partial(a,b,c) = (b,c) - (a,c) + (a,b) if TT
        #               = (b,c) - 0 + (a,b) = (b,c) + (a,b) if NT
        #
        # For a 3-cycle a->b->c->a (NT since c->a not a->c):
        # partial(a,b,c) = (b,c) + (a,b)
        # This is NOT a 1-cycle! partial_1((b,c) + (a,b)) = (c)-(b) + (b)-(a) = (c)-(a) != 0.
        #
        # So triangle cycles don't directly come from d_2 of a single NT 2-path.
        # They come from Om_2 elements which are COMBINATIONS of 2-paths.

        # For a 3-cycle a->b->c->a:
        # The three 2-paths are (a,b,c), (b,c,a), (c,a,b).
        # Each is NT. Their boundaries:
        # partial(a,b,c) = (b,c) + (a,b) [since c->a]
        # partial(b,c,a) = (c,a) + (b,c) [since a->b, not b->a... wait]
        # Actually (b,c,a): d_0=(c,a), d_1=(b,a), d_2=(b,c).
        # Is (b,a) allowed? b->a iff A[b][a]=1.
        # In 3-cycle a->b->c->a: we DON'T have b->a. So d_1(b,c,a) = 0.
        # partial(b,c,a) = (c,a) + (b,c)
        #
        # Similarly: partial(c,a,b) = (a,b) + (c,a)
        #
        # Sum of all three: 2(a,b) + 2(b,c) + 2(c,a)
        # Not zero! So the sum of the three NT paths is NOT a cycle.
        # But 1/2 of the sum is (a,b)+(b,c)+(c,a) which has boundary
        # (b)-(a)+(c)-(b)+(a)-(c) = 0. It IS a 1-cycle.
        #
        # But this 1-cycle is NOT the image of an Omega_2 element...
        # unless the alternating sum works.
        #
        # Actually, (a,b,c)-(b,c,a)+(c,a,b):
        # partial = [(b,c)+(a,b)] - [(c,a)+(b,c)] + [(a,b)+(c,a)]
        #         = 2(a,b) - 0 = 2(a,b)
        # Not a cycle either.
        #
        # The issue is that for NT paths, d_1 vanishes, so the alternating
        # sum doesn't telescope nicely.

        # Let me just compute im(d_2) directly and find its codimension in Z_1.
        # Use the library.
        paths0 = [(i,) for i in range(n)]
        paths1 = enumerate_allowed_paths(A, n, 1)
        paths2 = enumerate_allowed_paths(A, n, 2)

        omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
        dim_O1 = omega1.shape[1] if omega1.ndim == 2 else 0
        D1 = build_full_boundary_matrix([tuple(p) for p in paths1], paths0)
        D1_om = D1 @ omega1
        sv1 = np.linalg.svd(D1_om, compute_uv=False)
        rk_d1 = int(sum(s > 1e-8 for s in sv1))
        dim_Z1 = dim_O1 - rk_d1

        if paths2:
            omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
            dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
            if dim_O2 > 0:
                D2 = build_full_boundary_matrix([tuple(p) for p in paths2],
                                                [tuple(p) for p in paths1])
                D2_om = D2 @ omega2
                sv2 = np.linalg.svd(D2_om, compute_uv=False)
                rk_d2 = int(sum(s > 1e-8 for s in sv2))
            else:
                rk_d2 = 0
        else:
            rk_d2 = 0

        codim = dim_Z1 - rk_d2  # = b1
        codim_dist[codim] += 1

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): codim(im d_2 in Z_1)")
    for c, cnt in sorted(codim_dist.items(), key=lambda x: str(x)):
        print(f"  codim={c} (b1={c}): {cnt}")


# ============================================================
# Part 2: TT cycles — what do they span?
# ============================================================
print(f"\n{'=' * 70}")
print("TT CYCLE SPACE vs Z_1")
print("=" * 70)

# For a TT path (a,b,c): a->b->c and a->c.
# partial(a,b,c) = (b,c) - (a,c) + (a,b).
# This IS a 1-cycle: partial_1 = (c)-(b) - (c)+(a) + (b)-(a) = 0.
# So every TT boundary is a 1-cycle.
#
# Question: do TT boundaries span im(d_2)?
# Or equivalently: do TT elements span Om_2?

# Each TT triple (a,b,c) is an element of Om_2 (it satisfies d_1=0 since a->c exists).
# NT elements of Om_2 are combinations of NT paths that cancel under d_1.

# The boundary of a TT element is a "triangle cycle" in the usual sense.
# If TT elements span all of Om_2, then TT boundaries span im(d_2).

# Let's check: does im(d_2) = im(d_2 restricted to TT elements)?
for n in [5, 6]:
    total = 2 ** (n * (n - 1) // 2)
    tt_suffices = 0
    tt_fails = 0
    max_deficit = 0
    t0 = time.time()

    for bits in range(total):
        A = [[0] * n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        paths1 = enumerate_allowed_paths(A, n, 1)
        paths2 = enumerate_allowed_paths(A, n, 2)
        paths0 = [(i,) for i in range(n)]

        if not paths2:
            tt_suffices += 1
            continue

        # Find TT 2-paths
        tt_paths = []
        for p in paths2:
            a, b, c = p[0], p[1], p[2]
            if A[a][c]:  # a->c means TT
                tt_paths.append(p)

        # Full d_2 rank (restricted to Om_2)
        omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
        omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
        dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

        if dim_O2 == 0:
            tt_suffices += 1
            continue

        D2 = build_full_boundary_matrix([tuple(p) for p in paths2],
                                        [tuple(p) for p in paths1])
        D2_om = D2 @ omega2
        sv = np.linalg.svd(D2_om, compute_uv=False)
        rk_full = int(sum(s > 1e-8 for s in sv))

        # TT-only d_2 rank
        # Build the TT subspace of Om_2
        # Each TT path is a single basis vector in the path space.
        # Its projection onto Om_2 is... actually TT paths ARE in Om_2 individually.
        # So we just need the rank of d_2 restricted to TT paths.
        if tt_paths:
            path_to_idx = {tuple(p): i for i, p in enumerate(paths2)}
            tt_indices = [path_to_idx[tuple(p)] for p in tt_paths]
            # TT columns of D2 in path1 coords
            D2_tt = D2[:, tt_indices]
            sv_tt = np.linalg.svd(D2_tt, compute_uv=False)
            rk_tt = int(sum(s > 1e-8 for s in sv_tt))
        else:
            rk_tt = 0

        deficit = rk_full - rk_tt
        max_deficit = max(max_deficit, deficit)

        if deficit == 0:
            tt_suffices += 1
        else:
            tt_fails += 1

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s):")
    print(f"  TT boundaries span im(d_2): {tt_suffices}/{tt_suffices+tt_fails}")
    print(f"  Max deficit (TT vs full): {max_deficit}")


# ============================================================
# Part 3: b1=1 characterization — when does cokernel appear?
# ============================================================
print(f"\n{'=' * 70}")
print("b1=1 CHARACTERIZATION")
print("=" * 70)

# b1=1 iff rk(d_2) = dim(Z_1) - 1.
# dim(Z_1) = (n-1)(n-2)/2.
# rk(d_2) <= min(dim(Om_2), dim(Z_1)).
# So b1=1 requires dim(Om_2) >= dim(Z_1) - 1 = (n-1)(n-2)/2 - 1.

# At n=5: dim(Z_1) = 6. Need dim(Om_2) >= 5 and rk(d_2) = 5.
# All b1=1 at n=5 have dim(Om_2) in {8, 9, 10} >> 5.
# So the rank deficiency is NOT from too few Om_2 dimensions.
# It's from a specific dependency in the image.

# Key question: b1=1 iff T is strongly connected (at small n).
# At n=5: NOT SC => b1=0 (480 tournaments). SC => b1 in {0,1}.
# SC + b1=0: 240. SC + b1=1: 304.
# So most SC tournaments at n=5 have b1=1.

print("\nn=5 breakdown:")
n = 5
sc_b1 = Counter()
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
    sc = is_sc(A, n)
    b1 = compute_b1(A, n)
    sc_b1[(sc, b1)] += 1

for (s, b), cnt in sorted(sc_b1.items()):
    print(f"  SC={s}, b1={b}: {cnt}")

print("\nn=6 breakdown:")
n = 6
sc_b1 = Counter()
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
    sc = is_sc(A, n)
    b1 = compute_b1(A, n)
    sc_b1[(sc, b1)] += 1

for (s, b), cnt in sorted(sc_b1.items()):
    print(f"  SC={s}, b1={b}: {cnt}")


# ============================================================
# Part 4: Condensation gives b1=0 for non-SC
# ============================================================
print(f"\n{'=' * 70}")
print("b1=0 FOR NON-SC TOURNAMENTS — ALGEBRAIC PROOF")
print("=" * 70)

# A non-SC tournament T has a condensation into SC components C_1, ..., C_k
# with a transitive ordering. The path homology of T decomposes:
# H_1(T) = direct sum of H_1(C_i)? No, this is not generally true.
#
# But we know: NOT SC => b1=0 (verified exhaustively at n<=6).
# This should follow from the structure of Om_1 and Z_1.
#
# For a non-SC tournament: there exists a partition V = S cup S^c such that
# all edges from S^c to S are absent (no edge from S^c to S).
# This means: the boundary map d_1 has a "block" structure.
#
# Actually, the key fact is:
# If T is NOT SC, then T is a "Rees product" of its SC components.
# The 1-cycles in T are generated by 1-cycles within each component,
# and all such cycles are boundaries of 2-chains.
#
# Formal argument: If T is not SC, it has an acyclic condensation.
# The 2-paths in T that cross component boundaries contribute to d_2,
# and their images span all of Z_1.
#
# Let me verify: for non-SC tournaments, is dim(Om_2) always >= dim(Z_1)?

print("\nNon-SC: dim(Om_2) vs dim(Z_1)")
n = 5
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
    if is_sc(A, n):
        continue

    paths0 = [(i,) for i in range(n)]
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths2 = enumerate_allowed_paths(A, n, 2)

    omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
    dim_O1 = omega1.shape[1] if omega1.ndim == 2 else 0
    D1 = build_full_boundary_matrix([tuple(p) for p in paths1], paths0)
    D1_om = D1 @ omega1
    sv1 = np.linalg.svd(D1_om, compute_uv=False)
    rk_d1 = int(sum(s > 1e-8 for s in sv1))
    dim_Z1 = dim_O1 - rk_d1

    if paths2:
        omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
        dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
    else:
        dim_O2 = 0

    if dim_O2 < dim_Z1:
        print(f"  bits={bits}: dim_O2={dim_O2} < dim_Z1={dim_Z1} !!!")
        break
else:
    print(f"  All non-SC at n=5: dim(Om_2) >= dim(Z_1)")


# ============================================================
# Part 5: For SC tournaments with b1=1: WHAT is the cokernel?
# ============================================================
print(f"\n{'=' * 70}")
print("b1=1 SC TOURNAMENTS: COKERNEL GENERATOR")
print("=" * 70)

# The cokernel of d_2 in Z_1 has dimension 1 when b1=1.
# This generator is a 1-cycle that is NOT a boundary.
# For the 3-cycle at n=3: the cokernel is the cycle itself (only cycle, Om_2 empty).
# For n=4: score (1,1,2,2) (4-cycle with a diagonal).
# For n=5: various SC tournaments.
#
# KEY QUESTION: Is the cokernel generator always related to a "long cycle"
# (Hamiltonian cycle or near-Hamiltonian cycle) in the tournament?

# At n=3: 3-cycle = Hamiltonian cycle. Cokernel = Hamiltonian cycle. YES.
# At n=4: The tournament has a 4-cycle (Hamiltonian). Is cokernel related?

# Let me compute the cokernel generator for a few b1=1 tournaments at n=5.
n = 5
count = 0
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
    if not is_sc(A, n):
        continue
    b1 = compute_b1(A, n)
    if b1 != 1:
        continue

    count += 1
    if count > 3:
        break

    paths0 = [(i,) for i in range(n)]
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths2 = enumerate_allowed_paths(A, n, 2)

    omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
    dim_O1 = omega1.shape[1]

    D1 = build_full_boundary_matrix([tuple(p) for p in paths1], paths0)
    D1_om = D1 @ omega1

    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

    D2 = build_full_boundary_matrix([tuple(p) for p in paths2],
                                    [tuple(p) for p in paths1])
    D2_om = D2 @ omega2

    # Find Z_1 = ker(D1_om) in Om_1 coordinates
    # Use SVD
    U, S_vals, Vt = np.linalg.svd(D1_om, full_matrices=True)
    rk = int(sum(s > 1e-8 for s in S_vals))
    Z1_basis = Vt[rk:]  # null space vectors, shape (dim_Z1, dim_O1)

    # Project d_2 image into Z_1
    D2_in_om1 = np.linalg.lstsq(omega1, D2_om, rcond=None)[0]
    D2_in_Z1 = Z1_basis @ D2_in_om1  # shape (dim_Z1, dim_O2)

    # Find cokernel of D2_in_Z1
    U2, S2, Vt2 = np.linalg.svd(D2_in_Z1, full_matrices=True)
    rk2 = int(sum(s > 1e-8 for s in S2))
    coker_in_Z1 = U2[:, rk2:]  # shape (dim_Z1, corank)

    if coker_in_Z1.shape[1] > 0:
        # Convert to path coordinates
        coker_in_Om1 = coker_in_Z1[:, 0] @ Z1_basis  # shape (dim_O1,)
        path_coords = omega1 @ coker_in_Om1  # shape (|paths1|,)

        nonzero_edges = []
        for j in range(len(paths1)):
            if abs(path_coords[j]) > 1e-8:
                nonzero_edges.append((tuple(paths1[j]), round(float(path_coords[j]), 4)))

        scores = tuple(sorted(sum(A[i]) for i in range(n)))
        print(f"\n  bits={bits}, scores={scores}")
        print(f"    Cokernel generator ({len(nonzero_edges)} edges):")
        for e, c in sorted(nonzero_edges, key=lambda x: -abs(x[1])):
            print(f"      {e}: {c}")


print("\n\nDone.")
