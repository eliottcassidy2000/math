#!/usr/bin/env python3
"""
beta2_deformation_retract.py - Test deformation retract approach for beta_2=0

Inspired by Tang-Yau (2026): For circulant tournaments with S={1,...,d},
the path complex deformation retracts to C_n^{1,2}, giving H_m=0 for m>=2.

Key idea: Replace "long steps" with "short steps" using a chain homotopy.
For circulant digraphs, a "step" from a to b has length b-a mod n.
The homotopy h splits a long step s>=3 into (1, s-1) by inserting
an intermediate vertex.

For general tournaments, we don't have cyclic structure, but we DO have:
- Every pair connected by an edge (completeness)
- For any 2-path (a,b,c), vertex b is "between" a and c

Can we define a "vertex insertion" homotopy that works for all tournaments?

Strategy:
1. For each allowed 2-path (a,b,c) with a->b->c:
   - If a->c (transitive triple, "long step" from a to c):
     b is an intermediate vertex. The path already has "short steps" a->b, b->c.
   - If c->a (3-cycle, "non-transitive"):
     The face (a,c) is non-allowed. This is already "short" in some sense.

2. Key observation: In a tournament, the "diameter" of allowed 2-paths
   is bounded. There's no notion of "long step" since every pair is connected.

3. Alternative: Think of the deformation retract as reducing Omega_3 to
   something simpler. The Tang-Yau retract pi = Id - (dh + hd) projects
   chains to a subcomplex with only "1-step" or "2-step" transitions.

4. For tournaments: maybe define h using a TOPOLOGICAL SORT. Order
   vertices v_1 < v_2 < ... < v_n. Then for (a,b,c) with a < c, insert
   the vertex b' = (a + c) / 2... no, that doesn't make sense.

5. Actually, the key insight from swap cycles (delta-injectivity analysis):
   The only 2-cycles that matter are "swap cycles" (a,b,v)-(v,a,b).
   These involve a vertex v with v->a->b->v (3-path cycle).
   If we can fill each swap cycle with a 3-chain, we're done.

Let me test a SPECIFIC filling strategy for swap cycles.

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time, random
import numpy as np
from collections import defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis, path_betti_numbers
)
sys.stdout = _saved

random.seed(42)


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


# ============================================================
# PART 1: Explicit swap cycle filling test
# ============================================================
print("=" * 70)
print("SWAP CYCLE FILLING TEST")
print("=" * 70)

# The swap cycle sigma_{a,b}^v = (a,b,v) - (v,a,b)
# where v->a, a->b, b->v (so a in N+(v), b in N-(v), a->b)
#
# boundary: d_2(sigma) = (b,v)-(a,v)+(a,b) - [(a,b)-(v,b)+(v,a)]
#                       = (b,v)+(v,b) - (a,v) - (v,a)
# (the (a,b) terms cancel!)
#
# For sigma to be in Omega_2: it automatically is, since both
# (a,b,v) and (v,a,b) are allowed paths.
# But for a LINEAR COMBINATION to be in Omega_2, we need
# regularity: for each non-allowed face pair, coefficients sum to 0.
#
# For non-allowed (v,b) [since b->v, not v->b]:
#   Coefficient in d_2(sigma_{a,b}^v) = +1 (from face of (a,b,v))
#   So sum_a M_{ab} * 1 = 0 needed -> column sums = 0
#
# For non-allowed (a,v) [since v->a, not a->v]:
#   Coefficient in d_2(sigma_{a,b}^v) = -1 (from face of (a,b,v))
#   and -1 (from -(v,a) in d_2(v,a,b))
# Hmm wait, (v,a) IS an arc since v->a. And (a,v) is NOT an arc.
# Let me re-examine.
#
# d_2(a,b,v) = (b,v) - (a,v) + (a,b)
#   (b,v): ALLOWED (b->v since b in Q)
#   (a,v): NOT ALLOWED (v->a, so direction is wrong)
#   (a,b): ALLOWED (a->b)
#
# d_2(v,a,b) = (a,b) - (v,b) + (v,a)
#   (a,b): ALLOWED
#   (v,b): NOT ALLOWED (b->v, so wrong direction)
#   (v,a): ALLOWED (v->a)
#
# So sigma = (a,b,v) - (v,a,b):
# d_2(sigma) = [(b,v) - (a,v) + (a,b)] - [(a,b) - (v,b) + (v,a)]
#            = (b,v) - (a,v) - (v,a) + (v,b)
# Wait: (a,v) is NOT allowed, (v,a) is allowed. And (v,b) is NOT allowed.
#
# Hmm, d_2 maps to A_1 (all sequences), but Omega constraint says
# the non-allowed parts must vanish. The BOUNDARY d_2(sigma) as an
# A_1 element has:
#   Allowed: (b,v) with coeff +1, (v,a) with coeff -1
#   Non-allowed: (a,v) with coeff -1, (v,b) with coeff +1
#
# For sum M_{ab} sigma_{ab} to be in Omega_2:
#   We need the non-allowed parts of each individual path to be
#   absorbed by the Omega constraints. Actually, individual paths
#   (a,b,v) and (v,a,b) are already allowed 2-paths. The Omega_2
#   condition is on the CHAIN, not individual paths.
#   A chain sum alpha_p * e_p is in Omega_2 iff for each non-allowed
#   1-path (x,y): sum_{p with (x,y) as face} sign * alpha_p = 0.
#
# For non-allowed (a,v): paths with (a,v) as face are exactly
#   (a,b,v) with face (a,v) = delete b, sign = -(-1)^1 = +1?
#   Actually: d_2(a,b,v) = e_{b,v} - e_{a,v} + e_{a,b}
#   So (a,v) appears with coefficient -1 in d_2(a,b,v).
# WAIT. The Omega condition is NOT about d_2. It's about the
# constraint that projects A_2 onto Omega_2.
#
# Omega_2 = {f in A_2 : for each non-allowed (x,y),
#   sum_{z: (x,z,y) allowed} f(x,z,y) = 0}
# where the sum is over z forming allowed 2-paths with MIDDLE face (x,y).
#
# Hmm, actually the Omega condition at level 2 involves non-allowed
# 1-paths as the deletion of the middle vertex. Let me re-derive.
#
# The Omega_p complex is defined as:
# Omega_p = A_p / (d_{p+1}^{-1}(A_{p-1}) intersect A_{p+1} -> ...)
# Actually, Omega_p = {f in A_p : projections onto non-allowed directions vanish}
#
# More precisely: Omega_p = A_p cap (direct complement of im(d_{p+1}^hat))
# where d_{p+1}^hat is the "non-regular" part of the boundary.

# Let me just test computationally whether swap cycles are always
# boundaries.

n = 5
total_tests = 0
all_boundaries = True

for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)

    # For each vertex v, define P = N+(v), Q = N-(v)
    for v in range(n):
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]

        # Bipartite arcs from P to Q
        arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]
        if len(arcs_PQ) < 2:
            continue  # No swap cycle possible

        # Build the swap cycle basis (zero row/col sums on bipartite graph)
        # Variables: M[a,b] for each arc (a,b) from P to Q
        # Constraints: sum_b M[a,b] = 0 for each a in P
        #              sum_a M[a,b] = 0 for each b in Q

        m = len(arcs_PQ)
        # Build constraint matrix
        rows = []
        # Row sums
        for a in P:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if a2 == a:
                    row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)
        # Column sums
        for b in Q:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if b2 == b:
                    row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)

        if not rows:
            continue

        C = np.array(rows, dtype=float)
        Sc = np.linalg.svd(C, compute_uv=False)
        rank_C = sum(s > 1e-8 for s in Sc)
        ker_dim = m - rank_C  # dimension of swap cycle space

        if ker_dim == 0:
            continue

        # Build ker basis
        _, _, Vt = np.linalg.svd(C, full_matrices=True)
        ker_basis = Vt[rank_C:]  # rows = ker vectors

        # For each ker vector, build the swap cycle as an A_2 chain
        # and check if it's a boundary of some Omega_3 element
        paths2 = enumerate_allowed_paths(A, n, 2)
        paths3 = enumerate_allowed_paths(A, n, 3)
        paths1 = enumerate_allowed_paths(A, n, 1)
        omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
        dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

        path2_idx = {p: i for i, p in enumerate(paths2)}

        for ki in range(ker_basis.shape[0]):
            M_vec = ker_basis[ki]

            # Build swap cycle in A_2 coords
            z = np.zeros(len(paths2))
            for j, (a, b) in enumerate(arcs_PQ):
                coeff = M_vec[j]
                if abs(coeff) < 1e-12:
                    continue
                # sigma_{a,b}^v = (a,b,v) - (v,a,b)
                if (a,b,v) in path2_idx and (v,a,b) in path2_idx:
                    z[path2_idx[(a,b,v)]] += coeff
                    z[path2_idx[(v,a,b)]] -= coeff

            if np.max(np.abs(z)) < 1e-12:
                continue

            total_tests += 1

            # Check if z is in im(d_3|Omega_3)
            if dim_O3 > 0:
                D3 = build_full_boundary_matrix(paths3, paths2)
                D3_omega = D3 @ omega3
                w, res, _, _ = np.linalg.lstsq(D3_omega, z, rcond=None)
                err = np.max(np.abs(D3_omega @ w - z))

                if err > 1e-6:
                    all_boundaries = False
                    scores = tuple(sorted([sum(row) for row in A]))
                    print(f"  FAIL: T#{bits} v={v} scores={scores}: err={err:.2e}")
                    print(f"    P={P}, Q={Q}, arcs_PQ={arcs_PQ}")
                    print(f"    swap cycle: {z[:10]}...")
            else:
                # No Omega_3 at all
                if np.max(np.abs(z)) > 1e-12:
                    all_boundaries = False
                    scores = tuple(sorted([sum(row) for row in A]))
                    print(f"  FAIL (no O3): T#{bits} v={v} scores={scores}")

print(f"\nn=5: {total_tests} swap cycles tested")
print(f"  All boundaries? {all_boundaries}")


# ============================================================
# PART 2: Same for n=6 (exhaustive, may be slow)
# ============================================================
print(f"\n{'='*70}")
print("n=6: SWAP CYCLE FILLING TEST")
print("=" * 70)

n = 6
total_6 = 1 << (n*(n-1)//2)
t0 = time.time()
total_tests_6 = 0
all_boundaries_6 = True
tested_bits = 0

for bits in range(total_6):
    A = build_adj(n, bits)
    tested_bits += 1

    for v in range(n):
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]
        arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]

        if len(arcs_PQ) < 2:
            continue

        m = len(arcs_PQ)
        rows = []
        for a in P:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if a2 == a: row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)
        for b in Q:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if b2 == b: row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)

        if not rows:
            continue

        C = np.array(rows, dtype=float)
        Sc = np.linalg.svd(C, compute_uv=False)
        rank_C = sum(s > 1e-8 for s in Sc)
        ker_dim = m - rank_C

        if ker_dim == 0:
            continue

        _, _, Vt = np.linalg.svd(C, full_matrices=True)
        ker_basis = Vt[rank_C:]

        paths2 = enumerate_allowed_paths(A, n, 2)
        paths3 = enumerate_allowed_paths(A, n, 3)
        omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
        dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

        path2_idx = {p: i for i, p in enumerate(paths2)}

        for ki in range(min(ker_basis.shape[0], 3)):  # Test up to 3 per vertex
            M_vec = ker_basis[ki]
            z = np.zeros(len(paths2))
            for j, (a, b) in enumerate(arcs_PQ):
                coeff = M_vec[j]
                if abs(coeff) < 1e-12: continue
                if (a,b,v) in path2_idx and (v,a,b) in path2_idx:
                    z[path2_idx[(a,b,v)]] += coeff
                    z[path2_idx[(v,a,b)]] -= coeff

            if np.max(np.abs(z)) < 1e-12:
                continue

            total_tests_6 += 1

            if dim_O3 > 0:
                D3 = build_full_boundary_matrix(paths3, paths2)
                D3_omega = D3 @ omega3
                w, res, _, _ = np.linalg.lstsq(D3_omega, z, rcond=None)
                err = np.max(np.abs(D3_omega @ w - z))

                if err > 1e-6:
                    all_boundaries_6 = False
                    scores = tuple(sorted([sum(row) for row in A]))
                    print(f"  FAIL: T#{bits} v={v} scores={scores}: err={err:.2e}")
            else:
                if np.max(np.abs(z)) > 1e-12:
                    all_boundaries_6 = False
                    print(f"  FAIL (no O3): T#{bits} v={v}")

    if (bits+1) % 5000 == 0:
        elapsed = time.time() - t0
        print(f"  {bits+1}/{total_6} ({elapsed:.0f}s) tests={total_tests_6} all_ok={all_boundaries_6}")

elapsed = time.time() - t0
print(f"\nn=6 ({elapsed:.0f}s): {total_tests_6} swap cycles tested")
print(f"  All boundaries? {all_boundaries_6}")


# ============================================================
# PART 3: Check if there even ARE nonzero swap cycles at n=5
# ============================================================
print(f"\n{'='*70}")
print("SWAP CYCLE EXISTENCE AND DIMENSION")
print("=" * 70)

n = 5
swap_dim_data = defaultdict(list)

for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))

    total_swap_dim = 0
    for v in range(n):
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]
        arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]

        if len(arcs_PQ) < 2:
            continue

        m = len(arcs_PQ)
        rows = []
        for a in P:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if a2 == a: row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)
        for b in Q:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if b2 == b: row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)

        if not rows:
            continue

        C = np.array(rows, dtype=float)
        Sc = np.linalg.svd(C, compute_uv=False)
        rank_C = sum(s > 1e-8 for s in Sc)
        ker_dim = m - rank_C
        total_swap_dim += ker_dim

    swap_dim_data[scores].append(total_swap_dim)

print(f"n=5 swap cycle dimensions by score:")
for scores in sorted(swap_dim_data.keys()):
    dims = swap_dim_data[scores]
    dim_set = sorted(set(dims))
    print(f"  scores={scores}: swap_dims={dim_set} (count={len(dims)})")


print("\n\nDone.")
