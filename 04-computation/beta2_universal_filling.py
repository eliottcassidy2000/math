#!/usr/bin/env python3
"""
beta2_universal_filling.py - Search for a universal filling formula

From the explicit fillings at n=5:
- v appears only at positions 0 and 3 in the filling 3-paths
- The 3-paths are (v, a1, a2, a3) and (a1, a2, a3, v)
  where (a1,a2,a3) are ALLOWED 3-element paths in T' = T\{v}

Hypothesis: The filling for a swap cycle z = sum M_{ab}[sigma_{ab}^v]
uses EXACTLY the allowed 2-paths of T' = T\{v}, coned from both sides.

More precisely, define:
  K_front(z') = sum over allowed 2-paths (a,b,c) in T': alpha * (v,a,b,c)
  K_back(z') = sum over allowed 2-paths (a,b,c) in T': beta * (a,b,c,v)

where alpha, beta depend on the swap cycle z.

Since d_3(v,a,b,c) = (a,b,c) - (v,b,c) + (v,a,c) - (v,a,b)
  and d_3(a,b,c,v) = (b,c,v) - (a,c,v) + (a,b,v) - (a,b,c):

For a swap cycle z = sum M_{ab}[(a,b,v) - (v,a,b)]:
We need d_3(w) to produce (a,b,v) terms and -(v,a,b) terms.

From K_back: d_3(a,b,c,v) produces +(a,b,v) as the (delete c) face.
From K_front: d_3(v,a,b,c) produces -(v,a,b) as the (delete c) face.

So the natural construction is:
  w = sum_{(a,b,c) allowed in T', b->v (c in Q)} gamma_{abc} * (a,b,c,v)
    + sum_{(a,b,c) allowed in T', v->a (a in P)} delta_{abc} * (v,a,b,c)

where we choose gamma, delta to match the swap cycle.

This is essentially the "cone from T' through v" approach.

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time
import numpy as np
from collections import defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


# ============================================================
# PART 1: Analyze the actual fillings at n=5
# ============================================================
print("=" * 70)
print("FILLING ANALYSIS: WHAT 3-PATHS ARE USED?")
print("=" * 70)

n = 5

for bits in [12]:  # The example from before
    A = build_adj(n, bits)
    v = 0
    P = [a for a in range(n) if a != v and A[v][a] == 1]
    Q = [b for b in range(n) if b != v and A[b][v] == 1]

    print(f"T#{bits}, v={v}, P={P}, Q={Q}")

    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0
    path2_idx = {p: i for i, p in enumerate(paths2)}

    # Build swap cycle
    arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]
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

    C = np.array(rows, dtype=float)
    _, _, Vt = np.linalg.svd(C, full_matrices=True)
    Sc = np.linalg.svd(C, compute_uv=False)
    rank_C = sum(s > 1e-8 for s in Sc)
    M_vec = Vt[rank_C]

    z = np.zeros(len(paths2))
    for j, (a, b) in enumerate(arcs_PQ):
        coeff = M_vec[j]
        if abs(coeff) < 1e-12: continue
        if (a,b,v) in path2_idx: z[path2_idx[(a,b,v)]] += coeff
        if (v,a,b) in path2_idx: z[path2_idx[(v,a,b)]] -= coeff

    # Get the actual filling
    D3 = build_full_boundary_matrix(paths3, paths2)
    D3_omega = D3 @ omega3
    w, _, _, _ = np.linalg.lstsq(D3_omega, z, rcond=None)
    w_A3 = omega3 @ w

    # Analyze: which 3-paths in the filling involve v?
    print(f"\n  3-paths in T' (not involving v):")
    # Build T' = T \ {v}
    A_prime = [[A[i][j] for j in range(n) if j != v] for i in range(n) if i != v]
    V_prime = [i for i in range(n) if i != v]
    n_prime = n - 1
    # Allowed 2-paths in T'
    paths2_prime = enumerate_allowed_paths(A_prime, n_prime, 2)
    # Map back to original vertex labels
    paths2_prime_orig = [tuple(V_prime[x] for x in p) for p in paths2_prime]
    print(f"    Allowed 2-paths in T': {paths2_prime_orig}")

    print(f"\n  Filling 3-paths:")
    for j, path in enumerate(paths3):
        if abs(w_A3[j]) > 1e-10:
            # Is this v + (2-path in T') or (2-path in T') + v?
            if path[0] == v:
                rest = path[1:]
                rest_orig = tuple(rest)
                in_T_prime = rest_orig in paths2_prime_orig
                print(f"    {path}: {w_A3[j]:+.6f} = v + {rest_orig} [in T'? {in_T_prime}]")
            elif path[3] == v:
                rest = path[:3]
                rest_orig = tuple(rest)
                in_T_prime = rest_orig in paths2_prime_orig
                print(f"    {path}: {w_A3[j]:+.6f} = {rest_orig} + v [in T'? {in_T_prime}]")
            else:
                print(f"    {path}: {w_A3[j]:+.6f} [v at middle position]")


# ============================================================
# PART 2: The cone construction from T' perspective
# ============================================================
print(f"\n\n{'='*70}")
print("CONE FROM T' THROUGH v: d_3 ANALYSIS")
print("=" * 70)

# d_3(v,a,b,c) = (a,b,c) - (v,b,c) + (v,a,c) - (v,a,b)
#   = [T' part: (a,b,c)] + [v-parts: -(v,b,c) + (v,a,c) - (v,a,b)]
#
# d_3(a,b,c,v) = (b,c,v) - (a,c,v) + (a,b,v) - (a,b,c)
#   = [v-parts: (b,c,v) - (a,c,v) + (a,b,v)] + [T' part: -(a,b,c)]
#
# So if we define:
#   K_f(a,b,c) = (v,a,b,c) [cone from front]
#   K_b(a,b,c) = (a,b,c,v) [cone from back]
#
# Then d_3 K_f(a,b,c) = (a,b,c) + v-stuff
#      d_3 K_b(a,b,c) = -(a,b,c) + v-stuff
#
# So d_3(K_f + K_b)(a,b,c) = 0 + (only v-paths)
#
# The v-paths from K_f: -(v,b,c) + (v,a,c) - (v,a,b)
# The v-paths from K_b: (b,c,v) - (a,c,v) + (a,b,v)
#
# Total v-paths: (a,b,v) + (b,c,v) - (a,c,v) - (v,a,b) + (v,a,c) - (v,b,c)
#
# If we sum over ALL allowed 2-paths (a,b,c) in T' with weights alpha_{abc}:
# d_3(sum alpha * (K_f + K_b)(a,b,c)) = sum alpha * [v-path terms]
#
# For this to equal the swap cycle z = sum M_{ab}[(a,b,v)-(v,a,b)]:
# We need:
#   sum_{(a,b,c): c=fixed} alpha_{abc} * [(a,b,v)] = M_{ab} * (a,b,v)
#   Wait, the (a,b,v) term comes from K_b(a,b,c) face at position delete-c.
#   So: coefficient of (a,b,v) = sum_{c: b->c, c!=v} alpha_{abc}
#   And: coefficient of (v,a,b) = sum_{c: b->c, c!=v} (-alpha_{abc})
#   So: coefficient of [(a,b,v) - (v,a,b)] = 2 * sum_c alpha_{abc}
#   Wait, that's not right either. Let me be more careful.

# From K_b(a,b,c) = (a,b,c,v):
#   d_3(a,b,c,v) = (b,c,v) - (a,c,v) + (a,b,v) - (a,b,c)
#   Coefficient of (a,b,v): +1 (always, regardless of c)
#   Coefficient of (v,a,b): 0 (doesn't appear)

# From K_f(a,b,c) = (v,a,b,c):
#   d_3(v,a,b,c) = (a,b,c) - (v,b,c) + (v,a,c) - (v,a,b)
#   Coefficient of (v,a,b): -1 (always, regardless of c)
#   Coefficient of (a,b,v): 0 (doesn't appear)

# So K_f + K_b on (a,b,c) gives:
#   (a,b,v) - (v,a,b) + (b,c,v) - (a,c,v) + (v,a,c) - (v,b,c) + 0*(a,b,c)
#
# The (a,b,v) - (v,a,b) = sigma_{a,b}^v (the swap!)
# The extra terms: (b,c,v) - (a,c,v) + (v,a,c) - (v,b,c)

# So if I apply K_f + K_b to a SINGLE 2-path (a,b,c) in T':
#   I get the swap sigma_{a,b}^v PLUS extra terms.
# The extra terms involve c and would need to cancel when summing over
# multiple 2-paths with appropriate weights.

# This means: to get the swap cycle z = sum M_{ab} sigma_{ab}^v,
# I need to solve:
#   sum_c alpha_{abc} = M_{ab}  (for each a,b with a->b)
# AND:
#   sum of extra terms = 0 (cancellation condition)

# The extra terms at a 2-path (x,y) (v-type):
# (b,c,v): coefficient of (b,c,v) from K_f+K_b applied to (a,b,c) = +1
# (a,c,v): coefficient of (a,c,v) from K_f+K_b applied to (a,b,c) = -1
# (v,a,c): coefficient of (v,a,c) from K_f+K_b applied to (a,b,c) = +1
# (v,b,c): coefficient of (v,b,c) from K_f+K_b applied to (a,b,c) = -1

# So, for a given v-type 2-path (x,y,v) [x,y != v]:
# It appears as (b,c,v) when (a,b,c) has b=x, c=y
# and as (a,c,v) when (a,b,c) has a=x, c=y
# Total coefficient of (x,y,v) in extra terms:
#   sum_{a: (a,x,y) in T'} alpha_{axy} - sum_{b: (x,b,y) in T'} alpha_{xby}

# For v-type (v,x,y):
# It appears as (v,a,c) when a=x, c=y
# and as (v,b,c) when b=x, c=y
# Total coefficient of (v,x,y) in extra terms:
#   sum_{b: (x,b,y) in T'} alpha_{xby} - sum_{a: (a,x,y) in T'} alpha_{axy}
# = NEGATIVE of the (x,y,v) coefficient

# So if we denote the extra coefficient of (x,y,v) as E(x,y):
#   E(x,y) = sum_{a: (a,x,y) in T'} alpha_{axy} - sum_{b: (x,b,y) in T'} alpha_{xby}
# And the extra coefficient of (v,x,y) = -E(x,y).

# For the extra terms to cancel completely (giving pure swap cycle):
#   E(x,y) = 0 for all (x,y) with x,y != v and x->y
#   AND -E(x,y) = 0, which is the same condition.

# So: E(x,y) = 0 means:
# sum_{a: a->x in T', (a,x,y) allowed} alpha_{axy} = sum_{b: y->b?? No...}

# Wait: (x,b,y) is a 2-path in T' meaning x->b->y in T'.
# sum_{b: (x,b,y) in T'} alpha_{xby} means: sum over all mediators b
# between x and y.

# E(x,y) = [inflow to (x,y) from position 0] - [flow through (x,y) at position 0-1]
# This is like a "divergence" condition on alpha as a function on 2-paths.

# This is actually the d_2(alpha) condition!
# E(x,y) = sum_a alpha_{axy} - sum_b alpha_{xby}
#         = [d_2(alpha)]_{xy} component (part of it)
# More precisely, d_2(alpha) at 1-path (x,y):
# = sum_a alpha_{axy} - sum_b alpha_{xby} + sum_c alpha_{xyc}

# So E(x,y) = [d_2(alpha)]_{xy} - sum_c alpha_{xyc}

# Hmm, this is getting complicated. Let me just test numerically.

# Numerical test: find alpha_{abc} for all allowed 2-paths in T' such that:
# (1) sum_{c: b->c in T', c != v} alpha_{abc} = M_{ab} for each P->Q arc
# (2) Extra terms cancel

print(f"\nNumerical test: cone construction for T#12, v=0")

n = 5
bits = 12
A = build_adj(n, bits)
v = 0

# T' vertices
V_prime = [i for i in range(n) if i != v]
A_prime = [[A[i][j] for j in V_prime] for i in V_prime]
n_prime = len(V_prime)

# Allowed 2-paths in T' (in original vertex labels)
paths2_prime = []
for a_idx in range(n_prime):
    for b_idx in range(n_prime):
        for c_idx in range(n_prime):
            a, b, c = V_prime[a_idx], V_prime[b_idx], V_prime[c_idx]
            if a != b and b != c and a != c:
                if A[a][b] == 1 and A[b][c] == 1:
                    # Check if allowed in T'
                    # This is an allowed 2-path if either:
                    # - a->c (transitive triple in T')
                    # - OR it's part of a 3-cycle, with Omega constraint satisfied
                    paths2_prime.append((a, b, c))

# Actually, I should use the enumerate function properly
paths2_Tp = enumerate_allowed_paths(A_prime, n_prime, 2)
paths2_Tp_orig = [tuple(V_prime[x] for x in p) for p in paths2_Tp]

print(f"  T' vertices: {V_prime}")
print(f"  Allowed 2-paths in T' ({len(paths2_Tp_orig)}): {paths2_Tp_orig}")

paths2 = enumerate_allowed_paths(A, n, 2)
paths3 = enumerate_allowed_paths(A, n, 3)
path2_idx = {p: i for i, p in enumerate(paths2)}
path3_idx = {p: i for i, p in enumerate(paths3)}

P = [a for a in range(n) if a != v and A[v][a] == 1]
Q = [b for b in range(n) if b != v and A[b][v] == 1]
arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]

# Build the system: for each T' 2-path (a,b,c), the variable is alpha_{abc}.
# The 3-chain w = sum alpha_{abc} * [(v,a,b,c) + (a,b,c,v)]
# and we need d_3(w) = swap cycle z.

# Variables: alpha[j] for each j-th T' 2-path
m = len(paths2_Tp_orig)

# Build the boundary map from these 3-chains to paths2
# For each (a,b,c) in T': (v,a,b,c) and (a,b,c,v) are potential 3-paths
# d_3(v,a,b,c) = (a,b,c) - (v,b,c) + (v,a,c) - (v,a,b)
# d_3(a,b,c,v) = (b,c,v) - (a,c,v) + (a,b,v) - (a,b,c)

# Let's build the combined boundary matrix
B = np.zeros((len(paths2), m))

for j, (a, b, c) in enumerate(paths2_Tp_orig):
    # d_3(v,a,b,c) contributions:
    for face, sign in [((a,b,c), +1), ((v,b,c), -1), ((v,a,c), +1), ((v,a,b), -1)]:
        if face in path2_idx:
            B[path2_idx[face], j] += sign

    # d_3(a,b,c,v) contributions:
    for face, sign in [((b,c,v), +1), ((a,c,v), -1), ((a,b,v), +1), ((a,b,c), -1)]:
        if face in path2_idx:
            B[path2_idx[face], j] += sign

# The (a,b,c) terms cancel: +1 from front cone, -1 from back cone.
# So B should only have v-path entries.

print(f"\n  Boundary matrix B shape: {B.shape}")
print(f"  B nonzero pattern check (T' paths should cancel):")
for j, (a,b,c) in enumerate(paths2_Tp_orig):
    if (a,b,c) in path2_idx and abs(B[path2_idx[(a,b,c)], j]) > 1e-10:
        print(f"    T' path ({a},{b},{c}) has nonzero coeff in column {j}: {B[path2_idx[(a,b,c)], j]}")
print(f"  (Empty = good, means T' parts cancel)")

# Build the swap cycle z (same as before)
_, _, Vt_C = np.linalg.svd(np.array([[1 if a2==a else 0 for j, (a2, b2) in enumerate(arcs_PQ)]
                                     for a in P] +
                                    [[1 if b2==b else 0 for j, (a2, b2) in enumerate(arcs_PQ)]
                                     for b in Q], dtype=float), full_matrices=True)
Sc_C = np.linalg.svd(np.array([[1 if a2==a else 0 for j, (a2, b2) in enumerate(arcs_PQ)]
                                for a in P] +
                               [[1 if b2==b else 0 for j, (a2, b2) in enumerate(arcs_PQ)]
                                for b in Q], dtype=float), compute_uv=False)
rank_C = sum(s > 1e-8 for s in Sc_C)
M_vec = Vt_C[rank_C]

z = np.zeros(len(paths2))
for j, (a, b) in enumerate(arcs_PQ):
    coeff = M_vec[j]
    if abs(coeff) < 1e-12: continue
    if (a,b,v) in path2_idx: z[path2_idx[(a,b,v)]] += coeff
    if (v,a,b) in path2_idx: z[path2_idx[(v,a,b)]] -= coeff

# Solve B @ alpha = z
alpha, res, _, _ = np.linalg.lstsq(B, z, rcond=None)
err = np.max(np.abs(B @ alpha - z))

print(f"\n  Solve B @ alpha = z: err = {err:.2e}")

if err < 1e-6:
    print(f"  SUCCESS! The cone-from-T' construction works!")
    print(f"  Alpha coefficients:")
    for j, (a,b,c) in enumerate(paths2_Tp_orig):
        if abs(alpha[j]) > 1e-10:
            print(f"    alpha[({a},{b},{c})] = {alpha[j]:+.6f}")

    # Check: is the resulting 3-chain in Omega_3?
    w_3 = np.zeros(len(paths3))
    for j, (a,b,c) in enumerate(paths2_Tp_orig):
        if abs(alpha[j]) < 1e-12: continue
        path_f = (v,a,b,c)
        path_b = (a,b,c,v)
        if path_f in path3_idx:
            w_3[path3_idx[path_f]] += alpha[j]
        if path_b in path3_idx:
            w_3[path3_idx[path_b]] += alpha[j]

    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    proj3, _, _, _ = np.linalg.lstsq(omega3, w_3, rcond=None)
    w_3_proj = omega3 @ proj3
    omega3_err = np.max(np.abs(w_3 - w_3_proj))
    print(f"  w in Omega_3? err = {omega3_err:.2e}")
else:
    print(f"  FAIL: err = {err:.2e}")
    print(f"  The cone-from-T' approach does not directly work.")
    print(f"  Rank of B: {np.linalg.matrix_rank(B, tol=1e-8)}")
    print(f"  dim(im B) = {np.linalg.matrix_rank(B, tol=1e-8)}, dim target = {sum(abs(z[i]) > 1e-10 for i in range(len(z)))}")


# ============================================================
# PART 3: Test for ALL n=5 tournaments
# ============================================================
print(f"\n{'='*70}")
print("CONE-FROM-T' TEST FOR ALL n=5")
print("=" * 70)

success = 0
fail = 0
total_tested = 0

for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)
    path2_idx_local = {p: i for i, p in enumerate(paths2)}
    path3_idx_local = {p: i for i, p in enumerate(paths3)}

    for v in range(n):
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]
        arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]

        if len(arcs_PQ) < 2:
            continue

        m_arc = len(arcs_PQ)
        rows = []
        for a in P:
            row = [0] * m_arc
            for j, (a2, b2) in enumerate(arcs_PQ):
                if a2 == a: row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)
        for b in Q:
            row = [0] * m_arc
            for j, (a2, b2) in enumerate(arcs_PQ):
                if b2 == b: row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)

        if not rows:
            continue

        C_mat = np.array(rows, dtype=float)
        Sc = np.linalg.svd(C_mat, compute_uv=False)
        rank_C = sum(s > 1e-8 for s in Sc)
        _, _, Vt = np.linalg.svd(C_mat, full_matrices=True)
        ker_dim = m_arc - rank_C

        if ker_dim == 0:
            continue

        # T' paths
        V_prime = [i for i in range(n) if i != v]
        A_prime = [[A[i][j] for j in V_prime] for i in V_prime]
        n_prime = len(V_prime)
        paths2_Tp = enumerate_allowed_paths(A_prime, n_prime, 2)
        paths2_Tp_orig = [tuple(V_prime[x] for x in p) for p in paths2_Tp]

        # Build boundary matrix B
        m_Tp = len(paths2_Tp_orig)
        B_local = np.zeros((len(paths2), m_Tp))
        for j, (a, b, c) in enumerate(paths2_Tp_orig):
            for face, sign in [((a,b,c), +1), ((v,b,c), -1), ((v,a,c), +1), ((v,a,b), -1)]:
                if face in path2_idx_local:
                    B_local[path2_idx_local[face], j] += sign
            for face, sign in [((b,c,v), +1), ((a,c,v), -1), ((a,b,v), +1), ((a,b,c), -1)]:
                if face in path2_idx_local:
                    B_local[path2_idx_local[face], j] += sign

        for ki in range(ker_dim):
            M_vec = Vt[rank_C + ki]
            z = np.zeros(len(paths2))
            for j, (a, b) in enumerate(arcs_PQ):
                coeff = M_vec[j]
                if abs(coeff) < 1e-12: continue
                if (a,b,v) in path2_idx_local: z[path2_idx_local[(a,b,v)]] += coeff
                if (v,a,b) in path2_idx_local: z[path2_idx_local[(v,a,b)]] -= coeff

            if np.max(np.abs(z)) < 1e-12:
                continue

            total_tested += 1

            alpha, _, _, _ = np.linalg.lstsq(B_local, z, rcond=None)
            err = np.max(np.abs(B_local @ alpha - z))

            if err < 1e-6:
                success += 1
            else:
                fail += 1
                if fail <= 3:
                    scores = tuple(sorted([sum(row) for row in A]))
                    print(f"  FAIL: T#{bits} v={v} scores={scores} err={err:.2e}")

print(f"\nn=5: {total_tested} swap cycles, {success} success, {fail} fail")
print(f"  Cone-from-T' works for ALL? {'YES' if fail == 0 else 'NO'}")


print("\n\nDone.")
