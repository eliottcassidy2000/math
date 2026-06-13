#!/usr/bin/env python3
"""beta2_codimension_analysis.py - Analyze codimension structure of W_v = Z_1(T\v)

KEY ALGEBRAIC FACTS:
- V = Z_1(T) has dimension (n-1)(n-2)/2
- W_v = Z_1(T\v) has dimension (n-2)(n-3)/2
- W_v + W_w should have codimension 1 in V
- The 1-dim complement L_{vw} is a cycle visiting both v and w

PLAN:
1. Verify codimension formula: dim(W_v + W_w) = (n-1)(n-2)/2 - 1
2. Compute L_{vw} explicitly
3. How many vertices does L_{vw} visit?
4. Check: W_v + W_w + W_x = V iff L_{vw} visits x
5. Does this constrain the number of bad vertices?

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random, time
import numpy as np
from collections import Counter
from itertools import combinations
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


def compute_Z1_basis(A, n):
    """Compute basis for Z_1(T) in the full arc coordinate space.
    Returns: Z1_basis (matrix whose columns span Z_1), paths1 list."""
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths0 = [(i,) for i in range(n)]
    path1_list = [tuple(p) for p in paths1]

    if not paths1:
        return None, path1_list

    omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
    dim_O1 = omega1.shape[1] if omega1.ndim == 2 else 0
    if dim_O1 == 0:
        return None, path1_list

    D1 = build_full_boundary_matrix(path1_list, paths0)
    D1_om = D1 @ omega1

    U, S, Vt = np.linalg.svd(D1_om, full_matrices=True)
    rk_d1 = int(sum(s > 1e-8 for s in S))
    Z1_omega = Vt[rk_d1:].T
    Z1 = omega1 @ Z1_omega

    return Z1, path1_list


def embed_Z1_subgraph(Z1_sub, paths1_sub, paths1_full):
    """Embed Z1 of a subgraph into the full arc space."""
    m = len(paths1_full)
    k = Z1_sub.shape[1]
    Z1_emb = np.zeros((m, k))
    for i, p in enumerate(paths1_sub):
        if p in paths1_full:
            j = paths1_full.index(p)
            Z1_emb[j] = Z1_sub[i]
    return Z1_emb


# ============================================================
# Part 1: Verify codimension formula
# ============================================================
print("=" * 70)
print("CODIMENSION OF W_v + W_w IN V")
print("=" * 70)

for n in [5, 6, 7]:
    if n <= 6:
        num_tests = min(200, 2 ** (n * (n - 1) // 2))
        test_all = (n <= 5)
    else:
        num_tests = 100
        test_all = False

    codim_dist = Counter()
    random.seed(42)
    t0 = time.time()

    for trial in range(num_tests):
        if test_all:
            bits = trial
            A = [[0] * n for _ in range(n)]
            idx = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if (bits >> idx) & 1:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                    idx += 1
        else:
            A = random_tournament(n)

        Z1, paths1 = compute_Z1_basis(A, n)
        if Z1 is None:
            continue
        dim_V = Z1.shape[1]

        # Pick a random pair of vertices
        for v in range(min(3, n)):
            for w in range(v + 1, min(v + 3, n)):
                # Compute W_v = Z_1(T\v)
                others_v = [x for x in range(n) if x != v]
                B_v, vlist_v = get_induced(A, n, others_v)
                Z1_v_sub, paths1_v = compute_Z1_basis(B_v, len(others_v))
                if Z1_v_sub is None:
                    continue
                # Relabel paths to original vertices
                paths1_v_orig = [tuple(vlist_v[x] for x in p) for p in paths1_v]
                W_v = embed_Z1_subgraph(Z1_v_sub, paths1_v_orig, paths1)

                # Compute W_w = Z_1(T\w)
                others_w = [x for x in range(n) if x != w]
                B_w, vlist_w = get_induced(A, n, others_w)
                Z1_w_sub, paths1_w = compute_Z1_basis(B_w, len(others_w))
                if Z1_w_sub is None:
                    continue
                paths1_w_orig = [tuple(vlist_w[x] for x in p) for p in paths1_w]
                W_w = embed_Z1_subgraph(Z1_w_sub, paths1_w_orig, paths1)

                # dim(W_v + W_w)
                combined = np.column_stack([W_v, W_w])
                sv = np.linalg.svd(combined, compute_uv=False)
                dim_sum = int(sum(s > 1e-8 for s in sv))
                codim = dim_V - dim_sum
                codim_dist[codim] += 1

        if (trial + 1) % 50 == 0 and not test_all:
            print(f"  n={n}: {trial+1}/{num_tests} ({time.time()-t0:.0f}s)")

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s):")
    print(f"  dim V = (n-1)(n-2)/2 = {(n-1)*(n-2)//2}")
    print(f"  dim W_v = (n-2)(n-3)/2 = {(n-2)*(n-3)//2}")
    print(f"  Expected codim(W_v+W_w) = 1")
    for cd in sorted(codim_dist.keys()):
        print(f"  codim={cd}: {codim_dist[cd]}")


# ============================================================
# Part 2: The complement L_{vw} and its vertex support
# ============================================================
print(f"\n{'=' * 70}")
print("L_{{vw}} VERTEX SUPPORT")
print("=" * 70)

n = 5
total = 2 ** (n * (n - 1) // 2)

support_sizes = Counter()
examples_shown = 0

for bits in range(min(total, 200)):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    Z1, paths1 = compute_Z1_basis(A, n)
    if Z1 is None:
        continue
    dim_V = Z1.shape[1]

    for v in range(n):
        for w in range(v + 1, n):
            others_v = [x for x in range(n) if x != v]
            B_v, vlist_v = get_induced(A, n, others_v)
            Z1_v_sub, paths1_v = compute_Z1_basis(B_v, len(others_v))
            if Z1_v_sub is None:
                continue
            paths1_v_orig = [tuple(vlist_v[x] for x in p) for p in paths1_v]
            W_v = embed_Z1_subgraph(Z1_v_sub, paths1_v_orig, paths1)

            others_w = [x for x in range(n) if x != w]
            B_w, vlist_w = get_induced(A, n, others_w)
            Z1_w_sub, paths1_w = compute_Z1_basis(B_w, len(others_w))
            if Z1_w_sub is None:
                continue
            paths1_w_orig = [tuple(vlist_w[x] for x in p) for p in paths1_w]
            W_w = embed_Z1_subgraph(Z1_w_sub, paths1_w_orig, paths1)

            combined = np.column_stack([W_v, W_w])
            sv_comb = np.linalg.svd(combined, compute_uv=False)
            dim_sum = int(sum(s > 1e-8 for s in sv_comb))

            if dim_V - dim_sum != 1:
                continue

            # Find L_{vw}: orthogonal complement of W_v + W_w in V
            # Project V onto orthogonal complement of (W_v + W_w)
            # V is in arc coordinates, W_v + W_w columns span the subspace
            U_comb, S_comb, _ = np.linalg.svd(combined, full_matrices=True)
            # im = U_comb[:, :dim_sum]
            # The projection of Z1 onto complement of im
            proj = np.eye(len(paths1)) - U_comb[:, :dim_sum] @ U_comb[:, :dim_sum].T
            Z1_proj = proj @ Z1
            norms = np.linalg.norm(Z1_proj, axis=0)
            max_idx = np.argmax(norms)
            if norms[max_idx] < 1e-8:
                continue
            L_vw = Z1_proj[:, max_idx]
            L_vw = L_vw / L_vw[np.argmax(np.abs(L_vw))]

            # Which vertices does L_vw visit? (which arcs have nonzero coeff?)
            visited = set()
            for i, p in enumerate(paths1):
                if abs(L_vw[i]) > 1e-8:
                    visited.add(p[0])
                    visited.add(p[1])

            support_sizes[len(visited)] += 1

            if examples_shown < 5 and len(visited) <= 4:
                examples_shown += 1
                print(f"\n  L_{{{v},{w}}}: visited={visited}")
                for i, p in enumerate(paths1):
                    if abs(L_vw[i]) > 1e-8:
                        print(f"    {p}: {L_vw[i]:.4f}")

print(f"\nL_{{vw}} vertex support sizes:")
for s in sorted(support_sizes.keys()):
    print(f"  {s} vertices: {support_sizes[s]}")


# ============================================================
# Part 3: How many vertices x have W_v + W_w + W_x = V?
# ============================================================
print(f"\n{'=' * 70}")
print("W_v + W_w + W_x = V? HOW MANY x?")
print("=" * 70)

n = 5
num_fills = Counter()
t0 = time.time()

for bits in range(min(total, 200)):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    Z1, paths1 = compute_Z1_basis(A, n)
    if Z1 is None:
        continue
    dim_V = Z1.shape[1]

    # Precompute all W_x
    W = {}
    for x in range(n):
        others = [x2 for x2 in range(n) if x2 != x]
        B, vl = get_induced(A, n, others)
        Z1_sub, p1_sub = compute_Z1_basis(B, len(others))
        if Z1_sub is None:
            W[x] = np.zeros((len(paths1), 0))
            continue
        p1_orig = [tuple(vl[xx] for xx in p) for p in p1_sub]
        W[x] = embed_Z1_subgraph(Z1_sub, p1_orig, paths1)

    for v in range(n):
        for w in range(v + 1, n):
            fills_V = []
            for x in range(n):
                if x == v or x == w:
                    continue
                combined_3 = np.column_stack([W[v], W[w], W[x]])
                sv3 = np.linalg.svd(combined_3, compute_uv=False)
                dim3 = int(sum(s > 1e-8 for s in sv3))
                if dim3 == dim_V:
                    fills_V.append(x)
            num_fills[len(fills_V)] += 1

elapsed = time.time() - t0
print(f"\nn={n} ({elapsed:.0f}s):")
print(f"  #x with W_v+W_w+W_x = V, distribution:")
for k in sorted(num_fills.keys()):
    print(f"    {k} vertices fill: {num_fills[k]}")


# ============================================================
# Part 4: Triple codimension
# ============================================================
print(f"\n{'=' * 70}")
print("CODIM(W_v + W_w + W_x) IN V")
print("=" * 70)

codim3_dist = Counter()

for bits in range(min(total, 100)):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    Z1, paths1 = compute_Z1_basis(A, n)
    if Z1 is None:
        continue
    dim_V = Z1.shape[1]

    W = {}
    for x in range(n):
        others = [x2 for x2 in range(n) if x2 != x]
        B, vl = get_induced(A, n, others)
        Z1_sub, p1_sub = compute_Z1_basis(B, len(others))
        if Z1_sub is None:
            W[x] = np.zeros((len(paths1), 0))
            continue
        p1_orig = [tuple(vl[xx] for xx in p) for p in p1_sub]
        W[x] = embed_Z1_subgraph(Z1_sub, p1_orig, paths1)

    for v, w, x in combinations(range(n), 3):
        combined = np.column_stack([W[v], W[w], W[x]])
        sv = np.linalg.svd(combined, compute_uv=False)
        dim_sum = int(sum(s > 1e-8 for s in sv))
        codim = dim_V - dim_sum
        codim3_dist[codim] += 1

print(f"\nn={n}: codim(W_v + W_w + W_x) distribution:")
for cd in sorted(codim3_dist.keys()):
    print(f"  codim={cd}: {codim3_dist[cd]}")


print("\n\nDone.")
