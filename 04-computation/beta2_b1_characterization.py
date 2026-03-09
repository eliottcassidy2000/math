#!/usr/bin/env python3
"""beta2_b1_characterization.py - Characterize b1=0 vs b1=1

Goal: Find an EXPLICIT characterization of when b1(T)=1 vs b1(T)=0.

KEY FACTS:
- b1(T) in {0,1} for all tournaments (verified n<=20)
- b1 = dim(Z_1) - rk(d_2) = (n-1)(n-2)/2 - rk(d_2)
- b1=1 iff rk(d_2) = (n-1)(n-2)/2 - 1
- When b1=1: the cokernel of d_2|Z_1 is 1-dimensional

PLAN:
1. Compute the H_1 generator when b1=1
2. Look for patterns in which tournaments have b1=1
3. Try to express b1 in terms of tournament invariants
4. Understanding this is KEY to proving HYP-278

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


def compute_b1_detailed(A, n):
    """Compute b1 and return the H1 generator if b1=1."""
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths0 = [(i,) for i in range(n)]
    paths2 = enumerate_allowed_paths(A, n, 2)

    if not paths1:
        return 0, None, None

    omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
    dim_O1 = omega1.shape[1] if omega1.ndim == 2 else 0
    if dim_O1 == 0:
        return 0, None, None

    D1 = build_full_boundary_matrix([tuple(p) for p in paths1], paths0)
    D1_om = D1 @ omega1
    sv = np.linalg.svd(D1_om, compute_uv=False)
    rk_d1 = int(sum(s > 1e-8 for s in sv))

    # Compute Z1 basis (kernel of d1 restricted to Omega1)
    U, S, Vt = np.linalg.svd(D1_om, full_matrices=True)
    nullity = dim_O1 - rk_d1
    Z1_basis_omega = Vt[rk_d1:].T  # in Omega1 coordinates
    Z1_basis = omega1 @ Z1_basis_omega  # in path coordinates
    dim_Z1 = Z1_basis.shape[1]

    if paths2:
        omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
        dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
        if dim_O2 > 0:
            D2 = build_full_boundary_matrix([tuple(p) for p in paths2],
                                            [tuple(p) for p in paths1])
            D2_om = D2 @ omega2
            sv2 = np.linalg.svd(D2_om, compute_uv=False)
            rk_d2 = int(sum(s > 1e-8 for s in sv2))

            # B1 = im(d2) in path coordinates
            B1 = D2_om[:, :rk_d2] if rk_d2 > 0 else None  # don't need this

            # For H1 generator: project Z1 onto complement of im(d2)
            # H1 = Z1 / B1, so we need a representative of the coset
            if rk_d2 > 0:
                # im(d2) in path coordinates
                im_d2 = D2 @ omega2  # columns are boundary images
                # Project Z1 basis onto orthogonal complement of im(d2)
                U2, S2, _ = np.linalg.svd(im_d2, full_matrices=True)
                proj = np.eye(len(paths1)) - U2[:, :rk_d2] @ U2[:, :rk_d2].T
                Z1_proj = proj @ Z1_basis
                # Find nonzero columns
                norms = np.linalg.norm(Z1_proj, axis=0)
                h1_gens = Z1_proj[:, norms > 1e-8]
            else:
                h1_gens = Z1_basis
                rk_d2 = 0
        else:
            rk_d2 = 0
            h1_gens = Z1_basis
    else:
        rk_d2 = 0
        h1_gens = Z1_basis

    b1 = dim_Z1 - rk_d2
    path1_list = [tuple(p) for p in paths1]

    return b1, h1_gens, path1_list


# ============================================================
# Part 1: H1 generators at n=5
# ============================================================
print("=" * 70)
print("H1 GENERATORS WHEN b1=1, n=5")
print("=" * 70)

n = 5
total = 2 ** (n * (n - 1) // 2)
b1_1_examples = []

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

    b1, h1_gens, paths1 = compute_b1_detailed(A, n)
    if b1 == 1 and len(b1_1_examples) < 5:
        # Get the H1 generator
        gen = h1_gens[:, 0]
        # Normalize
        gen = gen / gen[np.argmax(np.abs(gen))]
        # Which arcs have nonzero coefficient?
        nonzero = [(paths1[i], gen[i]) for i in range(len(paths1)) if abs(gen[i]) > 1e-8]
        scores = tuple(sorted(sum(A[i]) for i in range(n)))
        c3 = 0
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    if A[i][j] and A[j][k] and A[k][i]:
                        c3 += 1
                    if A[i][k] and A[k][j] and A[j][i]:
                        c3 += 1

        b1_1_examples.append({
            'bits': bits, 'scores': scores, 'c3': c3,
            'generator': nonzero, 'A': [row[:] for row in A]
        })

print(f"First {len(b1_1_examples)} b1=1 tournaments at n=5:\n")
for ex in b1_1_examples:
    print(f"  bits={ex['bits']}, scores={ex['scores']}, c3={ex['c3']}")
    print(f"  H1 generator:")
    for arc, coeff in ex['generator']:
        print(f"    {arc}: {coeff:.4f}")
    print()


# ============================================================
# Part 2: b1 vs specific tournament invariants
# ============================================================
print("=" * 70)
print("b1 vs TOURNAMENT INVARIANTS (n=5,6)")
print("=" * 70)

for n in [5, 6]:
    total = 2 ** (n * (n - 1) // 2)
    # Track: b1, score, c3, and also:
    # - number of vertices with d_out >= (n-1)/2
    # - number of Hamiltonian cycles
    # - number of strong components
    data = {0: [], 1: []}

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

        b1, _, _ = compute_b1_detailed(A, n)
        scores = tuple(sorted(sum(A[i]) for i in range(n)))
        c3 = 0
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    if A[i][j] and A[j][k] and A[k][i]:
                        c3 += 1
                    if A[i][k] and A[k][j] and A[j][i]:
                        c3 += 1

        # Is the tournament strongly connected?
        # Simple BFS from 0
        def is_strongly_connected(A, n):
            for start in [0]:
                visited = {start}
                stack = [start]
                while stack:
                    u = stack.pop()
                    for v in range(n):
                        if A[u][v] and v not in visited:
                            visited.add(v)
                            stack.append(v)
                if len(visited) < n:
                    return False
                # Check reverse
                visited = {start}
                stack = [start]
                while stack:
                    u = stack.pop()
                    for v in range(n):
                        if A[v][u] and v not in visited:
                            visited.add(v)
                            stack.append(v)
                if len(visited) < n:
                    return False
            return True

        sc = is_strongly_connected(A, n)

        data[b1].append({
            'scores': scores, 'c3': c3, 'sc': sc
        })

    print(f"\nn={n}:")
    print(f"  b1=0: {len(data[0])}, b1=1: {len(data[1])}")

    # Is strong connectivity related to b1?
    sc_0 = sum(1 for d in data[0] if d['sc'])
    sc_1 = sum(1 for d in data[1] if d['sc'])
    print(f"  Strongly connected: b1=0: {sc_0}/{len(data[0])}, b1=1: {sc_1}/{len(data[1])}")

    # Check: is b1=1 equivalent to some simple condition?
    # Hypothesis: b1=1 iff strongly connected?
    print(f"  b1=1 and NOT SC: {sum(1 for d in data[1] if not d['sc'])}")
    print(f"  b1=0 and SC: {sum(1 for d in data[0] if d['sc'])}")

    # Score distribution
    score_b1 = {0: Counter(), 1: Counter()}
    for b in [0, 1]:
        for d in data[b]:
            score_b1[b][d['scores']] += 1

    print(f"\n  Scores that appear in BOTH b1=0 and b1=1:")
    both_scores = set(score_b1[0].keys()) & set(score_b1[1].keys())
    for s in sorted(both_scores):
        print(f"    {s}: b1=0: {score_b1[0][s]}, b1=1: {score_b1[1][s]}")

    print(f"\n  Scores ONLY in b1=0:")
    only_0 = set(score_b1[0].keys()) - set(score_b1[1].keys())
    for s in sorted(only_0):
        print(f"    {s}: {score_b1[0][s]}")

    print(f"\n  Scores ONLY in b1=1:")
    only_1 = set(score_b1[1].keys()) - set(score_b1[0].keys())
    for s in sorted(only_1):
        print(f"    {s}: {score_b1[1][s]}")


# ============================================================
# Part 3: Omega_1 and d_2 structure
# ============================================================
print(f"\n{'=' * 70}")
print("OMEGA DIMENSIONS AND RANK")
print("=" * 70)

n = 5
total = 2 ** (n * (n - 1) // 2)
dim_omega2_by_b1 = {0: Counter(), 1: Counter()}
rk_d2_by_b1 = {0: Counter(), 1: Counter()}

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

    omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
    dim_O1 = omega1.shape[1] if omega1.ndim == 2 else 0

    if paths2:
        omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
        dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
    else:
        dim_O2 = 0
        omega2 = np.zeros((0, 0))

    # Compute b1 from rank
    if dim_O1 == 0:
        b1 = 0
        rk_d2 = 0
    else:
        D1 = build_full_boundary_matrix([tuple(p) for p in paths1], paths0)
        D1_om = D1 @ omega1
        sv = np.linalg.svd(D1_om, compute_uv=False)
        rk_d1 = int(sum(s > 1e-8 for s in sv))
        dim_Z1 = dim_O1 - rk_d1

        if dim_O2 > 0:
            D2 = build_full_boundary_matrix([tuple(p) for p in paths2],
                                            [tuple(p) for p in paths1])
            D2_om = D2 @ omega2
            sv2 = np.linalg.svd(D2_om, compute_uv=False)
            rk_d2 = int(sum(s > 1e-8 for s in sv2))
        else:
            rk_d2 = 0

        b1 = dim_Z1 - rk_d2

    dim_omega2_by_b1[b1][dim_O2] += 1
    rk_d2_by_b1[b1][rk_d2] += 1

print(f"\nn=5:")
print(f"  dim(Omega_2) by b1:")
for b in [0, 1]:
    for d in sorted(dim_omega2_by_b1[b].keys()):
        print(f"    b1={b}, dim(Omega_2)={d}: {dim_omega2_by_b1[b][d]}")

print(f"\n  rk(d_2) by b1:")
for b in [0, 1]:
    for r in sorted(rk_d2_by_b1[b].keys()):
        print(f"    b1={b}, rk(d_2)={r}: {rk_d2_by_b1[b][r]}")


# ============================================================
# Part 4: b1 and strongly connected components
# ============================================================
print(f"\n{'=' * 70}")
print("b1=1 IFF STRONGLY CONNECTED? (n=5,6)")
print("=" * 70)

# Actually, b1 in graph homology relates to cycles.
# For path homology of a digraph, b1 counts independent 1-cycles
# not killed by 2-boundaries. For a tournament, what's the
# relationship to strong connectivity?

# Let's check: does b1=1 require strongly connected?
# From part 2, we'll see. Let me also check: for b1=0 tournaments
# with ALL sub-tournaments having b1=1, are they special?

# Check b1 of all n=4 tournaments
print("\nb1 at n=4:")
n = 4
total = 2 ** (n * (n - 1) // 2)
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
    b1, _, _ = compute_b1_detailed(A, n)
    if b1 > 0:
        scores = tuple(sorted(sum(A[i]) for i in range(n)))
        c3 = sum(1 for i in range(n) for j in range(i+1, n) for k in range(j+1, n)
                 if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]))
        print(f"  bits={bits}, scores={scores}, c3={c3}")

# Check b1 of all n=3 tournaments
print("\nb1 at n=3:")
n = 3
total = 2 ** (n * (n - 1) // 2)
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
    b1, _, _ = compute_b1_detailed(A, n)
    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    c3 = sum(1 for i in range(n) for j in range(i+1, n) for k in range(j+1, n)
             if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]))
    print(f"  bits={bits}, scores={scores}, c3={c3}, b1={b1}")


# ============================================================
# Part 5: Can b1 be expressed as a function of the score sequence?
# ============================================================
print(f"\n{'=' * 70}")
print("IS b1 DETERMINED BY SCORE SEQUENCE?")
print("=" * 70)

for n in [5, 6]:
    total = 2 ** (n * (n - 1) // 2)
    score_b1_values = {}

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

        b1, _, _ = compute_b1_detailed(A, n)
        scores = tuple(sorted(sum(A[i]) for i in range(n)))
        if scores not in score_b1_values:
            score_b1_values[scores] = set()
        score_b1_values[scores].add(b1)

    print(f"\nn={n}:")
    mixed = 0
    for s in sorted(score_b1_values.keys()):
        vals = score_b1_values[s]
        if len(vals) > 1:
            mixed += 1
            print(f"  MIXED: {s} -> {vals}")

    if mixed == 0:
        print(f"  b1 IS determined by score sequence!")
        print(f"  Score -> b1 mapping:")
        for s in sorted(score_b1_values.keys()):
            print(f"    {s} -> {list(score_b1_values[s])[0]}")
    else:
        print(f"  b1 is NOT determined by score sequence ({mixed} mixed scores)")


print("\n\nDone.")
