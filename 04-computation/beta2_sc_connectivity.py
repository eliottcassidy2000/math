#!/usr/bin/env python3
"""beta2_sc_connectivity.py - Strong connectivity and b1 relationship

KEY INSIGHT: b1=1 implies SC (strong connectivity). NOT-SC implies b1=0.

PROOF STRATEGY:
1. If T is NOT SC: delete vertex from boundary of condensation -> not SC -> b1=0.
   THIS CASE IS SOLVED.
2. If T IS SC with b1=0: harder.
   Subcase 2a: kappa(T) = 1 (vertex connectivity 1).
     Delete the cut vertex -> NOT SC -> b1=0. SOLVED.
   Subcase 2b: kappa(T) >= 2 (2-connected, all T\v are SC).
     Need: exists v with b1(T\v)=0 even though T\v is SC.

THIS SCRIPT: Check which case applies at n=5,6.
Is case 2b ever needed? Or does every SC+b1=0 tournament have kappa=1?

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
    """Check if tournament A on n vertices is strongly connected."""
    if n <= 1:
        return True
    # BFS forward from 0
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
    # BFS backward from 0
    visited = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in range(n):
            if A[v][u] and v not in visited:
                visited.add(v)
                stack.append(v)
    return len(visited) == n


def vertex_connectivity(A, n):
    """Compute vertex connectivity of tournament.
    Returns the minimum number of vertices whose removal disconnects it.
    If the tournament is not SC, returns 0."""
    if not is_sc(A, n):
        return 0
    if n <= 2:
        return n - 1
    # Try removing subsets of increasing size
    for k in range(1, n):
        # Try all subsets of size k
        from itertools import combinations
        for subset in combinations(range(n), k):
            remaining = [v for v in range(n) if v not in subset]
            B, _ = get_induced(A, n, remaining)
            if not is_sc(B, len(remaining)):
                return k
    return n - 1  # complete tournament


def compute_b1(A, n):
    """Compute beta_1 of tournament."""
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


# ============================================================
# Part 1: Connectivity of SC+b1=0 tournaments
# ============================================================
print("=" * 70)
print("VERTEX CONNECTIVITY OF SC + b1=0 TOURNAMENTS")
print("=" * 70)

for n in [5, 6]:
    total = 2 ** (n * (n - 1) // 2)
    kappa_dist = Counter()
    sc_b1_0_count = 0
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

        if not is_sc(A, n):
            continue
        b1 = compute_b1(A, n)
        if b1 != 0:
            continue

        sc_b1_0_count += 1
        kappa = vertex_connectivity(A, n)
        kappa_dist[kappa] += 1

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): {sc_b1_0_count} SC + b1=0 tournaments")
    for k in sorted(kappa_dist.keys()):
        print(f"  kappa={k}: {kappa_dist[k]} ({kappa_dist[k]/sc_b1_0_count*100:.1f}%)")

    if 1 in kappa_dist:
        print(f"  Case 2a (kappa=1, cut vertex exists): {kappa_dist[1]}")
    high_kappa = sum(v for k, v in kappa_dist.items() if k >= 2)
    print(f"  Case 2b (kappa>=2, NO cut vertex): {high_kappa}")


# ============================================================
# Part 2: For kappa >= 2 + SC + b1=0: does good vertex exist?
# ============================================================
print(f"\n{'=' * 70}")
print("KAPPA >= 2 + SC + b1=0: GOOD VERTEX EXISTS?")
print("=" * 70)

n = 5
total = 2 ** (n * (n - 1) // 2)
checked = 0
all_good = 0
no_good = 0

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

    if not is_sc(A, n):
        continue
    b1 = compute_b1(A, n)
    if b1 != 0:
        continue

    kappa = vertex_connectivity(A, n)
    if kappa < 2:
        continue

    checked += 1
    has_good = False
    for v in range(n):
        others = [x for x in range(n) if x != v]
        B, _ = get_induced(A, n, others)
        if compute_b1(B, n - 1) == 0:
            has_good = True
            break
    if has_good:
        all_good += 1
    else:
        no_good += 1
        print(f"  COUNTEREXAMPLE: bits={bits}")

print(f"\nn={n}: {checked} SC + b1=0 + kappa>=2 tournaments")
print(f"  Has good vertex: {all_good}/{checked}")
print(f"  No good vertex: {no_good}/{checked}")


# Check at n=6 too
print(f"\nn=6: checking...")
n = 6
total = 2 ** (n * (n - 1) // 2)
checked = 0
all_good = 0
no_good = 0
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

    if not is_sc(A, n):
        continue
    b1 = compute_b1(A, n)
    if b1 != 0:
        continue

    kappa = vertex_connectivity(A, n)
    if kappa < 2:
        continue

    checked += 1
    has_good = False
    for v in range(n):
        others = [x for x in range(n) if x != v]
        B, _ = get_induced(A, n, others)
        if compute_b1(B, n - 1) == 0:
            has_good = True
            break
    if has_good:
        all_good += 1
    else:
        no_good += 1
        print(f"  COUNTEREXAMPLE: bits={bits}")

    if checked % 1000 == 0:
        print(f"  checked {checked} ({time.time()-t0:.0f}s)")

elapsed = time.time() - t0
print(f"\nn=6 ({elapsed:.0f}s): {checked} SC + b1=0 + kappa>=2 tournaments")
print(f"  Has good vertex: {all_good}/{checked}")
print(f"  No good vertex: {no_good}/{checked}")


# ============================================================
# Part 3: Breakdown by kappa for good vertex existence
# ============================================================
print(f"\n{'=' * 70}")
print("BREAKDOWN: SC + b1=0 - HOW DOES GOOD VERTEX ARISE?")
print("=" * 70)

for n in [5, 6]:
    total = 2 ** (n * (n - 1) // 2)
    categories = Counter()
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

        b1 = compute_b1(A, n)
        sc = is_sc(A, n)

        if b1 != 0:
            categories[('b1=1', 'trivially good')] += 1
            continue

        if not sc:
            # Good vertex via not-SC argument
            categories[('b1=0,not-SC', 'good via condensation')] += 1
            continue

        # SC + b1=0
        # Check: is there v with T\v not SC?
        has_notsc_deletion = False
        has_sc_b1_0_deletion = False
        for v in range(n):
            others = [x for x in range(n) if x != v]
            B, _ = get_induced(A, n, others)
            sc_v = is_sc(B, n - 1)
            if not sc_v:
                has_notsc_deletion = True
                break  # good vertex found (not SC -> b1=0)

        if has_notsc_deletion:
            categories[('SC+b1=0', 'good via not-SC deletion')] += 1
            continue

        # All deletions are SC. Check b1 of each.
        has_good = False
        for v in range(n):
            others = [x for x in range(n) if x != v]
            B, _ = get_induced(A, n, others)
            if compute_b1(B, n - 1) == 0:
                has_good = True
                break

        if has_good:
            categories[('SC+b1=0,all-del-SC', 'good via SC+b1=0 deletion')] += 1
        else:
            categories[('SC+b1=0,all-del-SC', 'NO GOOD VERTEX')] += 1

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s):")
    for cat, count in sorted(categories.items()):
        print(f"  {cat[0]}: {cat[1]}: {count}")


print("\n\nDone.")
