"""
beta3_les_analysis.py — LES approach to proving beta_3 <= 1

Using the long exact sequence of (T, T\\v):
  H_3(T\\v) -> H_3(T) -> H_3(T,T\\v) -> H_2(T\\v) = 0

Since H_2(T\\v) = 0 (proved), the map H_3(T) -> H_3(T,T\\v) is surjective.
This gives: beta_3(T) <= beta_3(T\\v) + dim H_3(T,T\\v).

Proof strategy for beta_3 <= 1 by induction on n:
1. Base case: n <= 5, beta_3 = 0 always (exhaustive). DONE.
2. Induction step: Assume beta_3 <= 1 for (n-1)-vertex tournaments.
   If there EXISTS v with beta_3(T\\v) = 0, then:
   beta_3(T) <= 0 + dim H_3(T,T\\v).
   So we need: dim H_3(T,T\\v) <= 1 AND such v exists.

Questions to answer:
Q1: For beta_3(T) = 1 at n=6, is beta_3(T\\v) = 0 for all v? (Should be yes since n=5 always has beta_3=0)
Q2: For beta_3(T) = 1 at n=7, does there exist v with beta_3(T\\v) = 0?
Q3: What is dim H_3(T,T\\v)?

Also check: "beta_3 fragility" — does beta_3 > 0 vanish under some vertex deletion?

Author: kind-pasteur-S46 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
from collections import Counter, defaultdict
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def enumerate_allowed_paths(A, n, p):
    if p < 0: return []
    if p == 0: return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1: adj[i].append(j)
    paths = []
    stack = []
    for start in range(n):
        stack.append(([start], 1 << start))
        while stack:
            path, visited = stack.pop()
            if len(path) == p + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))
    return paths

def boundary_coeffs(path):
    return [((-1)**i, path[:i] + path[i+1:]) for i in range(len(path))]

def compute_omega_basis(ap, p):
    if not ap.get(p, []):
        return np.zeros((0, 0)), 0
    if p == 0:
        return np.eye(len(ap[p])), len(ap[p])

    apm1_set = set(ap.get(p-1, []))
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(ap[p]):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in apm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
    if na_count == 0:
        return np.eye(len(ap[p])), len(ap[p])

    P = np.zeros((na_count, len(ap[p])))
    for j, path in enumerate(ap[p]):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed:
                P[non_allowed[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = int(sum(s > 1e-10 for s in S))
    ns = Vt[rank:].T
    if ns.shape[1] > 0:
        return ns, ns.shape[1]
    else:
        return np.zeros((len(ap[p]), 0)), 0


def compute_beta3(A, n):
    """Compute beta_3 for tournament A."""
    ap = {}
    for p in range(min(6, n)):
        ap[p] = enumerate_allowed_paths(A, n, p)

    omega_bases = {}
    omega_dims = {}
    for p in range(min(6, n)):
        if not ap.get(p, []):
            omega_bases[p] = np.zeros((0, 0))
            omega_dims[p] = 0
        else:
            omega_bases[p], omega_dims[p] = compute_omega_basis(ap, p)

    dim3 = omega_dims.get(3, 0)
    if dim3 == 0:
        return 0

    # Build d_3
    bd3 = np.zeros((len(ap.get(2, [])), len(ap[3])))
    idx2 = {path: i for i, path in enumerate(ap.get(2, []))}
    for j, path in enumerate(ap[3]):
        for sign, face in boundary_coeffs(path):
            if face in idx2:
                bd3[idx2[face], j] += sign

    O3 = omega_bases[3]
    d3_om = bd3 @ O3
    sv3 = np.linalg.svd(d3_om, compute_uv=False)
    rank_d3 = int(sum(s > 1e-8 for s in sv3))
    ker_d3 = dim3 - rank_d3

    if ker_d3 == 0:
        return 0

    # Build d_4
    dim4 = omega_dims.get(4, 0)
    if dim4 == 0:
        return ker_d3

    bd4 = np.zeros((len(ap[3]), len(ap.get(4, []))))
    idx3 = {path: i for i, path in enumerate(ap[3])}
    for j, path in enumerate(ap.get(4, [])):
        for sign, face in boundary_coeffs(path):
            if face in idx3:
                bd4[idx3[face], j] += sign

    O4 = omega_bases[4]
    d4_om = bd4 @ O4

    O3pinv = np.linalg.pinv(O3)
    d4_omega3 = O3pinv @ d4_om
    sv4 = np.linalg.svd(d4_omega3, compute_uv=False)
    rank_d4 = int(sum(s > 1e-8 for s in sv4))

    return ker_d3 - rank_d4


def deletion_adj(A, n, v):
    """Return adjacency matrix of T\\v as (n-1)x(n-1)."""
    vertices = [i for i in range(n) if i != v]
    n1 = len(vertices)
    A1 = np.zeros((n1, n1), dtype=int)
    for i, vi in enumerate(vertices):
        for j, vj in enumerate(vertices):
            A1[i][j] = A[vi][vj]
    return A1, n1


def main():
    print("=" * 70)
    print("BETA_3 LES ANALYSIS — Towards proving beta_3 <= 1")
    print("=" * 70)

    # Part 1: n=6, beta_3=1 — check beta_3(T\\v) for each v
    print("\n--- Part 1: n=6 beta_3 deletion fragility ---")
    n = 6
    b3_fragile_count = 0
    b3_survive_count = 0
    total_b3 = 0

    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        b3 = compute_beta3(A, n)
        if b3 == 0:
            continue
        total_b3 += 1

        # Check each deletion
        max_b3_del = 0
        min_b3_del = n
        for v in range(n):
            A1, n1 = deletion_adj(A, n, v)
            b3v = compute_beta3(A1, n1)
            max_b3_del = max(max_b3_del, b3v)
            min_b3_del = min(min_b3_del, b3v)

        if max_b3_del == 0:
            b3_fragile_count += 1
        else:
            b3_survive_count += 1

    print(f"  beta_3>0 tournaments: {total_b3}")
    print(f"  All deletions have beta_3=0 (fragile): {b3_fragile_count}")
    print(f"  Some deletion has beta_3>0 (survives): {b3_survive_count}")

    # Part 2: n=7 sampled — same test
    print("\n--- Part 2: n=7 beta_3 deletion fragility ---")
    n = 7
    rng = np.random.RandomState(42)
    N = 300
    b3_count = 0
    fragile = 0
    survives = 0
    del_b3_vals = Counter()

    for trial in range(N):
        A = random_tournament(n, rng)
        b3 = compute_beta3(A, n)
        if b3 == 0:
            continue
        b3_count += 1

        del_b3_list = []
        for v in range(n):
            A1, n1 = deletion_adj(A, n, v)
            b3v = compute_beta3(A1, n1)
            del_b3_list.append(b3v)
            del_b3_vals[b3v] += 1

        max_del = max(del_b3_list)
        min_del = min(del_b3_list)
        num_zero = del_b3_list.count(0)

        if max_del == 0:
            fragile += 1
        else:
            survives += 1
            scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
            print(f"  trial={trial}: beta_3={b3}, scores={scores}, del_b3={del_b3_list}")

        if (trial + 1) % 100 == 0:
            print(f"  n=7: {trial+1}/{N} done, {b3_count} with beta_3>0", flush=True)

    print(f"\n  beta_3>0: {b3_count}")
    print(f"  Fragile (all del beta_3=0): {fragile}")
    print(f"  Survives (some del beta_3>0): {survives}")
    print(f"  Deletion beta_3 values: {dict(sorted(del_b3_vals.items()))}")

    # Part 3: When beta_3 survives deletion, what's the profile?
    print("\n--- Part 3: Detailed deletion analysis at n=7 ---")
    n = 7
    rng2 = np.random.RandomState(123)
    N2 = 500
    b3_count2 = 0
    num_zero_del = Counter()
    del_profile = Counter()

    for trial in range(N2):
        A = random_tournament(n, rng2)
        b3 = compute_beta3(A, n)
        if b3 == 0:
            continue
        b3_count2 += 1

        del_list = []
        for v in range(n):
            A1, n1 = deletion_adj(A, n, v)
            b3v = compute_beta3(A1, n1)
            del_list.append(b3v)

        num_zero = del_list.count(0)
        num_zero_del[num_zero] += 1
        prof = tuple(sorted(del_list))
        del_profile[prof] += 1

        if (trial + 1) % 100 == 0:
            print(f"  n=7: {trial+1}/{N2} done, {b3_count2} with beta_3>0", flush=True)

    print(f"\n  beta_3>0: {b3_count2}/{N2}")
    print(f"  #zeros among n deletions:")
    for k in sorted(num_zero_del.keys()):
        print(f"    {k} zeros: {num_zero_del[k]}")
    print(f"  Deletion profiles (sorted):")
    for prof, cnt in del_profile.most_common():
        print(f"    {prof}: {cnt}")

    # Part 4: n=8 small sample — does beta_3 survive?
    print("\n--- Part 4: n=8 sampled ---")
    n = 8
    rng3 = np.random.RandomState(77)
    N3 = 100
    b3_count3 = 0
    fragile3 = 0
    survives3 = 0

    for trial in range(N3):
        A = random_tournament(n, rng3)
        try:
            b3 = compute_beta3(A, n)
        except:
            continue
        if b3 == 0:
            continue
        b3_count3 += 1

        max_del = 0
        has_zero = False
        for v in range(n):
            A1, n1 = deletion_adj(A, n, v)
            b3v = compute_beta3(A1, n1)
            max_del = max(max_del, b3v)
            if b3v == 0:
                has_zero = True

        if max_del == 0:
            fragile3 += 1
        else:
            survives3 += 1

        if (trial + 1) % 25 == 0:
            print(f"  n=8: {trial+1}/{N3} done, {b3_count3} with beta_3>0", flush=True)

    print(f"\n  beta_3>0: {b3_count3}")
    print(f"  Fragile: {fragile3}")
    print(f"  Survives: {survives3}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
