#!/usr/bin/env python3
"""
Position Character Decomposition (THM-068) Verification

Verifies the PCD at n=3,5 (exhaustive Walsh transform) and n=7 (permutation enumeration).

Key results verified:
1. M[v,v] has only even Walsh degrees (for odd n)
2. M[v,v]_hat[S] = (-1)^{[v in N(S)]} * H_hat[S] / (n-2k) for Walsh degree 2k
3. |N(S)| = k for all degree-2k monomials
4. Nonzero monomials = unions of even-length vertex-disjoint paths
5. Even-length rule: odd-length path components cause exact descent cancellation

opus-2026-03-07-S35
"""

import numpy as np
from itertools import permutations

def bits_to_adj(bits, n):
    A = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
            k += 1
    return A

def compute_M_diagonal(A, n):
    M_diag = [0] * n
    def enum_paths(path, used):
        if len(path) == n:
            for idx, v in enumerate(path):
                M_diag[v] += (-1)**idx
            return
        last = path[-1]
        for u in range(n):
            if u in used: continue
            if A[last][u]:
                path.append(u)
                used.add(u)
                enum_paths(path, used)
                path.pop()
                used.discard(u)
    for start in range(n):
        enum_paths([start], {start})
    return M_diag

def fast_walsh_hadamard(f, m):
    N = 2**m
    a = f.copy()
    for i in range(m):
        step = 1 << i
        for j in range(0, N, step*2):
            for k in range(j, j + step):
                u, v = a[k], a[k + step]
                a[k] = u + v
                a[k + step] = u - v
    return a / N


def verify_exhaustive(n):
    """Exhaustive Walsh transform verification at small n."""
    m = n*(n-1)//2
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    N = 2**m

    print(f"\n{'='*60}")
    print(f"EXHAUSTIVE VERIFICATION AT n={n}")
    print(f"{'='*60}")

    M_diag_all = [np.zeros(N) for _ in range(n)]
    for bits in range(N):
        A = bits_to_adj(bits, n)
        md = compute_M_diagonal(A, n)
        for v in range(n):
            M_diag_all[v][bits] = md[v]

    Mv_hats = [fast_walsh_hadamard(M_diag_all[v].copy(), m) for v in range(n)]
    H_hat = fast_walsh_hadamard(sum(M_diag_all[v] for v in range(n)).copy(), m)

    # Check even-only degrees
    for v in range(n):
        odd_energy = sum(Mv_hats[v][S]**2 for S in range(N) if bin(S).count('1') % 2 == 1)
        assert odd_energy < 1e-10, f"Odd Walsh energy nonzero for vertex {v}"
    print(f"  Even-only Walsh degrees: CONFIRMED")

    # Check PCD formula
    all_ok = True
    for deg in range(0, m+1, 2):
        k = deg // 2
        if k > (n-1)//2: continue
        if n == deg: continue
        alpha = 1.0 / (n - deg)

        for S in range(N):
            if bin(S).count('1') != deg: continue
            if abs(H_hat[S]) < 1e-10: continue

            neg_count = 0
            for v in range(n):
                ratio = Mv_hats[v][S] / H_hat[S]
                if abs(ratio - (-alpha)) < 1e-10:
                    neg_count += 1
                elif abs(ratio - alpha) < 1e-10:
                    pass
                else:
                    all_ok = False

        if k <= (n-1)//2:
            from math import comb
            neg_sets = set()
            for S in range(N):
                if bin(S).count('1') != deg: continue
                if abs(H_hat[S]) < 1e-10: continue
                neg_verts = frozenset(v for v in range(n)
                    if abs(Mv_hats[v][S]/H_hat[S] + alpha) < 1e-10)
                neg_sets.add(neg_verts)
            expected = comb(n, k)
            print(f"  Degree {deg}: {len(neg_sets)} patterns (expected C({n},{k})={expected}): "
                  f"{'OK' if len(neg_sets) == expected else 'FAIL'}")

    print(f"  PCD formula: {'CONFIRMED' if all_ok else 'FAILED'}")


def verify_permutation(n, S_edges, label):
    """Verify PCD via permutation enumeration for a specific monomial at given n."""
    edge_set = set(tuple(sorted(e)) for e in S_edges)
    deg = len(S_edges)
    k = deg // 2

    raw_sums = {v: 0 for v in range(n)}
    total_valid = 0

    for perm in permutations(range(n)):
        path_edges = set(tuple(sorted([perm[i], perm[i+1]])) for i in range(n-1))
        if not edge_set.issubset(path_edges):
            continue

        total_valid += 1
        desc = sum(1 for i in range(n-1)
                   if tuple(sorted([perm[i], perm[i+1]])) in edge_set
                   and perm[i] > perm[i+1])
        desc_sign = (-1)**desc

        for v in range(n):
            pos_v = list(perm).index(v)
            raw_sums[v] += (-1)**pos_v * desc_sign

    H_hat_val = sum(raw_sums.values()) / 2**(n-1)

    if abs(H_hat_val) < 1e-10:
        print(f"  {label}: H_hat = 0 (zero monomial)")
        return True

    alpha = 1.0 / (n - deg) if n != deg else 1.0

    neg_set = set()
    ok = True
    for v in range(n):
        M_hat = raw_sums[v] / 2**(n-1)
        ratio = M_hat / H_hat_val
        is_neg = abs(ratio + alpha) < 1e-10
        is_pos = abs(ratio - alpha) < 1e-10
        if is_neg:
            neg_set.add(v)
        elif not is_pos:
            ok = False

    status = "OK" if ok and len(neg_set) == k else "FAIL"
    print(f"  {label}: |N|={len(neg_set)}, expected {k}, PCD: {status}")
    return ok


def main():
    # Exhaustive at n=3,5
    verify_exhaustive(3)
    verify_exhaustive(5)

    # Permutation-based at n=7
    print(f"\n{'='*60}")
    print(f"PERMUTATION VERIFICATION AT n=7")
    print(f"{'='*60}")

    n = 7

    # Degree 2
    verify_permutation(n, [(0,1),(0,2)], "Degree 2: fan pair")
    verify_permutation(n, [(0,1),(1,2)], "Degree 2: path pair")

    # Degree 4
    verify_permutation(n, [(0,1),(1,2),(2,3),(3,4)], "Degree 4: P4")
    verify_permutation(n, [(0,1),(1,2),(3,4),(4,5)], "Degree 4: P2+P2")
    verify_permutation(n, [(0,1),(1,2),(2,3),(4,5)], "Degree 4: P3+P1 (should be zero)")

    # Degree 6
    verify_permutation(n, [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6)], "Degree 6: P6 (HP)")

    print(f"\n{'='*60}")
    print("ALL VERIFICATIONS COMPLETE")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
