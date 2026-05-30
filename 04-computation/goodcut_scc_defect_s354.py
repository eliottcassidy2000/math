#!/usr/bin/env python3
"""Verify THM-354: good-cut count equals n minus SCC count.

The theorem itself is structural. This script keeps a small reproducible audit:

1. Exhaust all labeled tournaments through n=6 and all Hamiltonian paths.
2. Exhaust fixed-base tilings through n=7.
3. Sample labeled n=7,8 tournaments and fixed-base n=8 tilings.
"""

from itertools import permutations
from random import Random


def tournament_from_bits(n, bits):
    A = [[0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            k += 1
    return A


def fixed_base_tiling(n, bits):
    """Tournament with base path 0->1->...->n-1.

    Bits are assigned to nonconsecutive pairs. A 1 means the arc points backward
    relative to the base path, so it contributes an upward crossing interval.
    """
    A = [[0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if j == i + 1:
                A[i][j] = 1
            else:
                if (bits >> k) & 1:
                    A[j][i] = 1
                else:
                    A[i][j] = 1
                k += 1
    return A


def hamiltonian_paths(A):
    n = len(A)
    for path in permutations(range(n)):
        if all(A[path[i]][path[i + 1]] for i in range(n - 1)):
            yield path


def good_cut_count_for_path(A, path):
    n = len(path)
    total = 0
    for k in range(1, n):
        good = False
        for late in range(k, n):
            for early in range(k):
                if A[path[late]][path[early]]:
                    good = True
                    break
            if good:
                break
        if good:
            total += 1
    return total


def scc_count(A):
    n = len(A)
    reach = [[bool(A[i][j]) or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                row_i = reach[i]
                row_k = reach[k]
                for j in range(n):
                    row_i[j] = row_i[j] or row_k[j]

    seen = [False] * n
    count = 0
    for i in range(n):
        if seen[i]:
            continue
        count += 1
        for j in range(n):
            if reach[i][j] and reach[j][i]:
                seen[j] = True
    return count


def exact_labeled(max_n=6):
    print("EXACT LABELED TOURNAMENTS")
    for n in range(3, max_n + 1):
        total = 1 << (n * (n - 1) // 2)
        violations = 0
        mixed_paths = 0
        bucket_counts = {}
        for bits in range(total):
            A = tournament_from_bits(n, bits)
            expected = n - scc_count(A)
            vals = {good_cut_count_for_path(A, path) for path in hamiltonian_paths(A)}
            if vals != {expected}:
                violations += 1
                if len(vals) > 1:
                    mixed_paths += 1
        print(f"n={n}: checked={total}, violations={violations}, mixed_path_values={mixed_paths}")
        assert violations == 0


def exact_fixed_base(max_n=7):
    print("\nEXACT FIXED-BASE TILINGS")
    for n in range(3, max_n + 1):
        m = (n - 1) * (n - 2) // 2
        total = 1 << m
        violations = 0
        buckets = {}
        base_path = tuple(range(n))
        for bits in range(total):
            A = fixed_base_tiling(n, bits)
            g = good_cut_count_for_path(A, base_path)
            expected = n - scc_count(A)
            buckets[g] = buckets.get(g, 0) + 1
            if g != expected:
                violations += 1
        print(f"n={n}: checked={total}, violations={violations}, buckets={dict(sorted(buckets.items()))}")
        assert violations == 0


def sampled(seed=20260530):
    rng = Random(seed)
    print("\nSAMPLED CHECKS")

    for n, samples in [(7, 2000), (8, 500)]:
        total = 1 << (n * (n - 1) // 2)
        violations = 0
        for _ in range(samples):
            A = tournament_from_bits(n, rng.randrange(total))
            expected = n - scc_count(A)
            vals = set()
            for path in hamiltonian_paths(A):
                vals.add(good_cut_count_for_path(A, path))
                if len(vals) > 1 or expected not in vals:
                    break
            if vals != {expected}:
                violations += 1
        print(f"labeled n={n}: samples={samples}, violations={violations}")
        assert violations == 0

    n = 8
    m = (n - 1) * (n - 2) // 2
    base_path = tuple(range(n))
    violations = 0
    samples = 5000
    for _ in range(samples):
        A = fixed_base_tiling(n, rng.randrange(1 << m))
        if good_cut_count_for_path(A, base_path) != n - scc_count(A):
            violations += 1
    print(f"fixed-base n={n}: samples={samples}, violations={violations}")
    assert violations == 0


if __name__ == "__main__":
    print("THM-354 good-cut SCC defect verification")
    exact_labeled()
    exact_fixed_base()
    sampled()
    print("\nPASS")
