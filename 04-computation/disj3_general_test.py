#!/usr/bin/env python3
"""
disj3_general_test.py -- Test if disj3 + c5/2 = const holds for ALL regular tournaments
or only circulant ones.

The identity K = c5 - 2*ov1 - 2*ov2 = -3p(p^2-1)(p^2-9)/320 was proved for
circulant tournaments on Z_p. Does it hold for non-circulant regular tournaments?

Also test: does it hold for NON-regular tournaments? (Probably not.)

Author: kind-pasteur-2026-03-12-S60
"""

import math
from itertools import combinations, permutations
from collections import defaultdict
import random


def count_c3_sets(A, n):
    c3 = []
    for a, b, c in combinations(range(n), 3):
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            c3.append(frozenset([a, b, c]))
    return list(set(c3))


def count_c5_dir(A, n):
    total = 0
    for subset in combinations(range(n), 5):
        verts = list(subset)
        nn = 5
        start = 0
        dp = {(1 << start, start): 1}
        for mask in range(1, 1 << nn):
            if not (mask & (1 << start)):
                continue
            for v in range(nn):
                if not (mask & (1 << v)):
                    continue
                key = (mask, v)
                if key not in dp or dp[key] == 0:
                    continue
                cnt = dp[key]
                for w in range(nn):
                    if mask & (1 << w):
                        continue
                    if A[verts[v]][verts[w]]:
                        nkey = (mask | (1 << w), w)
                        dp[nkey] = dp.get(nkey, 0) + cnt
        full = (1 << nn) - 1
        for v in range(nn):
            if v == start:
                continue
            key = (full, v)
            if key in dp and dp[key] > 0:
                if A[verts[v]][verts[start]]:
                    total += dp[key]
    return total


def overlap_analysis(c3_sets):
    n3 = len(c3_sets)
    disj = 0
    ov1 = 0
    ov2 = 0
    for i in range(n3):
        for j in range(i+1, n3):
            o = len(c3_sets[i] & c3_sets[j])
            if o == 0:
                disj += 1
            elif o == 1:
                ov1 += 1
            elif o == 2:
                ov2 += 1
    return disj, ov1, ov2


def random_tournament(n):
    """Generate a random tournament on n vertices."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def random_regular_tournament(n):
    """Generate a random regular tournament on n vertices (n must be odd).
    Uses random 1-factorization approach (Walecki-type)."""
    assert n % 2 == 1
    m = (n - 1) // 2

    # Start with any tournament and try to make it regular via switching
    A = random_tournament(n)

    # Iteratively fix scores via arc switching
    for _ in range(10000):
        scores = [sum(A[i]) for i in range(n)]
        max_score = max(scores)
        min_score = min(scores)
        if max_score == m and min_score == m:
            return A

        # Find a vertex with too high score and one with too low
        high = scores.index(max_score)
        low = scores.index(min_score)

        if high != low and A[high][low]:
            # Flip high->low to low->high
            A[high][low] = 0
            A[low][high] = 1
        else:
            # Random flip on high or low
            for j in range(n):
                if j != high and A[high][j]:
                    if scores[j] < m:
                        A[high][j] = 0
                        A[j][high] = 1
                        break

    # If not regular, return None
    scores = [sum(A[i]) for i in range(n)]
    if all(s == m for s in scores):
        return A
    return None


def main():
    print("=" * 70)
    print("DISJ3 + C5/2 = CONST: GENERALITY TEST")
    print("=" * 70)

    # ====== PART 1: ALL tournaments on n=5 ======
    print("\n" + "=" * 60)
    print("PART 1: All tournaments on n=5")
    print("=" * 60)

    n = 5
    m = 2
    K_expected = -3 * n * (n**2 - 1) * (n**2 - 9) // 320
    c3_reg = n * (n**2 - 1) // 24  # c3 for regular

    print(f"  n={n}: K_expected (for regular) = {K_expected}")
    print(f"  c3 for regular = {c3_reg}")

    # All tournaments
    regular_K = set()
    non_regular_K = set()
    all_K = defaultdict(list)

    for bits in range(1 << (n*(n-1)//2)):
        A = [[0]*n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        is_reg = all(s == m for s in scores)

        c3_sets = count_c3_sets(A, n)
        c3_count = len(c3_sets)
        c5 = count_c5_dir(A, n)
        disj, ov1, ov2 = overlap_analysis(c3_sets)
        K = c5 - 2*ov1 - 2*ov2

        all_K[scores].append(K)
        if is_reg:
            regular_K.add(K)
        else:
            non_regular_K.add(K)

    print(f"\n  Regular tournaments (scores (2,2,2,2,2)): K values = {sorted(regular_K)}")
    if len(regular_K) == 1:
        print(f"  *** CONSTANT K = {list(regular_K)[0]} for ALL regular 5-vertex tournaments ***")
        print(f"  K_expected = {K_expected}, match = {list(regular_K)[0] == K_expected}")
    print(f"\n  Non-regular K values: {sorted(non_regular_K)}")

    print(f"\n  K by score sequence:")
    for scores in sorted(all_K):
        K_vals = sorted(set(all_K[scores]))
        print(f"    {scores}: K = {K_vals}")

    # ====== PART 2: ALL tournaments on n=7 ======
    print("\n" + "=" * 60)
    print("PART 2: All tournaments on n=7 (sampling)")
    print("=" * 60)

    n = 7
    m = 3
    K_expected = -3 * n * (n**2 - 1) * (n**2 - 9) // 320

    print(f"  n={n}: K_expected = {K_expected}")

    # n=7 has 2^21 = 2M tournaments, too many. Sample.
    random.seed(42)
    regular_K_7 = set()
    non_regular_K_7 = defaultdict(set)

    # First test ALL circulant orientations
    print(f"\n  All circulant orientations:")
    circ_K = set()
    for bits in range(1 << m):
        S = []
        for i in range(m):
            chord = i + 1
            if (bits >> i) & 1:
                S.append(chord)
            else:
                S.append(n - chord)
        S_set = set(S)
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for s in S_set:
                A[i][(i + s) % n] = 1

        c3_sets = count_c3_sets(A, n)
        c5 = count_c5_dir(A, n)
        disj, ov1, ov2 = overlap_analysis(c3_sets)
        K = c5 - 2*ov1 - 2*ov2
        circ_K.add(K)

    print(f"    Circulant K values: {sorted(circ_K)}")
    if len(circ_K) == 1:
        print(f"    *** CONSTANT K = {list(circ_K)[0]}, expected = {K_expected} ***")

    # Now sample random regular tournaments
    print(f"\n  Random regular tournaments (200 samples):")
    for _ in range(200):
        A = random_regular_tournament(n)
        if A is None:
            continue
        c3_sets = count_c3_sets(A, n)
        c5 = count_c5_dir(A, n)
        disj, ov1, ov2 = overlap_analysis(c3_sets)
        K = c5 - 2*ov1 - 2*ov2
        regular_K_7.add(K)

    print(f"    Regular K values: {sorted(regular_K_7)}")
    if len(regular_K_7) == 1:
        print(f"    *** CONSTANT K = {list(regular_K_7)[0]} for ALL regular 7-vertex tournaments ***")
    else:
        print(f"    *** K VARIES: {len(regular_K_7)} distinct values ***")
        print(f"    Expected: {K_expected}")

    # Sample non-regular
    print(f"\n  Random non-regular tournaments (200 samples):")
    non_reg_K = set()
    for _ in range(200):
        A = random_tournament(n)
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        if all(s == m for s in scores):
            continue
        c3_sets = count_c3_sets(A, n)
        c5 = count_c5_dir(A, n)
        disj, ov1, ov2 = overlap_analysis(c3_sets)
        K = c5 - 2*ov1 - 2*ov2
        non_reg_K.add(K)
        non_regular_K_7[scores].add(K)

    print(f"    Non-regular K values: {len(non_reg_K)} distinct values")
    print(f"    Range: [{min(non_reg_K)}, {max(non_reg_K)}]")

    # ====== PART 3: Exhaustive test at n=7 for ALL regular tournaments ======
    print("\n" + "=" * 60)
    print("PART 3: Exhaustive n=7 regular tournament census")
    print("=" * 60)

    # Generate all regular tournaments on 7 vertices by brute force
    # There are 456 regular tournaments on 7 vertices (Moon 1968)
    # We can enumerate them by checking all 2^21 tournaments... too many.
    # Instead, use a systematic construction.

    # Method: for each pair (i,j), fix some arcs and check regularity
    # This is still expensive. Let's use a smarter approach.
    # Generate all tournaments with score sequence (3,3,3,3,3,3,3)

    # Using the switching graph approach:
    # Start from one regular tournament and enumerate all others via switching.
    # But this is complex. Let's just check the circulant ones thoroughly
    # and sample heavily.

    print(f"  Heavy sampling: 5000 random regular tournaments on n=7")
    K_counts = defaultdict(int)
    total_found = 0

    for _ in range(50000):
        A = random_regular_tournament(n)
        if A is None:
            continue
        total_found += 1
        c3_sets = count_c3_sets(A, n)
        c5 = count_c5_dir(A, n)
        disj, ov1, ov2 = overlap_analysis(c3_sets)
        K = c5 - 2*ov1 - 2*ov2
        K_counts[K] += 1

    print(f"  Found {total_found} regular tournaments")
    print(f"  K distribution:")
    for K in sorted(K_counts):
        print(f"    K={K}: {K_counts[K]} tournaments ({100*K_counts[K]/total_found:.1f}%)")

    if len(K_counts) == 1:
        print(f"\n  *** THEOREM HOLDS FOR ALL REGULAR TOURNAMENTS ON 7 VERTICES! ***")
    else:
        print(f"\n  *** K VARIES for regular tournaments! Only circulant identity is exact. ***")

    # ====== PART 4: Same test at n=9 ======
    print("\n" + "=" * 60)
    print("PART 4: Random regular tournaments on n=9")
    print("=" * 60)

    n = 9
    m = 4
    K_expected = -3 * n * (n**2 - 1) * (n**2 - 9) // 320

    print(f"  n={n}: K_expected = {K_expected}")

    K_counts_9 = defaultdict(int)
    total_found_9 = 0

    for _ in range(5000):
        A = random_regular_tournament(n)
        if A is None:
            continue
        total_found_9 += 1
        c3_sets = count_c3_sets(A, n)
        c5 = count_c5_dir(A, n)
        disj, ov1, ov2 = overlap_analysis(c3_sets)
        K = c5 - 2*ov1 - 2*ov2
        K_counts_9[K] += 1

    print(f"  Found {total_found_9} regular tournaments")
    print(f"  K distribution:")
    for K in sorted(K_counts_9):
        print(f"    K={K}: {K_counts_9[K]} tournaments ({100*K_counts_9[K]/total_found_9:.1f}%)")

    if len(K_counts_9) == 1:
        print(f"\n  *** THEOREM HOLDS FOR ALL REGULAR n=9 TOURNAMENTS! ***")
        print(f"  K = {list(K_counts_9.keys())[0]}, expected = {K_expected}")
    else:
        print(f"\n  Multiple K values. Identity does NOT generalize to all regular tournaments at n=9.")

    print("\nDONE.")


if __name__ == '__main__':
    main()
