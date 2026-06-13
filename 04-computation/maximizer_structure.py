"""
maximizer_structure.py — Structural properties that force the H-maximizer to be circulant

At n=7, the 240 global H-maximizers are ALL relabelings of the Paley tournament.
What structural property distinguishes them?

Key question: Is the H-maximizer always REGULAR (all vertices have the same out-degree)?
If so, and if it always has maximal automorphism group, the Paley conjecture follows.

Approach: analyze the n=7 data and the near-maximizers too.

Author: opus-2026-03-12-S60
"""
import sys
import time
import numpy as np
from collections import defaultdict, Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


def hamiltonian_paths_dp(A, n):
    dp = defaultdict(lambda: defaultdict(int))
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        if not dp[mask]:
            continue
        for v in dp[mask]:
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


def score_sequence(A, n):
    return tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))


def is_regular(A, n):
    degs = [sum(A[i][j] for j in range(n)) for i in range(n)]
    return len(set(degs)) == 1


def aut_group_size(A, n):
    """Count automorphisms by brute force (only for small n)."""
    from itertools import permutations
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(n):
                if A[i][j] != A[perm[i]][perm[j]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            count += 1
    return count


def count_3cycles(A, n):
    count = 0
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                for k in range(n):
                    if A[j][k] and A[k][i]:
                        count += 1
    return count // 3


def main():
    print("STRUCTURAL PROPERTIES OF H-MAXIMIZERS")
    print("=" * 75)

    # n=7: exhaustive
    n = 7
    print(f"\nn={n}: Exhaustive enumeration")
    t0 = time.time()

    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)

    # Collect top H values and their properties
    H_to_props = defaultdict(list)

    total = 0
    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                A[i][j] = 1
            else:
                A[j][i] = 1
        H = hamiltonian_paths_dp(A, n)
        if H >= 171:  # Top 3 values
            ss = score_sequence(A, n)
            reg = is_regular(A, n)
            c3 = count_3cycles(A, n)
            H_to_props[H].append({
                'score_seq': ss, 'regular': reg, 'C3': c3
            })
        total += 1
        if total % 500000 == 0:
            elapsed = time.time() - t0
            print(f"  ... {total}/{2**m} ({elapsed:.1f}s)")

    elapsed = time.time() - t0
    print(f"  {total} tournaments, {elapsed:.1f}s")

    for H in sorted(H_to_props.keys(), reverse=True):
        props = H_to_props[H]
        ss_dist = Counter(p['score_seq'] for p in props)
        reg_count = sum(1 for p in props if p['regular'])
        c3_dist = Counter(p['C3'] for p in props)
        print(f"\n  H = {H} ({len(props)} tournaments)")
        print(f"    Score sequences: {dict(ss_dist)}")
        print(f"    Regular: {reg_count}/{len(props)}")
        print(f"    C3 distribution: {dict(c3_dist)}")

    # Compute |Aut| for a few maximizers
    print(f"\n  Automorphism groups of H={max(H_to_props.keys())} maximizers (sampling 5):")
    max_H = max(H_to_props.keys())
    # Reconstruct a few maximizer adjacency matrices
    count = 0
    for bits in range(1 << m):
        if count >= 5:
            break
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                A[i][j] = 1
            else:
                A[j][i] = 1
        H = hamiltonian_paths_dp(A, n)
        if H == max_H:
            aut = aut_group_size(A, n)
            ss = score_sequence(A, n)
            print(f"    |Aut| = {aut}, score = {ss}")
            count += 1

    # Key question: is EVERY regular tournament with maximal C3 a Paley relabeling?
    print(f"\n  REGULAR tournaments with C3=14:")
    count_reg = 0
    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                A[i][j] = 1
            else:
                A[j][i] = 1
        if is_regular(A, n):
            c3 = count_3cycles(A, n)
            H = hamiltonian_paths_dp(A, n)
            if count_reg < 20:
                print(f"    C3={c3}, H={H}")
            count_reg += 1

    print(f"  Total regular: {count_reg}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
