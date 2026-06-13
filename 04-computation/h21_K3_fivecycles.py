#!/usr/bin/env python3
"""
Check: tournaments with exactly 3 pairwise-intersecting 3-cycles at n=6.
Do they have 5-cycles? If so, the full Omega K3 is NOT isolated.

kind-pasteur-2026-03-07-S31
"""

from itertools import combinations, permutations
from collections import defaultdict


def held_karp(n, adj_bits):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            c = dp[S][v]
            if c == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if adj_bits[v] & (1 << u):
                    dp[S | (1 << u)][u] += c
    return sum(dp[(1 << n) - 1])


def check_exactly_3_cycles_n6():
    n = 6
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)

    count = 0
    has_5cycle = 0
    no_5cycle = 0
    H_vals = defaultdict(int)

    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                A[j][i] = 1
            else:
                A[i][j] = 1

        cycles = []
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    if (A[a][b] and A[b][c] and A[c][a]) or \
                       (A[a][c] and A[c][b] and A[b][a]):
                        cycles.append(frozenset({a, b, c}))

        if len(cycles) != 3:
            continue

        # Check pairwise intersecting
        if not all(cycles[i] & cycles[j] for i in range(3) for j in range(i+1, 3)):
            continue  # not K3 (but we showed all triples at n=6 with 3 cycles ARE K3)

        count += 1

        # Check for 5-cycles
        found_5 = False
        num_5 = 0
        for verts5 in combinations(range(n), 5):
            v0 = verts5[0]
            rest = verts5[1:]
            for p in permutations(rest):
                seq = (v0,) + p
                ok = all(A[seq[q]][seq[(q+1)%5]] for q in range(5))
                if ok:
                    found_5 = True
                    num_5 += 1
                    break  # just need existence per 5-set

        if found_5:
            has_5cycle += 1
        else:
            no_5cycle += 1

        # Compute H
        adj_bits = [0]*n
        for i in range(n):
            for j in range(n):
                if i != j and A[i][j]:
                    adj_bits[i] |= (1 << j)
        H = held_karp(n, adj_bits)
        H_vals[H] += 1

    print(f"Tournaments with exactly 3 pairwise-intersecting 3-cycles: {count}")
    print(f"  With 5-cycles: {has_5cycle}")
    print(f"  Without 5-cycles: {no_5cycle}")
    print(f"  H values: {dict(sorted(H_vals.items()))}")

    if no_5cycle > 0:
        print("\nWARNING: Some have NO 5-cycles!")
        print("These would have I(Omega,2) = I(K3, 2) = 7")
        print("But we know H=7 is not achievable... contradiction?")
    else:
        print("\nALL have 5-cycles, so full Omega is larger than K3")
        print("I(full Omega, 2) != I(K3, 2) = 7")


if __name__ == '__main__':
    check_exactly_3_cycles_n6()
