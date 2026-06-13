#!/usr/bin/env python3
"""
H=21 mod analysis: what congruence classes are achievable?

H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...

H mod 4: H = 1 + 2*alpha_1 mod 4.
  alpha_1 even -> H = 1 mod 4
  alpha_1 odd  -> H = 3 mod 4

H = 21 = 1 mod 4, so alpha_1 must be even (alpha_1 = 10).

H mod 8: H = 1 + 2*alpha_1 + 4*alpha_2 mod 8.
  21 mod 8 = 5, so 2*alpha_1 + 4*alpha_2 = 4 mod 8.
  => alpha_1 + 2*alpha_2 = 2 mod 4.
  alpha_1 = 10 -> 10 + 2*alpha_2 = 2 mod 4 -> 2*alpha_2 = -8 = 0 mod 4 -> alpha_2 = 0 mod 2.
  So alpha_2 is even.

More precisely: H = 21 needs alpha_1 = 10 and alpha_2 = even.
20 + 4*alpha_2 = 20 -> alpha_2 = 0.
Wait: H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ... = 21
=> 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ... = 20
=> alpha_1 + 2*alpha_2 + 4*alpha_3 + ... = 10.

So the exact requirement is: sum_k alpha_k * 2^{k-1} = 10.
With alpha_k = number of independent sets of size k in Omega.

Possible decompositions (alpha_1 = total cycles, alpha_2 = disjoint pairs, ...):
alpha_1=10, rest=0: 10 cycles, all pairwise conflicting. Omega = K_10.
alpha_1=8, alpha_2=1: 8 cycles, 1 disjoint pair. Omega = K_8 - edge.
alpha_1=6, alpha_2=2: 6 cycles, 2 disjoint pairs.
alpha_1=4, alpha_2=3: 4 cycles, 3 disjoint pairs.
alpha_1=2, alpha_2=4: need alpha_2=4 but only C(2,2)=1 pair -> impossible.
alpha_1=6, alpha_2=0, alpha_3=1: 6+0+4=10. But alpha_3=1 requires 3 pairwise
  disjoint cycles. This forces alpha_2 >= 3 (any 2 from the triple are disjoint).
  Contradiction with alpha_2=0.
alpha_1=4, alpha_2=1, alpha_3=1: 4+2+4=10. alpha_3=1 means 3 pairwise disjoint
  cycles among 4 total. Choose 3 from 4: all C(4,3)=4 triples must have alpha_3 count.
  alpha_2=1 means only 1 disjoint pair among 4 cycles. But if 3 are pairwise disjoint,
  alpha_2 >= 3. Contradiction.

So the only valid decompositions are with alpha_3 = 0:
(10, 0): K_10
(8, 1): K_8 - edge
(6, 2): see specific graphs
(4, 3): complement of triangle-free graph on 4 vertices

Now let's check what (alpha_1, alpha_2) values are achievable at each n.

kind-pasteur-2026-03-07-S31
"""

from itertools import combinations
from collections import defaultdict
import time


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


def count_directed_cycles_per_set(n, A):
    """Count directed 3-cycles (vertex sets) and 5-cycles (vertex sets with multiplicity)."""
    three_cycles = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if (A[a][b] and A[b][c] and A[c][a]) or \
                   (A[a][c] and A[c][b] and A[b][a]):
                    three_cycles.append(frozenset({a, b, c}))
    return three_cycles


def analyze_alpha_distribution(n):
    """At each n, find the (alpha_1, alpha_2) distribution and check H=21 candidates."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    total = 1 << m

    # For each tournament, compute alpha_1 (total 3-cycles) and alpha_2 (disjoint pairs)
    # This uses ONLY 3-cycles, which underestimates alpha_1 at n>=5 (missing 5-cycles).
    # But let's first check with 3-cycles only.

    alpha_pairs = defaultdict(int)
    h_by_alpha = defaultdict(set)

    t0 = time.time()
    for bits in range(total):
        adj_bits = [0]*n
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                A[j][i] = 1
                adj_bits[j] |= (1 << i)
            else:
                A[i][j] = 1
                adj_bits[i] |= (1 << j)

        cycles = count_directed_cycles_per_set(n, A)
        alpha1 = len(cycles)

        # Count disjoint pairs
        alpha2 = 0
        for i in range(alpha1):
            for j in range(i+1, alpha1):
                if not (cycles[i] & cycles[j]):
                    alpha2 += 1

        alpha_pairs[(alpha1, alpha2)] += 1

        # H value
        H = held_karp(n, adj_bits)
        h_by_alpha[(alpha1, alpha2)].add(H)

    elapsed = time.time() - t0
    print(f"n={n} ({elapsed:.1f}s):")

    # What (alpha_1, alpha_2) could give H=21?
    # H = 1 + 2*alpha_1_full + 4*alpha_2_full + ...
    # where alpha_1_full includes 5-cycles etc.
    # For 3-cycle-only: H_3cycle = 1 + 2*alpha1 + 4*alpha2 (partial)
    print(f"  3-cycle-only alpha pairs near target alpha1+2*alpha2=10:")
    for (a1, a2) in sorted(alpha_pairs.keys()):
        target = a1 + 2*a2
        if 8 <= target <= 12:
            H_set = h_by_alpha[(a1, a2)]
            print(f"    (alpha1={a1}, alpha2={a2}): target={target}, count={alpha_pairs[(a1,a2)]}, "
                  f"H values={sorted(H_set)}")

    # Check if H=21 appears
    all_H = set()
    for h_set in h_by_alpha.values():
        all_H |= h_set
    print(f"  H=21 achievable: {21 in all_H}")

    # Show gap structure around 21
    for h in [19, 21, 23]:
        if h in all_H:
            # Find which (a1,a2) give this H
            sources = [(a1, a2) for (a1, a2), hs in h_by_alpha.items() if h in hs]
            print(f"  H={h}: from {sources}")
        else:
            print(f"  H={h}: NOT ACHIEVABLE")


if __name__ == '__main__':
    for n in [5, 6, 7]:
        analyze_alpha_distribution(n)
        print()
