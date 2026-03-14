"""
alpha4_exhaustive_n6.py -- kind-pasteur-2026-03-14-S66

Exhaustive check at n=6 for alpha_1=4 alpha_2 values.
Also: examine WHY alpha_2 can only be 0 or 4 for alpha_1=4.

Key hypothesis: alpha_1=4 forcing theorem (HYP-1047)
  For any tournament with alpha_1=4:
    alpha_2 in {0, 4} only (never 1,2,3,5,6).

Mechanism analysis: When alpha_1=4 has a disjoint pair,
the cycle structure forces a K_{2,2} bipartite arrangement
(two disjoint pairs of overlapping cycles, creating 2x2=4
cross-disjoint pairs).
"""

import numpy as np
from itertools import combinations
from collections import Counter

def count_directed_hamcycles(A, vertices):
    k = len(vertices)
    if k < 3 or k % 2 == 0:
        return 0
    vlist = list(vertices)
    sub = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i != j:
                sub[i][j] = int(A[vlist[i]][vlist[j]])
    full = (1 << k) - 1
    dp = [[0]*k for _ in range(1 << k)]
    dp[1][0] = 1
    for mask in range(1, 1 << k):
        for v in range(k):
            if dp[mask][v] == 0:
                continue
            for u in range(1, k):
                if mask & (1 << u):
                    continue
                if sub[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    total = 0
    for v in range(1, k):
        if dp[full][v] and sub[v][0]:
            total += dp[full][v]
    return total

def main():
    print("=" * 70)
    print("EXHAUSTIVE alpha_1=4 ANALYSIS AT n=6")
    print("=" * 70)

    n = 6
    total_edges = n * (n - 1) // 2  # 15
    total_tournaments = 1 << total_edges  # 32768

    a1_4_count = 0
    a2_values = Counter()
    examples = {}  # alpha_2 -> list of (bits, cycles)

    for bits in range(total_tournaments):
        A = np.zeros((n, n), dtype=np.int8)
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if bits & (1 << idx):
                    A[j][i] = 1
                else:
                    A[i][j] = 1
                idx += 1

        # Find all directed odd cycles
        cycles = []
        for size in range(3, n + 1, 2):
            for subset in combinations(range(n), size):
                cnt = count_directed_hamcycles(A, list(subset))
                if cnt > 0:
                    cycles.append((frozenset(subset), cnt, size))

        alpha_1 = sum(cnt for _, cnt, _ in cycles)

        if alpha_1 == 4:
            a1_4_count += 1
            # Compute alpha_2
            alpha_2 = 0
            for i in range(len(cycles)):
                for j in range(i + 1, len(cycles)):
                    if len(cycles[i][0] & cycles[j][0]) == 0:
                        alpha_2 += cycles[i][1] * cycles[j][1]
            a2_values[alpha_2] += 1

            if alpha_2 not in examples:
                examples[alpha_2] = []
            if len(examples[alpha_2]) < 3:
                examples[alpha_2].append((bits, cycles))

        if (bits + 1) % 10000 == 0:
            print(f"  {bits + 1}/{total_tournaments}...")

    print(f"\nTotal tournaments: {total_tournaments}")
    print(f"alpha_1=4: {a1_4_count}")
    print(f"\nalpha_2 distribution for alpha_1=4:")
    for a2 in sorted(a2_values.keys()):
        pct = 100 * a2_values[a2] / a1_4_count
        print(f"  alpha_2={a2}: {a2_values[a2]} ({pct:.1f}%)")

    for a2 in sorted(examples.keys()):
        print(f"\nExamples with alpha_2={a2}:")
        for bits, cycles in examples[a2]:
            print(f"  bits={bits}")
            for vs, cnt, sz in cycles:
                print(f"    Cycle: {sorted(vs)} (size {sz}, {cnt} directed)")
            # Check pairwise disjointness
            for i in range(len(cycles)):
                for j in range(i + 1, len(cycles)):
                    inter = cycles[i][0] & cycles[j][0]
                    status = "DISJOINT" if len(inter) == 0 else f"share {sorted(inter)}"
                    print(f"    {i}-{j}: {status}")

    # Also check: alpha_1 in {2,3,5,6} forcing patterns
    print("\n" + "=" * 70)
    print("ALPHA_2 VALUES FOR EACH alpha_1 AT n=6 (EXHAUSTIVE)")
    print("=" * 70)

    a1_a2_pairs = Counter()
    a1_counts = Counter()

    for bits in range(total_tournaments):
        A = np.zeros((n, n), dtype=np.int8)
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if bits & (1 << idx):
                    A[j][i] = 1
                else:
                    A[i][j] = 1
                idx += 1

        cycles = []
        for size in range(3, n + 1, 2):
            for subset in combinations(range(n), size):
                cnt = count_directed_hamcycles(A, list(subset))
                if cnt > 0:
                    cycles.append((frozenset(subset), cnt, size))

        alpha_1 = sum(cnt for _, cnt, _ in cycles)
        a1_counts[alpha_1] += 1

        alpha_2 = 0
        for i in range(len(cycles)):
            for j in range(i + 1, len(cycles)):
                if len(cycles[i][0] & cycles[j][0]) == 0:
                    alpha_2 += cycles[i][1] * cycles[j][1]

        a1_a2_pairs[(alpha_1, alpha_2)] += 1

    print("\nalpha_1 | alpha_2 values")
    print("-" * 40)
    for a1 in sorted(a1_counts.keys()):
        a2_vals = {a2: cnt for (a, a2), cnt in a1_a2_pairs.items() if a == a1}
        a2_str = ", ".join(f"{a2}({cnt})" for a2, cnt in sorted(a2_vals.items()))
        print(f"  {a1:4d}  | {a2_str}")

    # Compute H and T for each (alpha_1, alpha_2) pair
    print("\n" + "=" * 70)
    print("H AND T VALUES AT n=6 (EXHAUSTIVE)")
    print("=" * 70)
    print("\n(alpha_1, alpha_2) -> T, H, count")
    H_values = set()
    for (a1, a2), cnt in sorted(a1_a2_pairs.items()):
        T = a1 + 2 * a2
        H = 1 + 2 * a1 + 4 * a2  # Only if alpha_3=0
        # Need to check alpha_3 for full H
        # For now, compute T
        print(f"  ({a1:3d}, {a2:3d}) -> T={T:4d}, H_min={H:4d}, count={cnt}")
        H_values.add(H)

    # List achievable H values
    print("\nAchievable H values (alpha_3=0 assumption):")
    for h in sorted(H_values):
        print(f"  H={h}")

if __name__ == "__main__":
    main()
