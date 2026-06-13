"""
p4_topology_check.py -- kind-pasteur-2026-03-14-S66

Check if the P_4 topology for (4,3) is achievable at n=8,9.

P_4 = path graph on 4 vertices: disjoint pairs (1,2), (2,3), (3,4).
This requires 4 directed odd cycles where:
  - Cycles 1,2 vertex-disjoint
  - Cycles 2,3 vertex-disjoint
  - Cycles 3,4 vertex-disjoint
  - Cycles 1,3 share a vertex
  - Cycles 1,4 share a vertex
  - Cycles 2,4 share a vertex

With alpha_1 = EXACTLY 4 (no other cycles in the tournament).

At n=8 with all 3-cycles: needs exactly 4 3-cycle vertex sets and
no 5-cycles or 7-cycles. Very constrained.

Strategy: Exhaustive search at n=8 for alpha_1=4, alpha_2=3.
If not found, try sampling at n=9.
"""

import numpy as np
from itertools import combinations

def tournament_from_bits(n, bits):
    A = np.zeros((n, n), dtype=np.int8)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[j][i] = 1
            else:
                A[i][j] = 1
            idx += 1
    return A

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
    print("P_4 TOPOLOGY CHECK FOR (4,3)")
    print("=" * 70)

    # At n=8: exhaustive search for alpha_1=4, alpha_2=3
    # n=8 has 28 edges, 2^28 = 268M tournaments — too many for exhaustive.
    # Instead: sample with score-filtered search.

    # A tournament with c3=4 at n=8 needs sum C(s_i,2) = 56-4 = 52.
    # Score sum = 28. Need scores with large spread.
    # Let's find score sequences with sum C(s_i,2) = 52.

    print("\n--- Score sequences with c3=4 at n=8 ---")
    print("Need sum C(s_i,2) = 52, sum s_i = 28")

    valid_scores = []
    def gen_scores(pos, remaining, lo, acc, n=8):
        if pos == n:
            if remaining == 0:
                s = sum(si*(si-1)//2 for si in acc)
                c3 = n*(n-1)*(n-2)//6 - s
                if c3 == 4:
                    valid_scores.append(tuple(acc))
            return
        hi = min(remaining, n-1)
        for si in range(lo, hi+1):
            if remaining - si >= (n-1-pos) * si:
                gen_scores(pos+1, remaining-si, si, acc + [si], n)

    gen_scores(0, 28, 0, [])
    print(f"Score sequences with c3=4: {len(valid_scores)}")
    for s in valid_scores[:10]:
        print(f"  {s}: sum C(s_i,2)={sum(si*(si-1)//2 for si in s)}")

    # Targeted sampling: generate tournaments with specific score sequences
    print("\n--- Targeted search at n=8 ---")
    n = 8
    rng = np.random.default_rng(2026_03_14_43)

    found_43 = 0
    found_a4 = 0
    total_checked = 0

    for trial in range(500000):
        A = np.zeros((n, n), dtype=np.int8)
        for i in range(n):
            for j in range(i+1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        # Quick c3 count
        c3 = 0
        for a, b, c in combinations(range(n), 3):
            if (A[a][b] and A[b][c] and A[c][a]) or \
               (A[a][c] and A[c][b] and A[b][a]):
                c3 += 1

        if c3 > 6:  # Need c3 <= 4 for alpha_1=4, allow some slack for 5-cycles
            continue

        total_checked += 1

        # Compute full alpha_1 and alpha_2
        cycles = []
        for size in range(3, n+1, 2):
            for subset in combinations(range(n), size):
                cnt = count_directed_hamcycles(A, list(subset))
                if cnt > 0:
                    cycles.append((frozenset(subset), cnt, size))

        alpha_1 = sum(cnt for _, cnt, _ in cycles)
        if alpha_1 == 4:
            found_a4 += 1
            # Compute alpha_2
            alpha_2 = 0
            for i in range(len(cycles)):
                for j in range(i+1, len(cycles)):
                    if len(cycles[i][0] & cycles[j][0]) == 0:
                        alpha_2 += cycles[i][1] * cycles[j][1]

            T = alpha_1 + 2*alpha_2
            H = 1 + 2*alpha_1 + 4*alpha_2
            if alpha_2 == 3:
                found_43 += 1
                print(f"  FOUND (4,3)! H={H}, T={T}")
                # Check topology
                for ci, cj in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]:
                    if ci < len(cycles) and cj < len(cycles):
                        disj = len(cycles[ci][0] & cycles[cj][0]) == 0
                        print(f"    cycles {ci}-{cj}: {'DISJOINT' if disj else 'CONFLICT'}")
            else:
                if total_checked <= 20 or (total_checked % 100 == 0 and found_a4 <= 10):
                    pass  # Don't spam

        if (trial+1) % 100000 == 0:
            print(f"  {trial+1}/500000, checked={total_checked}, a1=4: {found_a4}, (4,3): {found_43}")

    print(f"\nResults at n=8:")
    print(f"  Total sampled: 500000")
    print(f"  Score-filtered (c3<=6): {total_checked}")
    print(f"  alpha_1=4: {found_a4}")
    print(f"  (4,3) found: {found_43}")

    if found_a4 > 0 and found_43 == 0:
        print(f"\n  alpha_1=4 EXISTS at n=8 but (4,3) does NOT.")
        print(f"  What alpha_2 values does alpha_1=4 have?")

    # Also check: at n=8, what is the minimum alpha_2 for alpha_1=4?
    print("\n--- alpha_1=4 analysis at n=8 ---")
    rng2 = np.random.default_rng(2026_03_14_44)
    a4_a2_values = []

    for trial in range(500000):
        A = np.zeros((n, n), dtype=np.int8)
        for i in range(n):
            for j in range(i+1, n):
                if rng2.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        c3 = 0
        for a, b, c in combinations(range(n), 3):
            if (A[a][b] and A[b][c] and A[c][a]) or \
               (A[a][c] and A[c][b] and A[b][a]):
                c3 += 1
        if c3 > 6:
            continue

        cycles = []
        for size in range(3, n+1, 2):
            for subset in combinations(range(n), size):
                cnt = count_directed_hamcycles(A, list(subset))
                if cnt > 0:
                    cycles.append((frozenset(subset), cnt, size))

        alpha_1 = sum(cnt for _, cnt, _ in cycles)
        if alpha_1 == 4:
            alpha_2 = 0
            for i in range(len(cycles)):
                for j in range(i+1, len(cycles)):
                    if len(cycles[i][0] & cycles[j][0]) == 0:
                        alpha_2 += cycles[i][1] * cycles[j][1]
            a4_a2_values.append(alpha_2)

        if (trial+1) % 100000 == 0:
            if a4_a2_values:
                print(f"  {trial+1}/500000, a1=4 found: {len(a4_a2_values)}, "
                      f"a2 range: [{min(a4_a2_values)}, {max(a4_a2_values)}]")

    if a4_a2_values:
        from collections import Counter
        a2_counts = Counter(a4_a2_values)
        print(f"\nalpha_1=4 at n=8: {len(a4_a2_values)} found")
        print("alpha_2 distribution:")
        for a2 in sorted(a2_counts.keys()):
            print(f"  alpha_2={a2}: {a2_counts[a2]}")
    else:
        print("\nNo alpha_1=4 found at n=8 in sample!")

if __name__ == "__main__":
    main()
