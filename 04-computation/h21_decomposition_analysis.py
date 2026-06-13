"""
h21_decomposition_analysis.py -- kind-pasteur-2026-03-14-S66

H=21 requires alpha_1 + 2*alpha_2 = 10.
Four possible decompositions: (10,0), (8,1), (6,2), (4,3).
All others ruled out structurally.

For each decomposition, check if ANY tournament at any n achieves it.
Focus: which alpha_1 values can pair with which alpha_2 values?

Exhaustive at n=6, sample at n=7,8.
"""

import numpy as np
from itertools import combinations, permutations

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

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=np.int8)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_directed_hamcycles(A, vertices):
    """Count directed Hamiltonian cycles in subtournament on vertices."""
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

def compute_alpha12(A, n):
    """Compute alpha_1 and alpha_2 for tournament A on n vertices."""
    # Find all odd directed cycles
    cycles = []  # (vertex_set, count_of_directed_cycles)
    for size in range(3, n+1, 2):
        for subset in combinations(range(n), size):
            cnt = count_directed_hamcycles(A, list(subset))
            if cnt > 0:
                cycles.append((frozenset(subset), cnt))

    alpha_1 = sum(cnt for _, cnt in cycles)
    alpha_2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i][0] & cycles[j][0]) == 0:
                alpha_2 += cycles[i][1] * cycles[j][1]
    return alpha_1, alpha_2

def count_ham_paths(A, n):
    """Count Hamiltonian paths using DP."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def main():
    print("=" * 70)
    print("H=21 DECOMPOSITION ANALYSIS")
    print("H = 1 + 2*a1 + 4*a2 = 21 => a1 + 2*a2 = 10")
    print("Decompositions: (10,0), (8,1), (6,2), (4,3)")
    print("=" * 70)

    # Part 1: Exhaustive at n=6 - find all achievable (a1, a2) pairs
    print("\n--- Part 1: EXHAUSTIVE n=6 ---")
    n = 6
    num_edges = n*(n-1)//2
    total = 1 << num_edges

    pair_counts = {}
    h21_count = 0

    for bits in range(total):
        A = tournament_from_bits(n, bits)
        a1, a2 = compute_alpha12(A, n)
        H = 1 + 2*a1 + 4*a2
        pair = (a1, a2)
        pair_counts[pair] = pair_counts.get(pair, 0) + 1
        if H == 21:
            h21_count += 1

    print(f"Total tournaments: {total}")
    print(f"Tournaments with H=21: {h21_count}")
    print(f"\nAll achievable (a1, a2) pairs with a1+2*a2 near 10:")
    for (a1, a2) in sorted(pair_counts.keys()):
        t = a1 + 2*a2
        if 7 <= t <= 13:
            H = 1 + 2*a1 + 4*a2
            print(f"  ({a1:2d},{a2:2d}): T={t:3d}, H={H:3d}, count={pair_counts[(a1,a2)]}")

    # Check which (a1,a2) pairs near T=10 exist
    print("\nDecompositions of T=10 at n=6:")
    for (a1, a2) in [(10,0), (8,1), (6,2), (4,3), (2,4), (0,5)]:
        if (a1, a2) in pair_counts:
            print(f"  ({a1},{a2}): {pair_counts[(a1,a2)]} tournaments -- EXISTS")
        else:
            print(f"  ({a1},{a2}): 0 tournaments -- ABSENT")

    # Part 2: Sample at n=7
    print(f"\n--- Part 2: n=7 TARGETED SEARCH ---")
    n = 7
    rng = np.random.default_rng(2026_03_14_21)
    num_samples = 50000

    pair_counts_7 = {}
    h21_count_7 = 0

    for trial in range(num_samples):
        A = random_tournament(n, rng)
        H = count_ham_paths(A, n)
        if H == 21:
            h21_count_7 += 1
            # Compute decomposition
            a1, a2 = compute_alpha12(A, n)
            print(f"  FOUND H=21 at trial {trial}: (a1={a1}, a2={a2})")
        # For efficiency, only compute alpha decomposition when H is near 21
        if 15 <= H <= 27:
            a1, a2 = compute_alpha12(A, n)
            pair = (a1, a2)
            pair_counts_7[pair] = pair_counts_7.get(pair, 0) + 1

        if (trial + 1) % 10000 == 0:
            print(f"  {trial+1}/{num_samples}, H=21 found: {h21_count_7}")

    print(f"\nH=21 at n=7: {h21_count_7} / {num_samples}")
    print(f"\nAchievable (a1, a2) pairs near T=10 at n=7:")
    for (a1, a2) in sorted(pair_counts_7.keys()):
        t = a1 + 2*a2
        if 7 <= t <= 13:
            H = 1 + 2*a1 + 4*a2
            print(f"  ({a1:2d},{a2:2d}): T={t:3d}, H={H:3d}, count={pair_counts_7[(a1,a2)]}")

    print("\nDecompositions of T=10 at n=7:")
    for (a1, a2) in [(10,0), (8,1), (6,2), (4,3)]:
        if (a1, a2) in pair_counts_7:
            print(f"  ({a1},{a2}): {pair_counts_7[(a1,a2)]} tournaments -- EXISTS")
        else:
            print(f"  ({a1},{a2}): 0 tournaments -- ABSENT")

    # Part 3: Sample at n=8 targeting H near 21
    print(f"\n--- Part 3: n=8 TARGETED SEARCH ---")
    n = 8
    rng8 = np.random.default_rng(2026_03_14_28)
    num_samples = 20000

    pair_counts_8 = {}
    h21_count_8 = 0

    for trial in range(num_samples):
        A = random_tournament(n, rng8)
        H = count_ham_paths(A, n)
        if H == 21:
            h21_count_8 += 1
            a1, a2 = compute_alpha12(A, n)
            print(f"  FOUND H=21 at trial {trial}: (a1={a1}, a2={a2})")

        if 15 <= H <= 27:
            a1, a2 = compute_alpha12(A, n)
            pair = (a1, a2)
            pair_counts_8[pair] = pair_counts_8.get(pair, 0) + 1

        if (trial + 1) % 5000 == 0:
            print(f"  {trial+1}/{num_samples}, H=21 found: {h21_count_8}")

    print(f"\nH=21 at n=8: {h21_count_8} / {num_samples}")
    print(f"\nAchievable (a1, a2) near T=10 at n=8:")
    for (a1, a2) in sorted(pair_counts_8.keys()):
        t = a1 + 2*a2
        if 7 <= t <= 13:
            H = 1 + 2*a1 + 4*a2
            print(f"  ({a1:2d},{a2:2d}): T={t:3d}, H={H:3d}, count={pair_counts_8[(a1,a2)]}")

    print("\nDecompositions of T=10 at n=8:")
    for (a1, a2) in [(10,0), (8,1), (6,2), (4,3)]:
        if (a1, a2) in pair_counts_8:
            print(f"  ({a1},{a2}): {pair_counts_8[(a1,a2)]} tournaments -- EXISTS")
        else:
            print(f"  ({a1},{a2}): 0 tournaments -- ABSENT")

    # Part 4: Analysis
    print(f"\n{'='*70}")
    print("ANALYSIS: WHY EACH DECOMPOSITION FAILS")
    print(f"{'='*70}")

    print("""
    (10,0): alpha_1=10, alpha_2=0
      10 odd cycles, ALL pairwise sharing a vertex.
      BLOCKED: Exhaustive data shows alpha_1=10 always has alpha_2>=2.
      Mechanism: too many cycles force disjoint pairs.

    (8,1): alpha_1=8, alpha_2=1
      8 cycles, exactly 1 disjoint pair.
      Need: 8 cycles where exactly one pair is vertex-disjoint.
      This is a very specific structural requirement.

    (6,2): alpha_1=6, alpha_2=2
      6 cycles, exactly 2 disjoint pairs.
      H = 1 + 12 + 8 = 21.

    (4,3): alpha_1=4, alpha_2=3
      4 cycles, 3 disjoint pairs.
      But C(4,2)=6 pairs total. 3 disjoint means half are disjoint.
      4 cycles with 3 disjoint pairs: this constrains the cycle structure heavily.
      Actually need: among 4 cycles, exactly 3 pairs are vertex-disjoint.
    """)

if __name__ == "__main__":
    main()
