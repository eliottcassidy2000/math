"""
allconflict_threshold_fast.py -- kind-pasteur-2026-03-14-S66

Fast computation of the all-conflicting threshold:
What is the minimum alpha_1 where (alpha_1, 0) is tournament-achievable?

Uses bitwise tournament encoding and fast cycle detection.
Exhaustive at n=3..6, sampled at n=7,8.
"""

import numpy as np
from itertools import combinations

def tournament_from_bits(n, bits):
    """Build adjacency matrix from bit string."""
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

def count_directed_3cycles(A, n):
    """Count directed 3-cycles: vertex sets and count."""
    cycle_sets = []
    for a, b, c in combinations(range(n), 3):
        # Check a->b->c->a
        if A[a][b] and A[b][c] and A[c][a]:
            cycle_sets.append(frozenset([a, b, c]))
        # Check a->c->b->a
        elif A[a][c] and A[c][b] and A[b][a]:
            cycle_sets.append(frozenset([a, b, c]))
    return cycle_sets

def has_hamcycle_subtournament(A, vertices):
    """Check if subtournament on vertices has a directed Hamiltonian cycle.
    Uses DP: dp[mask][v] = True if there's a path visiting vertices in mask ending at v.
    Then check if last vertex connects back to first.
    """
    k = len(vertices)
    if k < 3:
        return False
    if k % 2 == 0:
        return False  # even cycles not counted

    vlist = list(vertices)
    # Build sub-adjacency
    sub = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i != j:
                sub[i][j] = A[vlist[i]][vlist[j]]

    # DP over all subsets
    full = (1 << k) - 1
    # dp[mask][v] = can we have a path visiting exactly the vertices in mask, ending at v?
    dp = [[False]*k for _ in range(1 << k)]
    for v in range(k):
        dp[1 << v][v] = True

    for mask in range(1, 1 << k):
        for v in range(k):
            if not dp[mask][v]:
                continue
            for u in range(k):
                if mask & (1 << u):
                    continue
                if sub[v][u]:
                    dp[mask | (1 << u)][u] = True

    # Check for Hamiltonian cycle: path visiting all vertices ending at v, with edge v->start
    # Fix start = 0 (vertex 0 in the subset) to avoid counting rotations
    for v in range(k):
        if dp[full][v] and sub[v][0]:
            return True
    return False

def count_directed_hamcycles(A, vertices):
    """Count directed Hamiltonian cycles in subtournament.
    Returns count of distinct directed cycles (not just existence).
    """
    k = len(vertices)
    if k < 3 or k % 2 == 0:
        return 0

    vlist = list(vertices)
    sub = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i != j:
                sub[i][j] = A[vlist[i]][vlist[j]]

    full = (1 << k) - 1
    dp = [[0]*k for _ in range(1 << k)]
    dp[1][0] = 1  # Fix start at vertex 0

    for mask in range(1, 1 << k):
        for v in range(k):
            if dp[mask][v] == 0:
                continue
            for u in range(1, k):  # don't return to 0 via middle
                if mask & (1 << u):
                    continue
                if sub[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]

    total = 0
    for v in range(1, k):
        if dp[full][v] and sub[v][0]:
            total += dp[full][v]
    return total

def find_all_odd_cycles(A, n):
    """Find all directed odd cycles. Returns list of (vertex_set, num_directed_cycles)."""
    cycles = []
    for size in range(3, n+1, 2):
        for subset in combinations(range(n), size):
            vs = frozenset(subset)
            cnt = count_directed_hamcycles(A, list(subset))
            if cnt > 0:
                cycles.append((vs, cnt))
    return cycles

def compute_alpha12(cycles):
    """Compute alpha_1, alpha_2 from cycle list.
    alpha_1 = total number of directed odd cycles
    alpha_2 = number of pairs of vertex-disjoint directed odd cycles

    For alpha_2: if cycle sets S1 and S2 are disjoint, they contribute
    cnt1 * cnt2 pairs of disjoint directed cycles.
    """
    alpha_1 = sum(cnt for _, cnt in cycles)
    alpha_2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            vs_i, cnt_i = cycles[i]
            vs_j, cnt_j = cycles[j]
            if len(vs_i & vs_j) == 0:
                alpha_2 += cnt_i * cnt_j
    return alpha_1, alpha_2

def main():
    print("=" * 70)
    print("ALL-CONFLICTING THRESHOLD ANALYSIS")
    print("Finding min alpha_1 where (alpha_1, 0) is achievable")
    print("=" * 70)

    # Exhaustive for n=3..6
    for n in range(3, 7):
        num_edges = n * (n - 1) // 2
        total = 1 << num_edges

        a1_with_a2_zero = {}  # a1 -> count
        a1_a2_pairs = {}  # (a1, a2) -> count
        max_a1_allconflict = 0

        for bits in range(total):
            A = tournament_from_bits(n, bits)
            cycles = find_all_odd_cycles(A, n)
            a1, a2 = compute_alpha12(cycles)

            pair = (a1, a2)
            a1_a2_pairs[pair] = a1_a2_pairs.get(pair, 0) + 1

            if a2 == 0:
                a1_with_a2_zero[a1] = a1_with_a2_zero.get(a1, 0) + 1
                if a1 > max_a1_allconflict:
                    max_a1_allconflict = a1

        print(f"\nn={n}: {total} tournaments")
        print(f"  All achievable (alpha_1, alpha_2) pairs:")
        for (a1, a2) in sorted(a1_a2_pairs.keys()):
            H = 1 + 2*a1 + 4*a2
            print(f"    (a1={a1:3d}, a2={a2:3d}): {a1_a2_pairs[(a1,a2)]:6d} tournaments, H={H}")

        print(f"\n  All-conflicting (alpha_2=0):")
        for a1 in sorted(a1_with_a2_zero.keys()):
            print(f"    a1={a1}: {a1_with_a2_zero[a1]} tournaments, H={1+2*a1}")
        print(f"  Max alpha_1 with alpha_2=0: {max_a1_allconflict}")

    # n=7 sampling
    print(f"\n{'='*70}")
    print("n=7: SAMPLING")
    print(f"{'='*70}")

    rng = np.random.default_rng(2026_03_14_66)
    n = 7
    num_samples = 5000
    a1_with_a2_zero_n7 = {}
    min_a2_by_a1_n7 = {}

    for trial in range(num_samples):
        A = random_tournament(n, rng)
        cycles = find_all_odd_cycles(A, n)
        a1, a2 = compute_alpha12(cycles)

        if a1 not in min_a2_by_a1_n7 or a2 < min_a2_by_a1_n7[a1]:
            min_a2_by_a1_n7[a1] = a2

        if a2 == 0:
            a1_with_a2_zero_n7[a1] = a1_with_a2_zero_n7.get(a1, 0) + 1

        if (trial + 1) % 1000 == 0:
            print(f"  {trial+1}/{num_samples}")

    print(f"\n  All-conflicting (alpha_2=0) at n=7:")
    for a1 in sorted(a1_with_a2_zero_n7.keys()):
        print(f"    a1={a1}: {a1_with_a2_zero_n7[a1]} tournaments, H={1+2*a1}")

    print(f"\n  Min alpha_2 by alpha_1 at n=7:")
    for a1 in sorted(min_a2_by_a1_n7.keys()):
        print(f"    a1={a1:3d}: min_a2={min_a2_by_a1_n7[a1]}")

    # n=8 sampling (much larger cycle space)
    print(f"\n{'='*70}")
    print("n=8: SAMPLING")
    print(f"{'='*70}")

    rng8 = np.random.default_rng(2026_03_14_68)
    n = 8
    num_samples = 2000
    a1_with_a2_zero_n8 = {}
    min_a2_by_a1_n8 = {}

    for trial in range(num_samples):
        A = random_tournament(n, rng8)
        cycles = find_all_odd_cycles(A, n)
        a1, a2 = compute_alpha12(cycles)

        if a1 not in min_a2_by_a1_n8 or a2 < min_a2_by_a1_n8[a1]:
            min_a2_by_a1_n8[a1] = a2

        if a2 == 0:
            a1_with_a2_zero_n8[a1] = a1_with_a2_zero_n8.get(a1, 0) + 1

        if (trial + 1) % 500 == 0:
            print(f"  {trial+1}/{num_samples}")

    print(f"\n  All-conflicting (alpha_2=0) at n=8:")
    if a1_with_a2_zero_n8:
        for a1 in sorted(a1_with_a2_zero_n8.keys()):
            print(f"    a1={a1}: {a1_with_a2_zero_n8[a1]} tournaments, H={1+2*a1}")
    else:
        print("    NONE FOUND")

    print(f"\n  Min alpha_2 by alpha_1 at n=8:")
    for a1 in sorted(min_a2_by_a1_n8.keys()):
        print(f"    a1={a1:3d}: min_a2={min_a2_by_a1_n8[a1]}")

    # Analysis
    print(f"\n{'='*70}")
    print("ANALYSIS: ALL-CONFLICTING THRESHOLD")
    print(f"{'='*70}")
    print()
    print("At n=3: max 1 cycle (one 3-cycle), alpha_2=0 trivial")
    print("At n=4: max 2 3-cycles on 4 vertices, always share 2+ vertices")
    print("At n=5: 5-cycles enter. 10 3-cycles + 1-3 5-cycles for regular")
    print("At n=6: More cycles, disjoint pairs become easier")
    print()
    print("KEY QUESTION: Is alpha_2=0 achievable for alpha_1 >= 3?")
    print("If NOT: then the all-conflicting property is bounded by alpha_1 <= 2")
    print("for ALL n, meaning H=7 is the FIRST impossible (a1,0) value")
    print("and all (a1,0) with a1 >= 3 are impossible.")

if __name__ == "__main__":
    main()
