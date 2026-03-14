"""
allconflict_threshold.py — kind-pasteur-2026-03-14-S66

Find the exact threshold where (a1, 0) becomes tournament-achievable.
(a1, 0) means: a1 odd directed cycles, ALL pairwise sharing a vertex (clique in CG).

Known:
  a1=1: achievable (n=3, one 3-cycle)
  a1=2: achievable (n=4, two sharing a vertex)
  a1=3: IMPOSSIBLE for all n (H=7 gap, splicing lemma)

Question: what is the smallest a1 >= 4 with (a1, 0) achievable?
And: does (a1, 0) become permanently impossible for all a1 >= 3?
      Or does it become achievable again at some a1 >= 4?

Strategy:
  1. Exhaustive at small n: for each tournament, compute alpha_1, alpha_2.
     Record all (alpha_1, alpha_2) pairs with alpha_2 = 0.
  2. Sample at larger n.
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

def find_odd_cycles(A, n, max_len=None):
    """Find all directed odd cycles in tournament A."""
    if max_len is None:
        max_len = n
    cycles = []
    # Find directed cycles by DFS on each vertex subset of odd size
    for size in range(3, max_len+1, 2):
        for subset in combinations(range(n), size):
            S = list(subset)
            # Check all cyclic orderings of this subset
            # A directed cycle on S means a Hamiltonian cycle in the subtournament A[S]
            sub = A[np.ix_(S, S)]
            # Count directed Hamiltonian cycles in sub
            # For small size, enumerate permutations
            from itertools import permutations
            for perm in permutations(range(len(S))):
                if perm[0] != 0:  # fix first vertex to avoid counting rotations
                    continue
                is_cycle = True
                for k in range(len(S)):
                    if sub[perm[k]][perm[(k+1) % len(S)]] != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    cycle_verts = frozenset(S)
                    cycles.append(cycle_verts)
                    break  # one direction is enough to confirm cycle exists on this set
    return cycles

def cycles_conflict(c1, c2):
    """Two cycles conflict if they share a vertex."""
    return len(c1 & c2) > 0

def compute_alpha(cycles):
    """Compute alpha_1 (# cycles) and alpha_2 (# disjoint pairs)."""
    alpha_1 = len(cycles)
    alpha_2 = 0
    cycle_list = list(cycles)
    for i in range(len(cycle_list)):
        for j in range(i+1, len(cycle_list)):
            if not cycles_conflict(cycle_list[i], cycle_list[j]):
                alpha_2 += 1
    return alpha_1, alpha_2

def main():
    print("=" * 70)
    print("ALL-CONFLICTING THRESHOLD: min a1 with (a1, 0) achievable")
    print("=" * 70)

    # Exhaustive for small n
    for n in range(3, 8):
        num_edges = n * (n - 1) // 2
        total = 1 << num_edges

        allconflict_values = {}  # a1 -> count of tournaments
        min_a2_by_a1 = {}  # a1 -> min alpha_2

        for bits in range(total):
            A = tournament_from_bits(n, bits)
            cycles = find_odd_cycles(A, n)
            a1, a2 = compute_alpha(cycles)

            if a1 not in min_a2_by_a1 or a2 < min_a2_by_a1[a1]:
                min_a2_by_a1[a1] = a2

            if a2 == 0 and a1 > 0:
                allconflict_values[a1] = allconflict_values.get(a1, 0) + 1

        print(f"\nn={n}: {total} tournaments")
        print(f"  All-conflicting (alpha_2=0) values of alpha_1:")
        for a1 in sorted(allconflict_values.keys()):
            print(f"    a1={a1}: {allconflict_values[a1]} tournaments, H={1+2*a1}")

        print(f"  Min alpha_2 by alpha_1:")
        for a1 in sorted(min_a2_by_a1.keys()):
            print(f"    a1={a1}: min_a2={min_a2_by_a1[a1]}")

    # For n=8, sample
    print(f"\n{'='*70}")
    print(f"n=8: SAMPLING (too large for exhaustive)")
    print(f"{'='*70}")

    rng = np.random.default_rng(2026_03_14)
    n = 8
    num_samples = 20000
    allconflict_n8 = {}
    min_a2_by_a1_n8 = {}

    for trial in range(num_samples):
        # Random tournament
        A = np.zeros((n, n), dtype=np.int8)
        for i in range(n):
            for j in range(i+1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        cycles = find_odd_cycles(A, n)
        a1, a2 = compute_alpha(cycles)

        if a1 not in min_a2_by_a1_n8 or a2 < min_a2_by_a1_n8[a1]:
            min_a2_by_a1_n8[a1] = a2

        if a2 == 0 and a1 > 0:
            allconflict_n8[a1] = allconflict_n8.get(a1, 0) + 1

    print(f"  {num_samples} random tournaments sampled")
    print(f"  All-conflicting (alpha_2=0) values of alpha_1:")
    for a1 in sorted(allconflict_n8.keys()):
        print(f"    a1={a1}: {allconflict_n8[a1]} tournaments")

    print(f"  Min alpha_2 by alpha_1:")
    for a1 in sorted(min_a2_by_a1_n8.keys()):
        print(f"    a1={a1}: min_a2={min_a2_by_a1_n8[a1]}")

    # Key question: what is max a1 with a2=0?
    print(f"\n{'='*70}")
    print("SUMMARY: Maximum alpha_1 with alpha_2=0 (all-conflicting)")
    print(f"{'='*70}")

    # Also: analyze the STRUCTURE of all-conflicting tournaments
    # A tournament with all cycles pairwise conflicting means CG(T) is a clique
    # This is a STRONG constraint. How many vertices can all cycles overlap?

    print("\nSTRUCTURAL ANALYSIS:")
    print("If all odd cycles pairwise share a vertex, consider the")
    print("'conflict hypergraph' where each cycle is a hyperedge.")
    print("A clique in the conflict graph = a set of cycles, every pair sharing a vertex.")
    print("This is like a 'sunflower' (or near-sunflower) structure.")
    print()
    print("Key theorem (Helly for cycles): If d+1 convex sets in R^d")
    print("are pairwise intersecting, they have a common point.")
    print("For vertex sets of size 3 (3-cycles): Helly dimension 2")
    print("=> 3 pairwise-intersecting 3-cycles SHOULD have a common vertex")
    print("... unless they violate convexity (which discrete sets always can).")
    print()
    print("But the splicing lemma (opus-S71e) shows that even without")
    print("a common vertex, extra cycles are forced, so a1=3 with a2=0")
    print("is impossible regardless.")

if __name__ == "__main__":
    main()
