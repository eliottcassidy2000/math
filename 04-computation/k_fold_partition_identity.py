#!/usr/bin/env python3
"""
MAJOR THEOREM: k-fold partition identity for tournaments.

For any tournament on n vertices and any 2k distinct vertices,
sum over (2k)! perms of A[P_1][P_2]*A[P_3][P_4]*...*A[P_{2k-1}][P_{2k}] = (2k)!/2^k.

Proof: decompose into (2k-1)!! perfect matchings. Each matching contributes
k! (from k! slot assignments) * 1 (from product of (A[a][b]+A[b][a])=1) = k!.
Total = (2k-1)!! * k! = (2k)!/2^k. QED.

Consequence: D_S = n!/2^k for any k pairwise non-adjacent positions.

Verify for k=1,2,3,4.
"""
from itertools import permutations
from math import factorial, comb
import random

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def test_kfold(k, num_trials=100):
    """Test k-fold partition identity on 2k vertices."""
    n = 2*k
    expected = factorial(n) // (2**k)
    failures = 0

    for _ in range(num_trials):
        A = random_tournament(n)
        total = 0
        for P in permutations(range(n)):
            prod = 1
            for i in range(k):
                prod *= A[P[2*i]][P[2*i+1]]
            total += prod

        if total != expected:
            failures += 1
            if failures <= 3:
                print(f"  k={k}: FAIL! got {total}, expected {expected}")

    return failures

print("=== k-fold Partition Identity ===")
print("sum_{P ∈ S_{2k}} prod_{i=1}^k A[P_{2i-1}][P_{2i}] = (2k)!/2^k")
print()

for k in range(1, 6):
    n = 2*k
    expected = factorial(n) // (2**k)
    trials = 200 if k <= 3 else 50
    fails = test_kfold(k, trials)
    print(f"  k={k} (2k={n}): expected={(n)}!/{2**k}={expected}, "
          f"verified {trials} tournaments, failures={fails}")

# ============ Consequence for D_S ============
print("\n\n=== Consequence: D_S = n!/2^k for non-adjacent positions ===")

def compute_DS(A, n, positions):
    """D_S = #{P: A[P_i][P_{i+1}]=1 for all i in positions}."""
    count = 0
    for P in permutations(range(n)):
        ok = all(A[P[pos]][P[pos+1]] for pos in positions)
        if ok:
            count += 1
    return count

# Test at n=7 with various k
for n in [5, 7]:
    print(f"\n  n={n}:")
    edge_positions = list(range(n-1))

    # Find non-adjacent position sets for various k
    from itertools import combinations as combos
    for k in range(1, (n+1)//2):
        non_adj_sets = []
        for S in combos(edge_positions, k):
            if all(abs(S[i]-S[j]) >= 2 for i in range(k) for j in range(i+1, k)):
                non_adj_sets.append(S)

        if not non_adj_sets:
            print(f"    k={k}: no non-adjacent sets of this size")
            continue

        expected = factorial(n) // (2**k)
        # Test with a random tournament
        A = random_tournament(n)

        results = set()
        for S in non_adj_sets:
            DS = compute_DS(A, n, S)
            results.add(DS)

        all_match = all(r == expected for r in results)
        print(f"    k={k}: {len(non_adj_sets)} non-adj sets, "
              f"expected D_S=n!/2^k={expected}, "
              f"results={results}, match={'YES' if all_match else 'NO'}")

# ============ PROOF sketch ============
print("\n\n=== PROOF of k-fold Partition Identity ===")
print("""
THEOREM: For any tournament T on m ≥ 2k vertices and any 2k distinct
vertices v_1,...,v_{2k}, we have:

  Σ_{σ ∈ S_{2k}} Π_{i=1}^k A[v_{σ(2i-1)}][v_{σ(2i)}] = (2k)! / 2^k

PROOF:
Consider the sum S = Σ_{σ ∈ S_{2k}} Π_{i=1}^k A[v_{σ(2i-1)}][v_{σ(2i)}].

Group the 2k positions into k slots: {1,2}, {3,4}, ..., {2k-1,2k}.

For any perfect matching M = {{a_1,b_1},...,{a_k,b_k}} of {v_1,...,v_{2k}}:
- There are k! ways to assign pairs to slots.
- For each assignment, there are 2^k orderings within slots.
- Summing over orderings within slot i:
    A[a_i][b_i] + A[b_i][a_i] = 1 (tournament axiom).
- Product over slots = 1^k = 1.
- So each (matching, slot-assignment) pair contributes 1.
- Each matching contributes k! × 1 = k!.

Number of perfect matchings = (2k)! / (k! × 2^k) = (2k-1)!!.

Total: S = (2k-1)!! × k! = (2k)! / 2^k.   QED.

COROLLARY: For k pairwise non-adjacent edge positions {s_1,...,s_k}
in a permutation of n vertices:

  D_S = n! / 2^k.

PROOF: The 2k vertices at positions {s_i, s_i+1} are all distinct (non-adjacency).
Summing over all n! permutations, condition on which 2k vertices appear at key
positions. For any choice: (n-2k)! arrangements of remaining vertices, and
(2k)!/2^k from the k-fold identity. Total: C(n,2k) × (2k)!/2^k × (n-2k)! = n!/2^k.
                                                                              QED.

This is the MASTER IDENTITY underlying the universal congruence theorem.
The tournament-dependent part of D_k comes ONLY from position sets
with at least one adjacent pair (shared vertex). These carry extra factors
of 2 from the constrained counting structure.
""")
