"""
three_power_bound.py -- kind-pasteur-2026-03-14-S66

THEOREM (3^k Lower Bound): If tournament T has k mutually vertex-disjoint
odd directed cycles, then H(T) >= 3^k.

PROOF: The k disjoint cycles form an independent set of size k in Omega(T).
Any j-element subset of these k cycles is an independent set of size j.
So alpha_j >= C(k,j) for all 0 <= j <= k.
Therefore: H(T) = I(Omega(T), 2) = sum alpha_j * 2^j
                >= sum_{j=0}^k C(k,j) * 2^j = (1+2)^k = 3^k.

GENERAL: For any graph G with independence number >= k and any x > 0:
I(G, x) >= (1+x)^k.

At x=2 (the tournament evaluation point): I(G, 2) >= 3^k.

This script verifies the bound and explores its tightness.
"""

import numpy as np
from itertools import combinations
from math import comb

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=np.int8)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_ham_paths(A, n):
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

def count_directed_hamcycles_sub(A, vertices):
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

def find_max_disjoint_3cycles(A, n):
    """Find maximum number of mutually vertex-disjoint 3-cycles."""
    # Find all 3-cycle vertex sets
    cycle_sets = []
    for a, b, c in combinations(range(n), 3):
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            cycle_sets.append(frozenset([a, b, c]))

    if not cycle_sets:
        return 0, []

    # Greedy + backtracking for max independent matching
    # For small counts, try all subsets
    best = 0
    best_collection = []

    def backtrack(idx, current, used_vertices):
        nonlocal best, best_collection
        if len(current) > best:
            best = len(current)
            best_collection = list(current)
        if idx >= len(cycle_sets):
            return
        # Pruning: remaining cycles + current can't beat best
        if len(current) + len(cycle_sets) - idx <= best:
            return
        for i in range(idx, len(cycle_sets)):
            cs = cycle_sets[i]
            if not (cs & used_vertices):
                current.append(cs)
                backtrack(i + 1, current, used_vertices | cs)
                current.pop()

    if len(cycle_sets) <= 50:
        backtrack(0, [], frozenset())
    else:
        # Greedy for large cases
        used = set()
        for cs in sorted(cycle_sets, key=lambda s: len(s)):
            if not (cs & used):
                best_collection.append(cs)
                used |= cs
        best = len(best_collection)

    return best, best_collection

def main():
    print("=" * 70)
    print("THE 3^k LOWER BOUND THEOREM")
    print("=" * 70)

    print("\nTHEOREM: If tournament T has k mutually vertex-disjoint odd")
    print("directed cycles, then H(T) >= 3^k.")
    print()
    print("PROOF: k disjoint cycles => alpha_j >= C(k,j) for 0<=j<=k")
    print("  H = I(Omega,2) >= sum C(k,j)*2^j = (1+2)^k = 3^k.")
    print()
    print("COROLLARY: The minimum H for maximum disjoint cycle count k:")
    for k in range(0, 8):
        print(f"  k={k}: H >= 3^{k} = {3**k}")

    print()
    print("APPLICATION TO PERMANENT GAPS:")
    print(f"  H=7:  T=(7-1)/2=3.  Since 3 < (3^2-1)/2=4, alpha_2=0.")
    print(f"        Only decomposition (3,0). Blocked by splicing lemma.")
    print(f"  H=21: T=(21-1)/2=10. Since 10 < (3^3-1)/2=13, alpha_3=0.")
    print(f"        Only quadratic decompositions. All 4 blocked.")
    print(f"  H=63: T=31 > 13. alpha_3>=1 possible. (21,3,1) works at n=9.")
    print(f"        NOT a permanent gap (found at n=8).")
    print()

    # Part 2: Verify the bound computationally
    print("=" * 70)
    print("VERIFICATION: H >= 3^k at small n")
    print("=" * 70)

    rng = np.random.default_rng(2026_03_14_66)

    for n in range(3, 10):
        num_samples = min(5000, 1 << (n*(n-1)//2))
        is_exhaustive = (1 << (n*(n-1)//2)) <= num_samples
        if is_exhaustive:
            num_samples = 1 << (n*(n-1)//2)

        violations = 0
        tight_cases = {}  # k -> count where H == 3^k exactly
        min_ratio = {}  # k -> min H/3^k

        for trial in range(num_samples):
            if is_exhaustive:
                bits = trial
                A = np.zeros((n, n), dtype=np.int8)
                idx = 0
                for i in range(n):
                    for j in range(i+1, n):
                        if bits & (1 << idx):
                            A[j][i] = 1
                        else:
                            A[i][j] = 1
                        idx += 1
            else:
                A = random_tournament(n, rng)

            H = count_ham_paths(A, n)
            k, _ = find_max_disjoint_3cycles(A, n)

            if k > 0:
                bound = 3 ** k
                if H < bound:
                    violations += 1
                if H == bound:
                    tight_cases[k] = tight_cases.get(k, 0) + 1
                if k not in min_ratio or H / bound < min_ratio[k]:
                    min_ratio[k] = H / bound

        mode = "EXHAUSTIVE" if is_exhaustive else f"SAMPLED ({num_samples})"
        print(f"\nn={n} ({mode}): violations = {violations}")
        if tight_cases:
            for k in sorted(tight_cases.keys()):
                print(f"  k={k}: H=3^{k}={3**k} achieved {tight_cases[k]} times"
                      f" (min ratio = {min_ratio.get(k, '?'):.4f})")
        if min_ratio:
            for k in sorted(min_ratio.keys()):
                if k not in tight_cases:
                    print(f"  k={k}: min H/3^{k} = {min_ratio[k]:.4f}"
                          f" (min H = {int(min_ratio[k] * 3**k)})")

    # Part 3: Connection to 2-3 number theory
    print(f"\n{'='*70}")
    print("CONNECTION TO 2-3 NUMBER THEORY")
    print(f"{'='*70}")
    print()
    print("The 3^k bound comes from (1+x)^k at x=2.")
    print("If we evaluated at x=1 instead: I(G,1) >= (1+1)^k = 2^k.")
    print("At x=3: I(G,3) >= 4^k.")
    print()
    print("The tournament evaluation x=2 gives EXACTLY the base-3 hierarchy!")
    print("  x=2: bound is 3^k = (2+1)^k")
    print("  This is WHY 3 is the key: 3 = evaluation_point + 1 = 2+1")
    print()
    print("The hierarchy of achievable H values:")
    print("  Level 0: H in {1} (alpha(Omega) = 0, no cycles)")
    print("  Level 1: H in [3, 8] odd (alpha(Omega) >= 1, at most 1 cycle island)")
    print("    Gap: H=7 (blocked at quadratic level)")
    print("  Level 2: H in [9, 26] odd (alpha(Omega) >= 2, can have disjoint pair)")
    print("    Gap: H=21 (blocked at cubic level)")
    print("  Level 3: H in [27, 80] odd (alpha(Omega) >= 3)")
    print("    No gaps! cubic I.P. provides enough flexibility")
    print("  Level 4+: no gaps.")
    print()
    print("The gap at level 1 is H = 2^2 + 2 + 1 = 7 = Phi_3(2)")
    print("The gap at level 2 is H = 3 * 7 = 21 = 3 * Phi_3(2)")
    print()
    print("Why 3*7 = 21 specifically?")
    print("  T = 10 is the largest value below (3^3-1)/2 = 13 that")
    print("  has NO valid quadratic decomposition.")
    print("  T=11 works via (11,0) or (9,1)")
    print("  T=12 works via (12,0)")
    print("  T=10: (10,0) blocked, (8,1) blocked, (6,2) blocked, (4,3) blocked")
    print()
    print("Connection to cyclotomic polynomials:")
    print(f"  Phi_1(2) = {1}")
    print(f"  Phi_2(2) = {3}")
    print(f"  Phi_3(2) = {7}")
    print(f"  Phi_4(2) = {5}")
    print(f"  Phi_6(2) = {3}")
    print(f"  H=7 = Phi_3(2)")
    print(f"  H=21 = Phi_2(2) * Phi_3(2) = 3 * 7")
    print(f"  H=7 is the value of x^2+x+1 at x=2")
    print(f"  H=21 is the value of (x+1)(x^2+x+1) = x^3+2x^2+2x+1 at x=2 = 8+8+4+1")
    print(f"  Note: x^3+2x^2+2x+1 = (x^3-1)/(x-1) + x^2 + x")
    print(f"  Actually: 21 = 2^4 + 2^2 + 2^0 = 16+4+1 (binary: 10101)")

    # Part 4: The overlap weight and why 10 is special
    print(f"\n{'='*70}")
    print("WHY T=10 IS THE UNIQUE LEVEL-2 GAP")
    print(f"{'='*70}")
    print()
    print("At level 2 (9 <= H < 27, i.e., 4 <= T <= 12):")
    print("  T=4: achievable as (4,0) at n=5 or (2,1) at n=6")
    print("  T=5: achievable as (5,0) at n=5")
    print("  T=6: achievable as (6,0)")
    print("  T=7: achievable as (7,0)")
    print("  T=8: achievable as (8,0) or (6,1)")
    print("  T=9: achievable as (9,0)")
    print("  T=10: BLOCKED -- all 4 decompositions impossible")
    print("  T=11: achievable as (11,0) or (9,1)")
    print("  T=12: achievable as (12,0)")
    print()
    print("The pattern: (a1, 0) is achievable for all a1 != 3, 10.")
    print("And (a1, a2) with a2 > 0 starts at T=4.")
    print("So T=10 is blocked because:")
    print("  1. (10, 0) requires alpha_1=10, alpha_2=0 -- but alpha_1=10 forces alpha_2>=2")
    print("  2. (8, 1) never occurs as an achievable (alpha_1, alpha_2) pair")
    print("  3. (6, 2) never occurs as an achievable pair")
    print("  4. (4, 3) never occurs as an achievable pair")
    print()
    print("The achievable alpha_1 values with alpha_2=0:")
    print("  n=5: {0,1,2,4,5,6,7} -- note 3 MISSING!")
    print("  n=6: {0,1,2,4,5,6,7,8,9,11,12,16} -- note 3,10 MISSING!")
    print("  n=7: {0,1,2,4,5,6,7,8,9,11,12,13,15,16,18,19,22,23,24,25,34,...}")
    print("        -- note 3 has min_a2=2, 10 has min_a2=2")
    print()
    print("CONJECTURE: alpha_1=3 and alpha_1=10 are the ONLY values")
    print("that PERMANENTLY force alpha_2 >= 2 at all n.")

if __name__ == "__main__":
    main()
