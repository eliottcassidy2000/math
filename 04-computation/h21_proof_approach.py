"""
h21_proof_approach.py -- kind-pasteur-2026-03-14-S67

Goal: Understand WHY H=21 is impossible at n=7 (and presumably all n).
H=21 requires alpha_1 + 2*alpha_2 = 10.

At n=7, alpha_1 can range widely. For each achievable alpha_1 value,
what alpha_2 values are achievable?

If no pair (a1, a2) with a1 + 2*a2 = 10 exists, we need to understand
the structural reason for each blocked case.

Approach:
1. Sample n=7 tournaments extensively, recording (a1, a2) pairs
2. For each a1 value that COULD give T=10, check what a2 values are achievable
3. Identify the gap in each case
"""

import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
import time

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

def compute_alpha_decomposed(A, n):
    """Compute alpha_1, alpha_2 and decompose by cycle type."""
    cycles_3 = []
    cycles_5 = []
    cycles_7 = []

    for size in [3, 5, 7]:
        if size > n:
            break
        for subset in combinations(range(n), size):
            cnt = count_directed_hamcycles(A, list(subset))
            if cnt > 0:
                entry = (frozenset(subset), cnt, size)
                if size == 3:
                    cycles_3.append(entry)
                elif size == 5:
                    cycles_5.append(entry)
                else:
                    cycles_7.append(entry)

    all_cycles = cycles_3 + cycles_5 + cycles_7

    # alpha_1 by type
    a1_3 = sum(cnt for _, cnt, _ in cycles_3)
    a1_5 = sum(cnt for _, cnt, _ in cycles_5)
    a1_7 = sum(cnt for _, cnt, _ in cycles_7)
    alpha_1 = a1_3 + a1_5 + a1_7

    # alpha_2 by type pair
    alpha_2 = 0
    a2_33 = 0  # pairs of disjoint 3-cycles
    a2_35 = 0  # 3-cycle disjoint from 5-cycle (impossible at n=7? 3+5=8>7)
    a2_other = 0

    for i in range(len(all_cycles)):
        for j in range(i+1, len(all_cycles)):
            if len(all_cycles[i][0] & all_cycles[j][0]) == 0:
                prod = all_cycles[i][1] * all_cycles[j][1]
                alpha_2 += prod
                si, sj = all_cycles[i][2], all_cycles[j][2]
                if si + sj == 6:  # 3+3
                    a2_33 += prod
                elif si + sj == 8:  # 3+5
                    a2_35 += prod
                else:
                    a2_other += prod

    return alpha_1, alpha_2, a1_3, a1_5, a1_7, a2_33, a2_35, a2_other

def main():
    print("=" * 70)
    print("H=21 PROOF APPROACH — n=7 (a1, a2) GAP STRUCTURE")
    print("=" * 70)

    n = 7
    rng = np.random.default_rng(2026_03_14_21)

    # Track (alpha_1, alpha_2) pairs and decomposition
    pair_counter = Counter()
    a1_to_a2 = defaultdict(set)  # alpha_1 -> set of achievable alpha_2
    a1_to_decomp = defaultdict(list)

    num_samples = 500000
    start = time.time()

    for trial in range(num_samples):
        A = np.zeros((n, n), dtype=np.int8)
        for i in range(n):
            for j in range(i+1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        a1, a2, a1_3, a1_5, a1_7, a2_33, a2_35, a2_other = compute_alpha_decomposed(A, n)
        pair_counter[(a1, a2)] += 1
        a1_to_a2[a1].add(a2)

        if len(a1_to_decomp[a1]) < 3:
            a1_to_decomp[a1].append((a2, a1_3, a1_5, a1_7, a2_33, a2_35, a2_other))

        if (trial+1) % 100000 == 0:
            elapsed = time.time() - start
            print(f"  {trial+1}/{num_samples} ({(trial+1)/elapsed:.0f}/sec)...", flush=True)

    elapsed = time.time() - start
    print(f"\nCompleted {num_samples} samples in {elapsed:.1f}s")

    # Focus on alpha_1 values where T=10 could be achieved
    print(f"\n{'='*70}")
    print("(alpha_1, alpha_2) NEAR T=10 DECOMPOSITION")
    print("=" * 70)
    print(f"\n  H=21 requires alpha_1 + 2*alpha_2 = 10")
    print(f"  Possible decompositions: (10,0), (8,1), (6,2), (4,3), (2,4), (0,5)")

    for a1_target in [0, 2, 4, 6, 8, 10]:
        a2_needed = (10 - a1_target) // 2
        a2_set = sorted(a1_to_a2.get(a1_target, set()))
        found = a2_needed in a2_set

        # Find nearby a1 values if this one doesn't exist
        if a1_target not in a1_to_a2:
            print(f"\n  alpha_1={a1_target}: NOT ACHIEVABLE at n=7")
            # Check if achievable at n=6
            continue

        # Find the gap
        print(f"\n  alpha_1={a1_target}: achievable alpha_2 = {a2_set}")
        print(f"    Need alpha_2 = {a2_needed}: {'EXISTS' if found else 'IN GAP'}")

        if not found and a2_set:
            below = [x for x in a2_set if x < a2_needed]
            above = [x for x in a2_set if x > a2_needed]
            if below:
                print(f"    Largest alpha_2 below: {max(below)}")
            if above:
                print(f"    Smallest alpha_2 above: {min(above)}")

        # Show decomposition examples
        if a1_target in a1_to_decomp:
            for a2, a1_3, a1_5, a1_7, a2_33, a2_35, a2_other in a1_to_decomp[a1_target][:2]:
                print(f"    Example: a1={a1_target} = {a1_3}(3c) + {a1_5}(5c) + {a1_7}(7c)")
                print(f"             a2={a2} = {a2_33}(3-3) + {a2_35}(3-5) + {a2_other}(other)")

    # Full (alpha_1, alpha_2) table for a1 in [0..12]
    print(f"\n{'='*70}")
    print("FULL (alpha_1, alpha_2) ACHIEVABILITY TABLE (a1 = 0..14)")
    print("=" * 70)
    for a1 in range(15):
        if a1 in a1_to_a2:
            a2_set = sorted(a1_to_a2[a1])
            a2_needed = (10 - a1) / 2
            marker = ""
            if a1 <= 10 and a1 % 2 == 0 and int(a2_needed) == a2_needed:
                a2n = int(a2_needed)
                if a2n not in a2_set and a2n >= 0:
                    marker = f"  *** NEED a2={a2n} FOR H=21"
            print(f"  a1={a1:3d}: a2 in {a2_set}{marker}")
        else:
            marker = ""
            if a1 <= 10 and a1 % 2 == 0:
                marker = f"  *** a1 NOT ACHIEVABLE (blocks H=21)"
            print(f"  a1={a1:3d}: NOT ACHIEVED{marker}")

    # Check a1=3 specifically (connected to H=7)
    print(f"\n{'='*70}")
    print("alpha_1 = 3 (H=7 ROOT CAUSE)")
    print("=" * 70)
    a1_3_exists = 3 in a1_to_a2
    print(f"  alpha_1 = 3 achievable at n=7? {a1_3_exists}")
    if a1_3_exists:
        print(f"  alpha_2 values: {sorted(a1_to_a2[3])}")
        print(f"  NOTE: H=7 needs (a1=3, a2=0). If a2=0 not in set, H=7 still blocked.")

    # Summary of the six-way block
    print(f"\n{'='*70}")
    print("SIX-WAY BLOCK SUMMARY FOR H=21 AT n=7")
    print("=" * 70)
    blocked = 0
    for a1 in [0, 2, 4, 6, 8, 10]:
        a2_needed = (10 - a1) // 2
        a2_set = a1_to_a2.get(a1, set())
        if a1 not in a1_to_a2:
            reason = "alpha_1 not achievable"
        elif a2_needed not in a2_set:
            reason = f"alpha_2={a2_needed} not achievable (have {sorted(a2_set)})"
        else:
            reason = "ACHIEVABLE — H=21 SHOULD EXIST!"
        status = "BLOCKED" if a2_needed not in a2_set or a1 not in a1_to_a2 else "OPEN"
        if status == "BLOCKED":
            blocked += 1
        print(f"  ({a1},{a2_needed}): {status} — {reason}")
    print(f"\n  Result: {blocked}/6 decompositions blocked")
    if blocked == 6:
        print(f"  H=21 IS BLOCKED AT n=7 by six-way mechanism!")
    else:
        print(f"  WARNING: some decomposition is achievable!")

if __name__ == "__main__":
    main()
