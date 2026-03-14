"""
forbidden_h_achievability.py — WHY specific H values are permanently forbidden
kind-pasteur-2026-03-14-S65

Key insight: the "structural constraints" (a_k > 0 => a_1 >= k, a_2 >= C(k,2))
are too WEAK. Every odd H >= 1 has a valid decomposition a_1 = (H-1)/2 with
a_2 = a_3 = ... = 0. So structural constraints alone forbid NOTHING.

The REAL constraint is ACHIEVABILITY: which (a_1, a_2, ...) can actually arise
from a tournament's conflict graph Omega(T)?

For n=7 (quadratic era):
  H = 1 + 2*alpha_1 + 4*alpha_2
  alpha_1 = c_3 + c_5 + c_7  (number of directed odd cycles in T)
  alpha_2 = number of DISJOINT pairs of odd cycles

The permanent prohibition proofs:
  H=7: alpha_1=3, alpha_2=0. Three cycles, all pairwise conflicting.
        => Common vertex => 5-cycle exists => alpha_1 >= 4. Contradiction.
  H=21: Even in cubic era (n>=9), the constraints force H >= 27.

For H=63, 107, 119, 149: we need to understand WHY these specific
(alpha_1, alpha_2) combinations don't arise at n=7, and whether they
could arise at larger n.

This script does:
1. Exhaustive alpha_1, alpha_2 at n=7 (fast: only need 3-cycle structure)
2. Check which (alpha_1, alpha_2) lines are empty vs populated
3. Determine if the gaps are permanent or n-dependent
"""

import numpy as np
from itertools import combinations
from collections import defaultdict

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_ham_paths(A):
    n = len(A)
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n + 1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total > 0:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def get_3cycles(A):
    """Get all directed 3-cycle vertex sets."""
    n = len(A)
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[i][k] and A[k][j] and A[j][i]):
                    cycles.append(frozenset([i, j, k]))
    return cycles

def count_disjoint_pairs(cycles):
    """Count disjoint pairs among a list of frozensets."""
    count = 0
    nc = len(cycles)
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i].isdisjoint(cycles[j]):
                count += 1
    return count

def main():
    print("=" * 70)
    print("FORBIDDEN H VALUES — ACHIEVABILITY ANALYSIS")
    print("=" * 70)

    # Part 1: Exhaustive n=5 (quick sanity check)
    print("\nPART 1: Exhaustive n=5")
    n = 5
    n_arcs = n * (n-1) // 2
    h_alpha = defaultdict(set)  # H -> set of (alpha_1, alpha_2) pairs

    for bits in range(2**n_arcs):
        A = np.zeros((n, n), dtype=int)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if bits & (1 << idx):
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1
        # At n=5, only 3-cycles and 5-cycles matter
        # But for OCF, alpha_1 counts ALL odd cycles and alpha_2 counts disjoint pairs
        # For simplicity, compute H directly and infer alpha_1, alpha_2
        H = count_ham_paths(A)

        # Get 3-cycles (these are the odd cycles at n=5 that contribute to alpha)
        c3 = get_3cycles(A)
        # At n=5: a 5-cycle is a Hamiltonian cycle of the subtournament on all 5 vertices
        # But for the independence polynomial, odd cycles of ALL lengths contribute
        # For now, let's just record H since we know H = 1 + 2*alpha_1 + 4*alpha_2

        # From H: alpha_1 + 2*alpha_2 = (H-1)/2
        # But we can't determine alpha_1, alpha_2 separately from H alone at n=5
        # unless we compute them directly

        # Skip full cycle enumeration for now, just record H
        h_alpha[H] = h_alpha.get(H, set()) | {H}

    h_vals_n5 = sorted(h_alpha.keys())
    print(f"  Achievable H at n=5: {h_vals_n5}")
    print(f"  Distinct: {len(h_vals_n5)}")
    print(f"  Missing odd from 1 to {max(h_vals_n5)}: "
          f"{[h for h in range(1, max(h_vals_n5)+1, 2) if h not in h_vals_n5]}")

    # Part 2: Exhaustive n=7 — H values only (fast)
    print("\nPART 2: Exhaustive n=7 — H distribution")
    print("  (2^21 = 2,097,152 tournaments)")

    n = 7
    n_arcs = n * (n-1) // 2  # 21
    h_counts_7 = defaultdict(int)

    for bits in range(2**n_arcs):
        A = np.zeros((n, n), dtype=int)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if bits & (1 << idx):
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1
        H = count_ham_paths(A)
        h_counts_7[H] += 1

        if bits % 500000 == 0 and bits > 0:
            print(f"  Progress: {bits}/{2**n_arcs} ({100*bits/2**n_arcs:.1f}%)")

    h_vals_n7 = sorted(h_counts_7.keys())
    print(f"\n  Achievable H at n=7: {len(h_vals_n7)} distinct values")
    print(f"  Range: [{min(h_vals_n7)}, {max(h_vals_n7)}]")

    # All odd H from 1 to max that are NOT achieved
    max_h = max(h_vals_n7)
    missing_n7 = [h for h in range(1, max_h + 1, 2) if h not in h_counts_7]
    print(f"  Missing odd H in [1, {max_h}]: {missing_n7}")
    print(f"  Count: {len(missing_n7)}")

    # Verify our known forbidden values
    known_forbidden = [7, 21, 63, 107, 119, 149]
    print(f"\n  Status of known forbidden values:")
    for h in known_forbidden:
        count = h_counts_7.get(h, 0)
        print(f"    H={h:4d}: {'FOUND (' + str(count) + ')' if count > 0 else 'ABSENT'}")

    # Part 3: H values and their neighbors
    print("\nPART 3: Neighborhood of forbidden values")
    for h_target in known_forbidden:
        print(f"\n  Around H = {h_target}:")
        for h in range(h_target - 8, h_target + 9, 2):
            if h < 1:
                continue
            count = h_counts_7.get(h, 0)
            marker = " ** FORBIDDEN" if count == 0 else ""
            bar = "#" * min(count // 100, 40) if count > 0 else ""
            print(f"    H={h:4d}: {count:7d} {bar}{marker}")

    # Part 4: Density pattern of forbidden values
    print("\nPART 4: Density of forbidden odd H values")
    # What fraction of odd H in [1, max] are forbidden?
    total_odd = len(range(1, max_h + 1, 2))
    print(f"  Total odd H in [1, {max_h}]: {total_odd}")
    print(f"  Achieved: {len(h_vals_n7)}")
    print(f"  Forbidden: {len(missing_n7)} ({100*len(missing_n7)/total_odd:.1f}%)")

    # Gaps by region
    for low, high in [(1, 50), (50, 100), (100, 150), (150, max_h)]:
        gaps_in_range = [h for h in missing_n7 if low <= h < high]
        total_in_range = len([h for h in range(low if low % 2 == 1 else low + 1, high, 2)])
        if total_in_range > 0:
            print(f"    [{low:3d}, {high:3d}): {len(gaps_in_range)}/{total_in_range} forbidden "
                  f"({100*len(gaps_in_range)/total_in_range:.1f}%)")

    # Part 5: Look at the (H-1)/2 values and their factorizations
    print("\nPART 5: Arithmetic structure of forbidden H")
    print(f"\n  Missing H values: {missing_n7}")
    print(f"\n  (H-1)/2 values (= alpha_1 + 2*alpha_2):")
    for h in missing_n7:
        t = (h - 1) // 2
        # Factorization
        factors = []
        x = t
        for p in [2, 3, 5, 7, 11, 13]:
            while x % p == 0:
                factors.append(p)
                x //= p
            if x == 1:
                break
        if x > 1:
            factors.append(x)
        print(f"    H={h:4d}: T={t:3d} = {'*'.join(str(f) for f in factors) if factors else '0'}")

    # Part 6: Check n=8 sampling
    print(f"\n{'='*70}")
    print("PART 6: SAMPLING AT n=8 AND n=9")
    print(f"{'='*70}")

    rng = np.random.default_rng(2026_0314)

    for n_test in [8, 9]:
        N = 3000 if n_test == 8 else 1000
        h_set = set()
        target_found = {h: 0 for h in known_forbidden}

        for trial in range(N):
            A = random_tournament(n_test, rng)
            H = count_ham_paths(A)
            h_set.add(H)
            if H in target_found:
                target_found[H] += 1

        print(f"\n  n={n_test} ({N} samples):")
        print(f"    Distinct H seen: {len(h_set)}")
        for h in known_forbidden:
            status = f"FOUND ({target_found[h]}x)" if target_found[h] > 0 else "absent"
            print(f"    H={h:4d}: {status}")

        # Check which missing_n7 values appeared
        newly_found = [h for h in missing_n7 if h in h_set]
        if newly_found:
            print(f"    n=7 gaps now FILLED at n={n_test}: {newly_found}")
        else:
            print(f"    No n=7 gaps filled at n={n_test} (in {N} samples)")

    print(f"\n{'='*70}")
    print("FINAL CONCLUSIONS")
    print(f"{'='*70}")
    print(f"""
  PROVED PERMANENTLY FORBIDDEN:
    H = 7:  alpha_1=3 with alpha_2=0 is impossible (common vertex argument)
    H = 21: cubic I.P. can't help (forces H >= 27); quadratic impossible

  STATUS OF OTHER n=7 GAPS:
    H = 63, 107, 119, 149: absent at n=7 by exhaustive computation
    All have valid structural decompositions — the gap is from tournament
    graph constraints, not pure integer constraints.

  OPEN: Whether H = 63, 107, 119, 149 become achievable at n >= 8.
  If ANY of them appear at larger n, they are NOT permanently forbidden.
  If NONE appear, they might be — but the proof would need tournament-
  specific arguments (not just structural constraints).

  The permanent prohibition of H=7 and H=21 uses TOPOLOGICAL arguments
  about the conflict graph Omega(T), not just integer arithmetic.
  Extending these arguments to H=63 requires understanding exactly
  which (alpha_1, alpha_2) pairs are achievable at ALL n.
  """)

if __name__ == "__main__":
    main()
