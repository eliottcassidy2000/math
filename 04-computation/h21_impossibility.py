"""
h21_impossibility.py — Why is H=21 impossible?
kind-pasteur-2026-03-14-S64

H=21 = 1 + 2*a1 + 4*a2 requires a1 + 2*a2 = 10.
The 6 candidate (a1,a2) pairs are (10,0), (8,1), (6,2), (4,3), (2,4), (0,5).
NONE of these are achievable at n=7. WHY?

Key insight: a1 counts ALL directed odd cycles (3, 5, 7),
and a2 counts pairs of vertex-disjoint odd cycles.

At n=7, disjoint pairs can only come from two 3-cycles (3+3=6 <= 7).
So a2 = number of disjoint pairs of 3-CYCLES only.

The constraint: c3 + c5 + c7 = a1, where c_k = # directed k-cycles.

We know: c3 = C(7,3) - #{transitive triples} = 35 - trans_triples
c7 in {0, 1} (at most one directed Hamiltonian cycle on all 7 vertices)
c5 depends on the tournament structure.

For a1 = 10 with a2 = 0:
Need exactly 10 total odd cycles with NO disjoint pairs of 3-cycles.
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
from math import comb

def all_tournaments(n):
    """Generate all tournaments on n vertices."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A

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

def count_ham_cycles_sub(A_sub, m):
    """Count directed Hamiltonian cycles in m-vertex sub-tournament."""
    if m < 3:
        return 0
    dp = {}
    dp[(1 << 0, 0)] = 1
    for mask_size in range(2, m + 1):
        for mask in range(1 << m):
            if bin(mask).count('1') != mask_size:
                continue
            if not (mask & 1):
                continue
            for v in range(m):
                if v == 0 and mask_size < m:
                    continue
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(m):
                    if (prev_mask & (1 << u)) and A_sub[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total > 0:
                    dp[(mask, v)] = dp.get((mask, v), 0) + total
    full = (1 << m) - 1
    count = 0
    for u in range(1, m):
        if A_sub[u][0]:
            count += dp.get((full, u), 0)
    return count

def get_cycle_data(A):
    """Return (c3, c5, c7, a1, a2, list_of_3cycle_sets)."""
    n = len(A)

    # 3-cycles
    cycles_3 = []
    for verts in combinations(range(n), 3):
        i, j, k = verts
        if A[i][j] and A[j][k] and A[k][i]:
            cycles_3.append(frozenset(verts))
        elif A[i][k] and A[k][j] and A[j][i]:
            cycles_3.append(frozenset(verts))
    c3 = len(cycles_3)

    # 5-cycles
    c5 = 0
    for verts in combinations(range(n), 5):
        sub = [[A[verts[i]][verts[j]] for j in range(5)] for i in range(5)]
        c5 += count_ham_cycles_sub(sub, 5)

    # 7-cycles (only if n=7)
    c7 = 0
    if n == 7:
        c7 = count_ham_cycles_sub(A, 7)

    a1 = c3 + c5 + c7

    # a2: disjoint pairs of odd cycles
    # At n=7, only 3-cycle pairs can be disjoint
    a2 = 0
    for i in range(len(cycles_3)):
        for j in range(i+1, len(cycles_3)):
            if cycles_3[i].isdisjoint(cycles_3[j]):
                a2 += 1

    return c3, c5, c7, a1, a2, cycles_3

def main():
    # ================================================================
    print("=" * 70)
    print("PART 1: Exhaustive analysis at n=5 and n=6")
    print("=" * 70)

    for n in [5, 6]:
        print(f"\n  --- n={n} ---")
        all_a1_a2 = defaultdict(int)
        h_values = defaultdict(int)
        c_decomp = defaultdict(set)

        edges = [(i, j) for i in range(n) for j in range(i+1, n)]
        m = len(edges)
        total = 1 << m

        for bits in range(total):
            A = [[0]*n for _ in range(n)]
            for k, (i, j) in enumerate(edges):
                if bits & (1 << k):
                    A[i][j] = 1
                else:
                    A[j][i] = 1

            c3, c5, c7, a1, a2, _ = get_cycle_data(A)
            H = 1 + 2*a1 + 4*a2
            H_actual = count_ham_paths(A)
            assert H == H_actual, f"OCF mismatch at n={n}: {H} != {H_actual}"

            all_a1_a2[(a1, a2)] += 1
            h_values[H] += 1
            c_decomp[a1].add((c3, c5, c7))

        # H=21 check
        print(f"  Total tournaments: {total}")
        if 21 in h_values:
            print(f"  H=21: {h_values[21]} tournaments")
        else:
            print(f"  H=21: IMPOSSIBLE at n={n}")

        # What (a1, a2) are achievable?
        # For H=21: need a1 + 2*a2 = 10
        print(f"  Checking H=21 candidates (a1 + 2*a2 = 10):")
        for a2_test in range(6):
            a1_test = 10 - 2*a2_test
            if (a1_test, a2_test) in all_a1_a2:
                print(f"    ({a1_test}, {a2_test}): ACHIEVABLE ({all_a1_a2[(a1_test, a2_test)]} tours)")
            else:
                print(f"    ({a1_test}, {a2_test}): IMPOSSIBLE")

        # What a1 values are achievable?
        a1_range = sorted(set(a for a, _ in all_a1_a2))
        print(f"  Achievable a1 range: {a1_range}")

        # Focus: is a1=10 achievable?
        a1_10_pairs = [(a, b) for (a, b) in all_a1_a2 if a == 10]
        if a1_10_pairs:
            print(f"  a1=10 achieved with a2 in: {[b for _, b in sorted(a1_10_pairs)]}")
        else:
            print(f"  a1=10: IMPOSSIBLE at n={n}")

        # What's the minimum achievable a1 at each a2?
        a2_range = sorted(set(b for _, b in all_a1_a2))
        print(f"  Achievable a2 range: {a2_range}")
        for a2_val in a2_range:
            a1_vals = sorted(a for a, b in all_a1_a2 if b == a2_val)
            print(f"    a2={a2_val}: a1 in {a1_vals}")

    # ================================================================
    print("\n" + "=" * 70)
    print("PART 2: The (c3, c5, c7) decomposition at n=7")
    print("=" * 70)

    print("\n  At n=7, a1 = c3 + c5 + c7.")
    print("  c3 = number of cyclic triples (3-cycles)")
    print("  c5 = number of directed 5-cycles")
    print("  c7 = 0 or 1 (existence of directed Hamiltonian cycle)")
    print()
    print("  For H=21, need a1 + 2*a2 = 10.")
    print("  a2 = disjoint pairs of 3-cycles.")
    print()
    print("  Sampling n=7 tournaments to find (a1, a2) near H=21...")

    rng = np.random.default_rng(42)
    n = 7

    # Collect comprehensive statistics
    all_data = defaultdict(list)  # (a1, a2) -> list of (c3, c5, c7)
    h_data = defaultdict(int)
    near_21 = []

    N = 5000
    for trial in range(N):
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(i+1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        c3, c5, c7, a1, a2, cycles_3 = get_cycle_data(A.tolist())
        H = 1 + 2*a1 + 4*a2
        h_data[H] += 1
        all_data[(a1, a2)].append((c3, c5, c7))

        if abs(H - 21) <= 10:
            near_21.append((H, a1, a2, c3, c5, c7))

        if (trial+1) % 1000 == 0:
            print(f"  {trial+1}/{N}")

    # Tournaments with H near 21
    print(f"\n  Tournaments with |H-21| <= 10:")
    near_21.sort()
    for H, a1, a2, c3, c5, c7 in near_21[:30]:
        print(f"    H={H:3d}: a1={a1:2d}, a2={a2:2d}, c3={c3:2d}, c5={c5:2d}, c7={c7}")

    # H gap analysis
    print(f"\n  H values near 21:")
    for h in range(11, 32):
        count = h_data.get(h, 0)
        marker = " <-- GAP" if count == 0 else ""
        print(f"    H={h:3d}: {count:5d} tournaments{marker}")

    # ================================================================
    print("\n" + "=" * 70)
    print("PART 3: Why a1=10 at n=7 is trapped")
    print("=" * 70)

    # For a1=10, what are the possible (c3, c5, c7) decompositions?
    print("\n  a1 = c3 + c5 + c7 = 10")
    print("  c7 in {0, 1}")
    print("  Possible decompositions: (c3, c5, c7) with c3+c5+c7=10")

    # At n=7: c3 ranges from 0 (transitive) to 35 (= C(7,3))
    # But c3 has constraints based on score sequence
    # c3 = C(7,3) - sum C(s_i, 2) where s_i are out-degrees

    print("\n  Score sequence constraint: c3 = 35 - sum C(s_i, 2)")
    print("  For c3=10: sum C(s_i, 2) = 25")
    print("  For c3=9: sum C(s_i, 2) = 26")
    print("  For c3=8: sum C(s_i, 2) = 27")

    # What score sequences give c3 = 8, 9, 10?
    from math import comb as C
    for target_c3 in [8, 9, 10]:
        target_sum = 35 - target_c3
        print(f"\n  c3={target_c3}: need sum C(s_i,2) = {target_sum}")
        # Enumerate score sequences (s_0 <= s_1 <= ... <= s_6, sum = C(7,2) = 21)
        found_scores = []
        for s in range(7):
            for s1 in range(s, 7):
                for s2 in range(s1, 7):
                    for s3 in range(s2, 7):
                        for s4 in range(s3, 7):
                            for s5 in range(s4, 7):
                                s6 = 21 - s - s1 - s2 - s3 - s4 - s5
                                if s6 < s5 or s6 >= 7:
                                    continue
                                scores = [s, s1, s2, s3, s4, s5, s6]
                                if sum(C(si, 2) for si in scores) == target_sum:
                                    found_scores.append(tuple(scores))
        if found_scores:
            for sc in found_scores:
                print(f"    Score sequence: {sc}, c3={target_c3}")
        else:
            print(f"    NO valid score sequence!")

    # ================================================================
    print("\n" + "=" * 70)
    print("PART 4: The c5 constraint — what limits 5-cycles?")
    print("=" * 70)

    # For tournaments near H=21, what c5 values appear?
    # Group by c3
    c3_to_c5 = defaultdict(list)
    for (a1, a2), decomps in all_data.items():
        for c3, c5, c7 in decomps:
            c3_to_c5[c3].append((c5, c7, a1, a2))

    print(f"\n  c5 statistics by c3 (n=7, {N} samples):")
    for c3_val in sorted(c3_to_c5.keys()):
        c5_vals = [c5 for c5, _, _, _ in c3_to_c5[c3_val]]
        if len(c5_vals) > 0:
            print(f"    c3={c3_val:2d}: c5 in [{min(c5_vals):2d}, {max(c5_vals):2d}], "
                  f"mean={np.mean(c5_vals):.1f}, count={len(c5_vals)}")

    # For c3 near 8-10, what a2 values appear?
    print(f"\n  a2 by c3 (for c3 in 8-12):")
    for c3_val in range(8, 13):
        if c3_val in c3_to_c5:
            a2_vals = [a2 for _, _, _, a2 in c3_to_c5[c3_val]]
            print(f"    c3={c3_val:2d}: a2 in {sorted(set(a2_vals))}, count={len(a2_vals)}")

    # ================================================================
    print("\n" + "=" * 70)
    print("PART 5: The RIGID constraint — c3 determines a2 range")
    print("=" * 70)

    # Key question: for c3=k, what a2 values are possible?
    # a2 = number of disjoint pairs of 3-cycles
    # If c3 cycles are {S_1, ..., S_c3}, a2 = number of vertex-disjoint pairs

    # For a1=10 and a2=0: c3+c5+c7=10, no disjoint 3-cycle pairs
    # This means ALL 3-cycles share at least one vertex
    # With c3 triples on 7 vertices, this is a CLIQUE condition on the conflict graph

    print("\n  a2=0 requires: every pair of 3-cycles shares a vertex.")
    print("  This means the 3-cycles form a CLIQUE in the conflict graph.")
    print()
    print("  Maximum clique in conflict graph of 3-cycles on 7 vertices:")
    print("  A 3-cycle uses 3 of 7 vertices.")
    print("  Two 3-cycles conflict iff they share >= 1 vertex.")
    print("  Two 3-cycles are DISJOINT only if they share 0 vertices (using 6 of 7).")
    print()
    print("  For c3 cycles with a2=0: ALL cycles pairwise conflict.")
    print("  Max # of pairwise conflicting 3-cycles on 7 vertices:")
    print("  If they all share vertex v, up to C(6,2)=15 cycles through v.")
    print("  But not all of those need to be 3-cycles in the tournament.")
    print()
    print("  KEY: if c3=10 and a2=0, all 10 cyclic triples share a common vertex.")
    print("  But with 10 triples on 7 vertices, Helly-type argument...")

    # Check: at c3 >= 8, does a2=0 still occur?
    c3_a2_zero = defaultdict(int)
    c3_total = defaultdict(int)
    for (a1, a2), decomps in all_data.items():
        for c3, c5, c7 in decomps:
            c3_total[c3] += 1
            if a2 == 0:
                c3_a2_zero[c3] += 1

    print(f"\n  Frequency of a2=0 by c3:")
    for c3_val in sorted(c3_total.keys()):
        zero_count = c3_a2_zero.get(c3_val, 0)
        total = c3_total[c3_val]
        pct = 100*zero_count/total if total > 0 else 0
        print(f"    c3={c3_val:2d}: a2=0 in {zero_count:4d}/{total:4d} ({pct:5.1f}%)")

    # ================================================================
    print("\n" + "=" * 70)
    print("PART 6: Exhaustive check — does (a1=10, a2=0) exist at n=7?")
    print("=" * 70)

    # This is hard to do exhaustively (2^21 = 2M tournaments)
    # But we can try MANY more random samples
    print(f"\n  Checking 20000 random tournaments for (a1=10, a2=0)...")
    rng2 = np.random.default_rng(12345)
    found_10_0 = 0
    found_h21 = 0
    h_near = defaultdict(int)

    for trial in range(20000):
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(i+1, n):
                if rng2.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        c3, c5, c7, a1, a2, _ = get_cycle_data(A.tolist())
        H = 1 + 2*a1 + 4*a2

        if a1 == 10 and a2 == 0:
            found_10_0 += 1
        if H == 21:
            found_h21 += 1
        if 15 <= H <= 25:
            h_near[H] += 1

        if (trial+1) % 5000 == 0:
            print(f"  {trial+1}/20000")

    print(f"\n  (a1=10, a2=0) found: {found_10_0}")
    print(f"  H=21 found: {found_h21}")
    print(f"\n  H distribution near 21 (20000 samples):")
    for h in range(15, 26):
        print(f"    H={h:3d}: {h_near.get(h, 0):5d}{'  <-- GAP' if h_near.get(h, 0) == 0 else ''}")

    # ================================================================
    print("\n" + "=" * 70)
    print("PART 7: The (a1, a2) constraint graph — what IS achievable?")
    print("=" * 70)

    # Merge all data
    all_merged = defaultdict(int)
    for trial in range(20000):
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(i+1, n):
                if np.random.RandomState(trial + 50000).random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        c3, c5, c7, a1, a2, _ = get_cycle_data(A.tolist())
        all_merged[(a1, a2)] += 1

    # Near the H=21 line: a1 + 2*a2 = 10
    print(f"\n  Points near the line a1 + 2*a2 = 10:")
    for target in [8, 9, 10, 11, 12]:
        points = [(a1, a2, c) for (a1, a2), c in all_merged.items() if a1 + 2*a2 == target]
        points.sort()
        label = " <-- H=21" if target == 10 else ""
        if points:
            print(f"    a1+2a2={target:2d}: {points}{label}")
        else:
            print(f"    a1+2a2={target:2d}: EMPTY{label}")

    print("\nDone.")

if __name__ == "__main__":
    main()
