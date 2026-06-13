"""
forbidden_boundary.py — Map the exact boundary of achievable (a1, a2) at n=7
kind-pasteur-2026-03-14-S64

Key questions:
1. What (a1, a2) pairs exist with a1+2*a2 = 9 and a1+2*a2 = 11?
2. What prevents a1+2*a2 = 10?
3. Is there a parity constraint on a1+2*a2?
4. What is the EXACT achievable set?

Strategy: Use H directly (much faster than cycle decomposition).
For targeted (a1,a2) analysis, use smaller samples with full decomposition.
"""

import numpy as np
from itertools import combinations
from collections import defaultdict

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

def get_full_data(A):
    """Get (c3, c5, c7, a1, a2)."""
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

    # 7-cycles
    c7 = count_ham_cycles_sub(A, 7) if n == 7 else 0

    a1 = c3 + c5 + c7

    # a2: disjoint pairs of 3-cycles (only possible disjoint odd cycle pairs at n=7)
    a2 = 0
    for i in range(len(cycles_3)):
        for j in range(i+1, len(cycles_3)):
            if cycles_3[i].isdisjoint(cycles_3[j]):
                a2 += 1

    return c3, c5, c7, a1, a2

def main():
    n = 7

    # ================================================================
    print("=" * 70)
    print("PART 1: Map achievable (a1, a2) with full decomposition")
    print("=" * 70)

    rng = np.random.default_rng(2026)
    N = 3000

    all_pairs = defaultdict(list)  # (a1, a2) -> [(c3, c5, c7), ...]
    k_to_pairs = defaultdict(set)  # k = a1+2*a2 -> set of (a1, a2)

    for trial in range(N):
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(i+1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        c3, c5, c7, a1, a2 = get_full_data(A.tolist())
        H = 1 + 2*a1 + 4*a2
        H_check = count_ham_paths(A)
        assert H == H_check, f"OCF mismatch: {H} != {H_check}"

        k = a1 + 2*a2
        all_pairs[(a1, a2)].append((c3, c5, c7))
        k_to_pairs[k].add((a1, a2))

        if (trial+1) % 500 == 0:
            print(f"  {trial+1}/{N}")

    # ================================================================
    print("\n" + "=" * 70)
    print("PART 2: The achievable k = a1+2*a2 values")
    print("=" * 70)

    print(f"\n  Achievable k values and their (a1, a2) fibers:")
    for k in sorted(k_to_pairs.keys()):
        pairs = sorted(k_to_pairs[k])
        h = 2*k + 1
        print(f"  k={k:3d} (H={h:3d}): {pairs}")

    # What k values are missing?
    all_k = sorted(k_to_pairs.keys())
    max_k = max(all_k)
    missing_k = [k for k in range(max_k+1) if k not in k_to_pairs]
    print(f"\n  Missing k values (up to {max_k}): {missing_k}")
    print(f"  Corresponding forbidden H: {[2*k+1 for k in missing_k]}")

    # ================================================================
    print("\n" + "=" * 70)
    print("PART 3: Parity structure of a1 and a2")
    print("=" * 70)

    # Check if a1 parity constrains a2 parity
    parity_dist = defaultdict(int)
    for (a1, a2), decomps in all_pairs.items():
        parity_dist[(a1 % 2, a2 % 2)] += len(decomps)

    print(f"\n  (a1 mod 2, a2 mod 2) distribution:")
    for (p1, p2), count in sorted(parity_dist.items()):
        print(f"    ({p1}, {p2}): {count}")

    # Check k = a1 + 2*a2 parity (always same parity as a1)
    print(f"\n  k = a1 + 2*a2 mod 2 = a1 mod 2 (since 2*a2 is even)")
    print(f"  So k is odd iff a1 is odd, k is even iff a1 is even.")

    # ================================================================
    print("\n" + "=" * 70)
    print("PART 4: Detailed structure near k=10 (H=21)")
    print("=" * 70)

    print(f"\n  For k=10 (H=21): need a1+2*a2=10")
    print(f"  Possible: (a1,a2) in {{(10,0), (8,1), (6,2), (4,3), (2,4), (0,5)}}")
    print(f"  All have even a1 (so k=10 even, consistent).")

    # Look at nearby achievable points
    for k_target in [8, 9, 10, 11, 12]:
        h = 2*k_target + 1
        if k_target in k_to_pairs:
            print(f"\n  k={k_target} (H={h}): ACHIEVABLE")
            for (a1, a2) in sorted(k_to_pairs[k_target]):
                decomps = all_pairs[(a1, a2)]
                c3_vals = [c3 for c3, c5, c7 in decomps]
                c5_vals = [c5 for c3, c5, c7 in decomps]
                c7_vals = [c7 for c3, c5, c7 in decomps]
                print(f"    (a1={a1}, a2={a2}): {len(decomps)} samples, "
                      f"c3 in [{min(c3_vals)},{max(c3_vals)}], "
                      f"c5 in [{min(c5_vals)},{max(c5_vals)}], "
                      f"c7 in [{min(c7_vals)},{max(c7_vals)}]")
        else:
            print(f"\n  k={k_target} (H={h}): FORBIDDEN")

    # ================================================================
    print("\n" + "=" * 70)
    print("PART 5: Why a1=10 with a2=0 is impossible")
    print("=" * 70)

    # a1 = c3 + c5 + c7 = 10, a2 = 0
    # a2 = 0 means no disjoint 3-cycle pairs
    # c5 is constrained by c3
    # Let's check: what c3 values allow a1=10?

    print(f"\n  a1 = c3 + c5 + c7 = 10 with a2 = 0")
    print(f"  Need: no two disjoint 3-cycles AND c3+c5+c7=10")
    print()

    # From data: what's the relationship between c3 and c5?
    c3_c5_map = defaultdict(list)
    c3_a2_map = defaultdict(list)
    for (a1, a2), decomps in all_pairs.items():
        for c3, c5, c7 in decomps:
            c3_c5_map[c3].append(c5)
            c3_a2_map[c3].append(a2)

    print(f"  c3 -> min c5 (determines min a1 for given c3):")
    for c3 in sorted(c3_c5_map.keys()):
        c5_list = c3_c5_map[c3]
        a2_list = c3_a2_map[c3]
        min_c5 = min(c5_list)
        max_a2_0 = sum(1 for a in a2_list if a == 0)
        total = len(c5_list)
        print(f"    c3={c3:2d}: c5 in [{min(c5_list):2d}, {max(c5_list):2d}], "
              f"min a1 = {c3+min_c5} (or +1 if c7=1), "
              f"a2=0 rate: {max_a2_0}/{total} = {100*max_a2_0/total:.0f}%")

    print(f"\n  For a1=10 with c7=0: need c3+c5=10")
    print(f"  c3=0: min c5=0, c3+c5=0 (too small)")
    print(f"  c3=1: min c5=0, need c5=9 (max c5 for c3=1?)")
    print(f"  c3=4: need c5=6 (max c5 for c3=4?)")
    print(f"  c3=6: need c5=4 (min c5 for c3=6?)")

    # Check: for c3=6, is c5=4 achievable?
    c5_at_c3_6 = c3_c5_map.get(6, [])
    if c5_at_c3_6:
        print(f"\n  c3=6: c5 values observed: {sorted(set(c5_at_c3_6))}")
        print(f"  c5=4 achievable for c3=6: {4 in c5_at_c3_6}")

    # For c3=5, need c5=5
    c5_at_c3_5 = c3_c5_map.get(5, [])
    if c5_at_c3_5:
        print(f"  c3=5: c5 values observed: {sorted(set(c5_at_c3_5))}")
        print(f"  c5=5 achievable for c3=5: {5 in c5_at_c3_5}")

    # ================================================================
    print("\n" + "=" * 70)
    print("PART 6: The c3-c5-a2 joint constraint")
    print("=" * 70)

    # For c3 where a2=0 is possible AND c3+c5 could = 10:
    # Check if c5 is LIMITED when a2=0

    print("\n  When a2=0 (no disjoint 3-cycle pairs), what c5 values occur?")
    c3_c5_when_a2_0 = defaultdict(list)
    for (a1, a2), decomps in all_pairs.items():
        if a2 == 0:
            for c3, c5, c7 in decomps:
                c3_c5_when_a2_0[c3].append(c5)

    for c3 in sorted(c3_c5_when_a2_0.keys()):
        c5_list = c3_c5_when_a2_0[c3]
        a1_list = [c3 + c5 for c5 in c5_list]  # ignoring c7 for now
        print(f"    c3={c3:2d} (a2=0): c5 in [{min(c5_list):2d}, {max(c5_list):2d}], "
              f"a1 = c3+c5 in [{min(a1_list)}, {max(a1_list)}], count={len(c5_list)}")

    print(f"\n  Can c3+c5 = 10 when a2=0?")
    for c3 in sorted(c3_c5_when_a2_0.keys()):
        needed_c5 = 10 - c3
        if needed_c5 >= 0:
            c5_list = c3_c5_when_a2_0[c3]
            achievable = needed_c5 in c5_list
            min_c5 = min(c5_list) if c5_list else "N/A"
            max_c5 = max(c5_list) if c5_list else "N/A"
            if needed_c5 <= max(c5_list) + 5:  # reasonable range
                print(f"    c3={c3:2d}: need c5={needed_c5}, range=[{min_c5},{max_c5}], "
                      f"achievable: {achievable}")

    # ================================================================
    print("\n" + "=" * 70)
    print("PART 7: Score sequence constraint for a1=10")
    print("=" * 70)

    print("""
  At n=7, c3 = 35 - sum C(s_i, 2) where s_i are out-degrees.
  sum s_i = C(7,2) = 21.

  For c3 = 5: sum C(s_i, 2) = 30
  For c3 = 6: sum C(s_i, 2) = 29
  For c3 = 7: sum C(s_i, 2) = 28
  For c3 = 8: sum C(s_i, 2) = 27

  a2=0 requires all 3-cycles to be pairwise conflicting (share a vertex).
  At c3=k with a2=0: the k cyclic triples form a clique in the conflict graph.

  Maximum clique: at most C(n-1, 2) = 15 cycles through a common vertex.
  But not all triples through a vertex need be cyclic.

  KEY QUESTION: for tournaments where a2=0 and c3 is small (5-8),
  what are the achievable c5 values?
  The RIGID relationship between c3, c5, and the score sequence
  might prevent c3+c5 = 10.
""")

    # Check: for each a2=0 tournament, what is the exact a1?
    a1_when_a2_0 = []
    for (a1, a2), decomps in all_pairs.items():
        if a2 == 0:
            a1_when_a2_0.extend([a1] * len(decomps))

    if a1_when_a2_0:
        print(f"  a1 values when a2=0: {sorted(set(a1_when_a2_0))}")
        print(f"  a1=10 when a2=0: {a1_when_a2_0.count(10)} occurrences")
        print(f"  Nearest to 10: {min(a1_when_a2_0, key=lambda x: abs(x-10))}")

    # What about a2 > 0?
    print(f"\n  For a1+2*a2=10 with a2>0:")
    for a2_val in range(1, 6):
        a1_val = 10 - 2*a2_val
        if a1_val >= 0:
            count = len(all_pairs.get((a1_val, a2_val), []))
            print(f"    (a1={a1_val}, a2={a2_val}): {count} occurrences")

    # ================================================================
    print("\n" + "=" * 70)
    print("PART 8: The achievable a1 values at each a2")
    print("=" * 70)

    a2_to_a1 = defaultdict(set)
    for (a1, a2), decomps in all_pairs.items():
        a2_to_a1[a2].add(a1)

    for a2_val in sorted(a2_to_a1.keys()):
        a1_set = sorted(a2_to_a1[a2_val])
        # Check gaps
        if len(a1_set) >= 2:
            gaps = [a1_set[i+1] - a1_set[i] for i in range(len(a1_set)-1)]
            max_gap = max(gaps)
        else:
            max_gap = 0
        print(f"  a2={a2_val:2d}: a1 in {a1_set[:15]}{'...' if len(a1_set)>15 else ''}, "
              f"count={len(a1_set)}, max_gap={max_gap}")

    print("\nDone.")

if __name__ == "__main__":
    main()
