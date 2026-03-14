#!/usr/bin/env python3
"""
forbidden_H_analysis.py -- Why are certain H values impossible?

The gaps {7, 21, 63, 107, 119, 127, 147, 149, 157} at n=7 cannot be
explained by any modular constraint (all odd residues mod 72 are achieved).
The gaps come from constraints on achievable (a1, a2) pairs.

H = 1 + 2*a1 + 4*a2. For given H, the achievable pairs are:
  a1 = (H-1-4*a2)/2 for each valid a2

So the question is: for which a2, does a1 = (H-1-4*a2)/2 actually occur?

This script exhaustively catalogs:
1. The full achievable (a1, a2) lattice at n=5,6,7
2. The fiber H -> set of (a1,a2) pairs
3. The gap structure and what makes a gap
4. The relationship between gaps and the a2/a1 ratio
5. Does the gap at H=7 persist for all n?
6. Does H=21 ever become achievable?

Author: kind-pasteur-2026-03-14-S63
"""

import numpy as np
from itertools import combinations
from collections import defaultdict


def bits_to_adj_list(bits, n):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
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
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_ham_cycles_exact(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp or dp[(mask, v)] == 0:
                continue
            cnt = dp[(mask, v)]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and dp[(full, v)] > 0:
            if A[verts[v]][verts[0]]:
                total += dp[(full, v)]
    return total


def get_a1_a2(A, n):
    cycles = []
    for k in range(3, n+1, 2):
        for subset in combinations(range(n), k):
            verts = list(subset)
            nc = count_ham_cycles_exact(A, verts)
            for _ in range(nc):
                cycles.append(frozenset(subset))
    a1 = len(cycles)
    a2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if not (cycles[i] & cycles[j]):
                a2 += 1
    return a1, a2


def main():
    # ============================================================
    # PART 1: Exhaustive at n=5
    # ============================================================
    print("=" * 70)
    print("PART 1: Exhaustive (a1, a2, H) catalog at n=5")
    print("=" * 70)

    n = 5
    total_bits = n*(n-1)//2

    all_data_5 = []
    for bits in range(2**total_bits):
        A = bits_to_adj_list(bits, n)
        H = count_ham_paths(A, n)
        a1, a2 = get_a1_a2(A, n)
        all_data_5.append((a1, a2, H))

    pairs_5 = sorted(set((a1, a2) for a1, a2, H in all_data_5))
    h_vals_5 = sorted(set(H for a1, a2, H in all_data_5))

    print(f"\n  Total tournaments: {len(all_data_5)}")
    print(f"  Distinct (a1,a2): {len(pairs_5)}")
    print(f"  Distinct H: {len(h_vals_5)}")
    print(f"  H values: {h_vals_5}")

    # Gaps
    max_H_5 = max(h_vals_5)
    gaps_5 = [h for h in range(1, max_H_5+1, 2) if h not in set(h_vals_5)]
    print(f"  Gaps (odd values not achieved): {gaps_5}")

    # Full (a1, a2) -> H table
    print(f"\n  (a1, a2) -> H:")
    for a1, a2 in pairs_5:
        H_val = 1 + 2*a1 + 4*a2
        count = sum(1 for d in all_data_5 if d[0]==a1 and d[1]==a2)
        print(f"    a1={a1:2d}, a2={a2}: H={H_val:3d} ({count} tournaments)")

    # ============================================================
    # PART 2: Exhaustive at n=6
    # ============================================================
    print("\n" + "=" * 70)
    print("PART 2: Exhaustive (a1, a2, H) catalog at n=6")
    print("=" * 70)

    n = 6
    total_bits = n*(n-1)//2
    print(f"  Total tournaments: {2**total_bits}")

    all_data_6 = []
    for bits in range(2**total_bits):
        A = bits_to_adj_list(bits, n)
        H = count_ham_paths(A, n)
        a1, a2 = get_a1_a2(A, n)
        all_data_6.append((a1, a2, H))
        if (bits + 1) % 10000 == 0:
            print(f"    {bits+1}/{2**total_bits}")

    pairs_6 = sorted(set((a1, a2) for a1, a2, H in all_data_6))
    h_vals_6 = sorted(set(H for a1, a2, H in all_data_6))

    print(f"\n  Distinct (a1,a2): {len(pairs_6)}")
    print(f"  Distinct H: {len(h_vals_6)}")
    print(f"  H range: [{min(h_vals_6)}, {max(h_vals_6)}]")

    max_H_6 = max(h_vals_6)
    gaps_6 = [h for h in range(1, max_H_6+1, 2) if h not in set(h_vals_6)]
    print(f"  Gaps: {gaps_6}")

    # (a1, a2) table
    print(f"\n  (a1, a2) -> H:")
    for a1, a2 in pairs_6:
        H_val = 1 + 2*a1 + 4*a2
        count = sum(1 for d in all_data_6 if d[0]==a1 and d[1]==a2)
        print(f"    a1={a1:2d}, a2={a2}: H={H_val:3d} ({count:5d} tournaments)")

    # ============================================================
    # PART 3: Sampled at n=7
    # ============================================================
    print("\n" + "=" * 70)
    print("PART 3: Sampled (a1, a2) at n=7 (5000 samples)")
    print("=" * 70)

    n = 7
    rng = np.random.default_rng(42)
    num_samples = 5000

    all_data_7 = []
    for i in range(num_samples):
        A = [[0]*n for _ in range(n)]
        for ii in range(n):
            for jj in range(ii+1, n):
                if rng.random() < 0.5:
                    A[ii][jj] = 1
                else:
                    A[jj][ii] = 1
        H = count_ham_paths(A, n)
        cycles = []
        for k in range(3, n+1, 2):
            for subset in combinations(range(n), k):
                verts = list(subset)
                nc = count_ham_cycles_exact(A, verts)
                for _ in range(nc):
                    cycles.append(frozenset(subset))
        a1 = len(cycles)
        a2 = 0
        for ii in range(len(cycles)):
            for jj in range(ii+1, len(cycles)):
                if not (cycles[ii] & cycles[jj]):
                    a2 += 1
        all_data_7.append((a1, a2, H))
        if (i+1) % 1000 == 0:
            print(f"    {i+1}/{num_samples}")

    pairs_7 = sorted(set((a1, a2) for a1, a2, H in all_data_7))
    h_vals_7 = sorted(set(H for a1, a2, H in all_data_7))
    h_set_7 = set(h_vals_7)

    print(f"\n  Distinct (a1,a2): {len(pairs_7)}")
    print(f"  Distinct H: {len(h_vals_7)}")
    print(f"  H range: [{min(h_vals_7)}, {max(h_vals_7)}]")

    max_H_7 = max(h_vals_7)
    gaps_7 = [h for h in range(1, max_H_7+1, 2) if h not in h_set_7]
    print(f"  Gaps: {gaps_7}")
    print(f"  Number of gaps: {len(gaps_7)}")

    # ============================================================
    # PART 4: Understanding the gaps
    # ============================================================
    print("\n" + "=" * 70)
    print("PART 4: Understanding H gaps at n=7")
    print("=" * 70)

    for h_gap in gaps_7:
        # For H = h_gap, what (a1,a2) would be needed?
        # H = 1 + 2*a1 + 4*a2
        # a1 = (h_gap - 1 - 4*a2) / 2
        print(f"\n  H = {h_gap}:")
        possible = []
        for a2 in range(0, (h_gap-1)//4 + 1):
            rem = h_gap - 1 - 4*a2
            if rem >= 0 and rem % 2 == 0:
                a1 = rem // 2
                possible.append((a1, a2))

        for a1, a2 in possible:
            seen = (a1, a2) in set(pairs_7)
            # Check if this a1 value appears with any a2
            a1_seen_any = any(p[0] == a1 for p in pairs_7)
            a2_seen_any = any(p[1] == a2 for p in pairs_7)
            # Closest pairs observed
            closest = sorted(pairs_7, key=lambda p: abs(p[0]-a1) + abs(p[1]-a2))[:3]
            print(f"    (a1={a1:2d}, a2={a2:2d}): "
                  f"seen={seen}, a1_exists={a1_seen_any}, a2_exists={a2_seen_any}, "
                  f"closest={closest}")

    # ============================================================
    # PART 5: The a1 achievability at n=7
    # ============================================================
    print("\n" + "=" * 70)
    print("PART 5: Achievable a1 values at n=7")
    print("=" * 70)

    a1_achieved = sorted(set(a1 for a1, a2 in pairs_7))
    a1_gaps = [a for a in range(max(a1_achieved)+1) if a not in set(a1_achieved)]
    print(f"  Achievable a1: {a1_achieved}")
    print(f"  Missing a1: {a1_gaps}")
    print(f"  Note: a1 = total directed odd cycles in T (with multiplicity)")

    # a1 mod 7
    print(f"\n  a1 mod 7: {sorted(set(a % 7 for a in a1_achieved))}")

    # a2 achievability
    a2_achieved = sorted(set(a2 for a1, a2 in pairs_7))
    print(f"\n  Achievable a2: {a2_achieved}")

    # a1 by a2
    print(f"\n  a1 range for each a2:")
    for a2 in a2_achieved:
        a1s = sorted(set(a1 for a1_, a2_ in pairs_7 if a2_ == a2))
        a1_gaps_for_a2 = [a for a in range(min(a1s), max(a1s)+1) if a not in set(a1s)]
        print(f"    a2={a2:2d}: a1 in [{min(a1s)}, {max(a1s)}], "
              f"{len(a1s)} values, gaps={a1_gaps_for_a2[:10]}")

    # ============================================================
    # PART 6: H=7 impossibility analysis
    # ============================================================
    print("\n" + "=" * 70)
    print("PART 6: H=7 impossibility -- structural constraint")
    print("=" * 70)

    print("""
  H=7 requires (a1=3, a2=0) -- the ONLY decomposition.

  a1=3 means exactly 3 directed odd cycles (with multiplicity).
  a2=0 means ALL pairs of cycles share a vertex (all conflicting).

  At n<=6: 3 directed cycles on <=6 vertices always share a vertex
  (by pigeonhole: each cycle uses >=3 vertices, 3*3=9 > 6).
  So a2=0 is automatic from a1=3 at n<=6.

  But we know H=7 is impossible for ALL n (THM-029).
  The proof: a1=3 with a2=0 forces a common vertex, which forces
  a 5-cycle, which gives a1>=4, contradiction.

  The key is that the conflict graph structure (all 3 cycles mutually
  conflicting) geometrically forces additional cycles.
""")

    # Verify at n=5,6: is a1=3 achievable?
    a1_3_n5 = [d for d in all_data_5 if d[0] == 3]
    a1_3_n6 = [d for d in all_data_6 if d[0] == 3]
    print(f"  a1=3 at n=5: {len(a1_3_n5)} tournaments, H values = {sorted(set(d[2] for d in a1_3_n5))}")
    print(f"  a1=3 at n=6: {len(a1_3_n6)} tournaments, H values = {sorted(set(d[2] for d in a1_3_n6))}")

    # a1=3 with a2=0?
    a1_3_a2_0_n5 = [d for d in all_data_5 if d[0] == 3 and d[1] == 0]
    a1_3_a2_0_n6 = [d for d in all_data_6 if d[0] == 3 and d[1] == 0]
    print(f"  a1=3, a2=0 at n=5: {len(a1_3_a2_0_n5)} tournaments")
    print(f"  a1=3, a2=0 at n=6: {len(a1_3_a2_0_n6)} tournaments")

    # ============================================================
    # PART 7: H=21 impossibility analysis
    # ============================================================
    print("\n" + "=" * 70)
    print("PART 7: H=21 impossibility -- deeper constraint")
    print("=" * 70)

    print("""
  H=21 requires (a1, a2) with 1 + 2*a1 + 4*a2 = 21.
  Possible decompositions:
    a2=0: a1=10
    a2=1: a1=8
    a2=2: a1=6
    a2=3: a1=4
    a2=4: a1=2
    a2=5: a1=0 (no, a1=0 means transitive, H=1)

  Each of these requires specific cycle/independence structure.
""")

    # Check each at n=6 (exhaustive)
    for a2_try in range(5):
        a1_try = (20 - 4*a2_try) // 2
        count = sum(1 for d in all_data_6 if d[0]==a1_try and d[1]==a2_try)
        print(f"  n=6: (a1={a1_try}, a2={a2_try}): {count} tournaments")

    # Check at n=7 (sampled)
    for a2_try in range(5):
        a1_try = (20 - 4*a2_try) // 2
        count = sum(1 for d in all_data_7 if d[0]==a1_try and d[1]==a2_try)
        print(f"  n=7: (a1={a1_try}, a2={a2_try}): {count} tournaments (sampled)")

    # Does H=21 appear at n=6?
    h21_n6 = [d for d in all_data_6 if d[2] == 21]
    print(f"\n  H=21 at n=6: {len(h21_n6)} tournaments")
    if h21_n6:
        print(f"    (a1,a2) values: {sorted(set((d[0],d[1]) for d in h21_n6))}")

    # Does H=21 appear at n=5?
    h21_n5 = [d for d in all_data_5 if d[2] == 21]
    print(f"  H=21 at n=5: {len(h21_n5)} tournaments")

    # ============================================================
    # PART 8: Gap persistence across n
    # ============================================================
    print("\n" + "=" * 70)
    print("PART 8: Gap persistence across n")
    print("=" * 70)

    # Which gaps at n=7 also appear at n=5 and n=6?
    h_set_5 = set(h_vals_5)
    h_set_6 = set(h_vals_6)

    print(f"\n  {'H':>5} {'n=5':>5} {'n=6':>5} {'n=7':>5}")
    for h in sorted(set(gaps_7 + gaps_5 + gaps_6)):
        if h <= 20:
            in5 = "GAP" if h not in h_set_5 and h <= max(h_vals_5) else ("yes" if h in h_set_5 else "n/a")
            in6 = "GAP" if h not in h_set_6 and h <= max(h_vals_6) else ("yes" if h in h_set_6 else "n/a")
            in7 = "GAP" if h not in h_set_7 and h <= max(h_vals_7) else ("yes" if h in h_set_7 else "n/a")
            print(f"  {h:>5} {in5:>5} {in6:>5} {in7:>5}")

    print("\nDone.")


if __name__ == '__main__':
    main()
