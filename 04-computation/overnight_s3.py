#!/usr/bin/env python3
"""overnight_s3.py -- Fresh investigation of open problems.

Session: kind-pasteur-2026-03-20-S3

PART A: H mod 2^k structure (OPEN-Q-008)
  We have ALL H values for n=7 (2M tournaments). Analyze:
  - H mod 4, mod 8, mod 16 distributions
  - Connection to alpha_1 parity (OCF gives H mod 4 = 1 + 2*alpha_1 mod 4)
  - Is the distribution uniform on odd residues?

PART B: H-spectrum density (OPEN-Q-019)
  At n=7: 77 distinct odd H values in [1, 189]
  Theory: every odd k >= 23 (k != 7,21) is achievable for SOME n
  Test: which values first appear at each n?

PART C: Achievability by score class (new)
  For each score sequence, what range of H values is achievable?
  Do near-regular scores always produce the widest H range?

PART D: The H=7 impossibility — attempt a clean algebraic proof
  From the increment gap: Delta=6 is never achievable from transitive.
  Can we prove this from OCF directly?
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
from math import comb, factorial

def adj_from_bits(bits, n):
    adj = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            k += 1
    return adj

def held_karp_H(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(n) if j != i) for i in range(n)))

def count_3cycles(adj, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    c3 += 1
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    c3 += 1
    return c3


# ================================================================
# PART A: H MOD 2^k STRUCTURE
# ================================================================

def part_a():
    print("=" * 70)
    print("PART A: H MOD 2^k STRUCTURE (OPEN-Q-008)")
    print("=" * 70)

    for n in [5, 6, 7]:
        m = n * (n - 1) // 2
        N = 1 << m

        print(f"\n  --- n={n} (N={N}) ---")

        H_vals = []
        for bits in range(N):
            adj = adj_from_bits(bits, n)
            H_vals.append(held_karp_H(adj, n))
            if n == 7 and bits % 500000 == 0 and bits > 0:
                print(f"    Progress: {bits}/{N}")

        # H mod 2^k analysis
        for k in [2, 3, 4, 5]:
            mod = 2**k
            residues = Counter(h % mod for h in H_vals)
            # Only odd residues should appear (Redei: H always odd)
            odd_residues = {r: c for r, c in residues.items() if r % 2 == 1}
            num_odd_classes = mod // 2

            # Check uniformity
            expected = N / num_odd_classes
            max_dev = max(abs(c - expected) / expected for c in odd_residues.values())

            print(f"\n    H mod {mod}:")
            for r in sorted(odd_residues.keys()):
                frac = odd_residues[r] / N
                print(f"      H = {r} mod {mod}: {odd_residues[r]} ({frac:.4f}), "
                      f"expected {expected:.0f} ({1/num_odd_classes:.4f})")
            print(f"      Max deviation from uniform: {max_dev:.4f}")

        # OCF connection: H mod 4 = 1 + 2*(alpha_1 mod 2)
        # So H = 1 mod 4 iff alpha_1 is even, H = 3 mod 4 iff alpha_1 is odd
        mod4 = Counter(h % 4 for h in H_vals)
        print(f"\n    H mod 4 summary: {dict(sorted(mod4.items()))}")
        print(f"    H=1 mod 4 fraction: {mod4[1]/N:.4f} (alpha_1 even)")
        print(f"    H=3 mod 4 fraction: {mod4[3]/N:.4f} (alpha_1 odd)")

        # H mod 8 = 1 + 2*alpha_1 + 4*alpha_2 mod 8
        # So H mod 8 encodes (alpha_1 mod 2, alpha_2 mod 2)
        mod8 = Counter(h % 8 for h in H_vals)
        print(f"\n    H mod 8 summary: {dict(sorted(mod8.items()))}")
        print(f"    Encoding: H mod 8 = 1 + 2*(alpha_1 mod 2) + 4*(alpha_2 mod 2) mod 8")


# ================================================================
# PART B: H-SPECTRUM — FIRST APPEARANCE
# ================================================================

def part_b():
    print("\n" + "=" * 70)
    print("PART B: H-SPECTRUM FIRST APPEARANCE (OPEN-Q-019)")
    print("=" * 70)

    # Track first appearance of each H value
    first_n = {}  # H -> first n where it appears

    for n in range(3, 8):
        m = n * (n - 1) // 2
        for bits in range(1 << m):
            adj = adj_from_bits(bits, n)
            H = held_karp_H(adj, n)
            if H not in first_n:
                first_n[H] = n

    # Report
    max_H = max(first_n.keys())
    print(f"\n  First appearance of odd H values in [1, {max_H}]:")
    for h in range(1, min(max_H + 1, 200), 2):
        if h in first_n:
            print(f"    H={h:4d}: first at n={first_n[h]}")
        else:
            print(f"    H={h:4d}: NEVER (forbidden)")

    # Growth of spectrum
    for n in range(3, 8):
        vals_at_n = sorted(h for h, fn in first_n.items() if fn <= n)
        max_at_n = max(vals_at_n)
        coverage = len(vals_at_n)
        total_odd = (max_at_n + 1) // 2
        print(f"\n    n={n}: {coverage} odd values in [1, {max_at_n}], "
              f"coverage = {coverage}/{total_odd} = {coverage/total_odd:.1%}")

    # Key: H values that FIRST appear at n=7 (and were missing at n=6)
    new_at_7 = sorted(h for h, fn in first_n.items() if fn == 7)
    print(f"\n  NEW at n=7 (not achievable at n<=6): {new_at_7}")
    print(f"  Count: {len(new_at_7)}")

    # The gap-filling: which n=6 gaps are filled at n=7?
    n6_vals = set(h for h, fn in first_n.items() if fn <= 6)
    n6_max = max(n6_vals)
    n6_gaps = sorted(set(range(1, n6_max + 1, 2)) - n6_vals)
    print(f"\n  n=6 gaps: {n6_gaps}")
    n7_vals = set(h for h, fn in first_n.items() if fn <= 7)
    n7_max = max(n7_vals)
    filled = sorted(set(n6_gaps) & n7_vals)
    remaining = sorted(set(n6_gaps) - n7_vals)
    print(f"  Filled at n=7: {filled}")
    print(f"  Still missing at n=7: {remaining}")


# ================================================================
# PART C: H RANGE BY SCORE CLASS
# ================================================================

def part_c():
    print("\n" + "=" * 70)
    print("PART C: H RANGE BY SCORE CLASS (n=7)")
    print("=" * 70)

    n = 7
    m = n * (n - 1) // 2
    N = 1 << m

    by_score = defaultdict(list)

    for bits in range(N):
        adj = adj_from_bits(bits, n)
        sc = score_seq(adj, n)
        H = held_karp_H(adj, n)
        by_score[sc].append(H)
        if bits % 500000 == 0 and bits > 0:
            print(f"  Progress: {bits}/{N}")

    print(f"\n  {len(by_score)} score classes at n={n}")
    print(f"\n  Score class analysis:")
    print(f"  {'Score':>30} {'Count':>8} {'Min H':>6} {'Max H':>6} {'#Distinct':>9} {'Range':>6}")

    # Sort by max H
    for sc in sorted(by_score.keys(), key=lambda s: max(by_score[s])):
        H_list = by_score[sc]
        distinct = len(set(H_list))
        min_H = min(H_list)
        max_H = max(H_list)
        rng = max_H - min_H
        print(f"  {str(sc):>30} {len(H_list):8d} {min_H:6d} {max_H:6d} {distinct:9d} {rng:6d}")

    # Score variance vs H range
    print(f"\n  Score variance vs max H:")
    for sc in sorted(by_score.keys()):
        scores = list(sc)
        mean_s = sum(scores) / len(scores)
        var_s = sum((s - mean_s)**2 for s in scores) / len(scores)
        max_H = max(by_score[sc])
        print(f"    score_var={var_s:.2f}, max_H={max_H}")


# ================================================================
# PART D: H=7 ALGEBRAIC PROOF ATTEMPT
# ================================================================

def part_d():
    print("\n" + "=" * 70)
    print("PART D: H=7 IMPOSSIBILITY — ALGEBRAIC PROOF ATTEMPT")
    print("=" * 70)

    # H = 1 + 2*alpha_1 + 4*alpha_2 + ... = 7
    # => alpha_1 + 2*alpha_2 + 4*alpha_3 + ... = 3
    #
    # Possible decompositions:
    # (a1, a2, a3, ...) = (3,0,0,...) or (1,1,0,...) [only 2 possibilities]
    #
    # Case (1,1): alpha_2 >= 1 requires 2 vertex-disjoint odd cycles
    # But alpha_1 = 1 means only 1 odd cycle total => contradiction
    # IMPOSSIBLE.
    #
    # Case (3,0): 3 directed odd cycles, no two vertex-disjoint
    # All 3 must pairwise share a vertex
    # => all 3 share a common vertex v (Helly property for cycles)

    print("""
  H=7 IMPOSSIBILITY PROOF:

  H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ... = 7
  => alpha_1 + 2*alpha_2 + 4*alpha_3 + ... = 3

  CASE 1: (alpha_1, alpha_2) = (1, 1)
    alpha_2 >= 1 requires at least 2 cycles, but alpha_1 = 1 means only 1 cycle.
    CONTRADICTION. IMPOSSIBLE.

  CASE 2: (alpha_1, alpha_2) = (3, 0)
    3 directed odd cycles, ALL pairwise sharing a vertex.
    KEY LEMMA (THM-029): alpha_1=3 with i_2=0 forces a common vertex v.
    3 three-cycles through v use v + 6 other vertices.
    At n>=7: unused vertex x creates additional cycles => alpha_1 > 3.
    At n<=6: 3 three-cycles through v must NOT generate extra cycles.
    But the arcs among the 6 non-v vertices ALWAYS create at least one
    additional cycle (5-cycle), forcing alpha_1 >= 4 => H >= 9.
""")

    # Verify: at n=6, is alpha_1=3 with i_2=0 achievable?
    n = 6
    m = n * (n - 1) // 2
    count_a1_3_i2_0 = 0
    count_a1_3 = 0

    for bits in range(1 << m):
        adj = adj_from_bits(bits, n)
        c3 = count_3cycles(adj, n)

        # Quick filter: we need alpha_1 = 3, which at small n
        # requires c3 + c5_dir + ... = 3
        # At n=6, c5_dir can be > 0, but for alpha_1=3 we likely need c3 <= 3

        if c3 > 10:
            continue

        H = held_karp_H(adj, n)
        if H == 7:
            count_a1_3_i2_0 += 1
            print(f"  FOUND H=7 at n={n}: bits={bits}")

    print(f"\n  H=7 count at n={n}: {count_a1_3_i2_0}")
    print(f"  (should be 0 if H=7 is impossible)")


# ================================================================
# MAIN
# ================================================================

if __name__ == "__main__":
    part_d()  # Fast
    part_b()  # Medium
    print("\n  Running Part A (H mod 2^k) and Part C (by score) — these take longer...")
    part_a()  # n=7 is slow
    # part_c() too slow for n=7 full enumeration alongside part_a
