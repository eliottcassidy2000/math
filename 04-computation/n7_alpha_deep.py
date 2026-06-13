#!/usr/bin/env python3
"""n7_alpha_deep.py -- Deep analysis of alpha_1, alpha_2 tradeoff at n=7 regular.

Session: kind-pasteur-2026-03-20-S1

At n=7 regular, ALL tournaments have c3=14 and alpha_3=0.
The 3 isomorphism classes have:
  H=189: alpha_2=7 (BIBD, minimum!), alpha_1=80
  H=175: alpha_2=14 (maximum!), alpha_1=59
  H=171: alpha_2=10, alpha_1=65

Key question: WHY does the minimizer of alpha_2 maximize alpha_1?

The constraint is: alpha_1 + 2*alpha_2 = S, and H = 1 + 2*S.
  H=189: S = 80 + 14 = 94
  H=175: S = 59 + 28 = 87
  H=171: S = 65 + 20 = 85

So the BIBD class maximizes S=94, giving max H.

In this script, we:
1. Verify these values by computing alpha for specific n=7 regular tournaments
2. Understand the cycle length distribution for each class
3. Find what makes BIBD tournaments have more total cycles
"""

from itertools import combinations, permutations
from collections import Counter
from math import comb

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

def count_ham_cycles_on_subset(adj, verts):
    """Count directed Hamiltonian cycles on vertex subset using Held-Karp."""
    k = len(verts)
    if k < 3:
        return 0
    sub = [[adj[verts[i]][verts[j]] for j in range(k)] for i in range(k)]
    dp = [[0]*k for _ in range(1 << k)]
    dp[1][0] = 1
    for mask in range(1, 1 << k):
        for v in range(k):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(k):
                if mask & (1 << u):
                    continue
                if sub[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << k) - 1
    cycles = 0
    for v in range(1, k):
        if dp[full][v] > 0 and sub[v][0]:
            cycles += dp[full][v]
    return cycles

def full_cycle_analysis(adj, n):
    """Complete directed cycle count by length."""
    by_length = Counter()
    by_vset = {}
    for size in range(3, n+1, 2):
        for verts in combinations(range(n), size):
            c = count_ham_cycles_on_subset(adj, verts)
            if c > 0:
                by_length[size] += c
                by_vset[frozenset(verts)] = c
    return by_length, by_vset

def compute_alpha2(by_vset):
    """Compute alpha_2 from vertex set data."""
    vsets = list(by_vset.keys())
    counts = [by_vset[s] for s in vsets]
    alpha_2 = 0
    for i in range(len(vsets)):
        for j in range(i+1, len(vsets)):
            if not (vsets[i] & vsets[j]):
                alpha_2 += counts[i] * counts[j]
    return alpha_2


def main():
    import random
    random.seed(42)

    n = 7
    m = n * (n - 1) // 2  # 21

    print("=" * 70)
    print("n=7 REGULAR TOURNAMENT ALPHA ANALYSIS")
    print("=" * 70)

    # Find representatives of each H class
    H_representatives = {}  # H -> bits
    attempts = 0
    while len(H_representatives) < 3 and attempts < 100000:
        bits = random.randint(0, (1 << m) - 1)
        adj = adj_from_bits(bits, n)
        sc = score_seq(adj, n)
        if sc != (3, 3, 3, 3, 3, 3, 3):
            attempts += 1
            continue
        H = held_karp_H(adj, n)
        if H not in H_representatives:
            H_representatives[H] = bits
            print(f"  Found H={H} representative (bits={bits}) in {attempts} attempts")
        attempts += 1

    # Analyze each representative
    for H_val in sorted(H_representatives.keys()):
        bits = H_representatives[H_val]
        adj = adj_from_bits(bits, n)

        print(f"\n  --- H={H_val} ---")
        by_length, by_vset = full_cycle_analysis(adj, n)

        alpha_1 = sum(by_length.values())
        alpha_2 = compute_alpha2(by_vset)

        print(f"  Directed cycles by length: {dict(sorted(by_length.items()))}")
        print(f"  Total directed cycles (alpha_1): {alpha_1}")
        print(f"  Disjoint pairs (alpha_2): {alpha_2}")
        print(f"  S = alpha_1 + 2*alpha_2 = {alpha_1 + 2*alpha_2}")
        print(f"  H = 1 + 2*S = {1 + 2*(alpha_1 + 2*alpha_2)}")
        print(f"  Verification: {1 + 2*(alpha_1 + 2*alpha_2) == H_val}")

        # Cycle density per vertex set size
        for size in sorted(by_length.keys()):
            num_vsets_with_cycles = sum(1 for v in by_vset if len(v) == size)
            total_vsets = comb(n, size)
            cycles_this_size = by_length[size]
            avg_per_set = cycles_this_size / max(num_vsets_with_cycles, 1)
            frac_cyclic = num_vsets_with_cycles / total_vsets
            print(f"    {size}-cycles: {cycles_this_size} on {num_vsets_with_cycles}/{total_vsets} sets "
                  f"({frac_cyclic:.1%}), avg {avg_per_set:.2f}/set")

    # KEY ANALYSIS: Why does BIBD (H=189) have more 5-cycles and 7-cycles?
    print(f"\n\n  ===== KEY QUESTION: WHY MORE CYCLES FOR BIBD? =====")
    print(f"  At n=7, the 3 classes all have c3=14 (same # 3-cycle vertex sets).")
    print(f"  The difference is in 5-cycles and 7-cycles.")
    print(f"  BIBD has MORE 5-cycles per 5-vertex subset because its")
    print(f"  3-cycle distribution is more UNIFORM across vertex pairs.")
    print(f"  This creates more regular sub-tournaments on 5 vertices,")
    print(f"  which have more Hamiltonian cycles (H_5 is maximized by regularity).")

    # Verify: check 5-vertex subtournament regularity
    print(f"\n  5-vertex subtournament analysis:")
    for H_val in sorted(H_representatives.keys()):
        bits = H_representatives[H_val]
        adj = adj_from_bits(bits, n)

        sub_H_vals = Counter()
        sub_score_vals = Counter()
        for verts in combinations(range(n), 5):
            sub_adj = [[adj[verts[i]][verts[j]] for j in range(5)] for i in range(5)]
            sub_H = held_karp_H(sub_adj, 5)
            sub_score = score_seq(sub_adj, 5)
            sub_H_vals[sub_H] += 1
            sub_score_vals[sub_score] += 1

        print(f"\n    H={H_val}: 5-vertex sub-tournament H distribution: {dict(sorted(sub_H_vals.items()))}")
        print(f"    Score distribution: {dict(sorted(sub_score_vals.items()))}")

        # Average H of 5-vertex subtournaments
        avg_H5 = sum(h*c for h,c in sub_H_vals.items()) / sum(sub_H_vals.values())
        print(f"    Average H_5 = {avg_H5:.2f}")


if __name__ == "__main__":
    main()
