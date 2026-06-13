#!/usr/bin/env python3
"""
2-adic Tower and THM-J Verification (kind-pasteur-S34).

Tests:
1. THM-J criterion: S mod 2^{n-1} universal iff s_2(n-3) <= 1
2. 2-adic tower: H(T) mod 2^k for k = 1, 2, 3, ...
3. Connection between v_2(H(T)-1) and tournament structure
4. Whether gaps in H-spectrum correspond to 2-adic obstructions

Novel hypothesis: The H-gap at 7 and 21 might be visible 2-adically.
H=7: 7 = 0b111, 7-1 = 6 = 2*3, v_2(6) = 1
H=21: 21 = 0b10101, 21-1 = 20 = 4*5, v_2(20) = 2
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations
from math import comb


def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj


def count_ham_paths(adj, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])


def count_3cycles(adj, n):
    count = 0
    for a, b, c in combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or \
           (adj[a][c] and adj[c][b] and adj[b][a]):
            count += 1
    return count


def signed_hp_permanent(adj, n):
    """Compute S(T) = sum_P prod B[P_i][P_{i+1}], B = 2A-1."""
    B = [[2*adj[i][j] - 1 if i != j else 0 for j in range(n)] for i in range(n)]
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)):
                    dp[mask | (1 << u)][u] += dp[mask][v] * B[v][u]
    return sum(dp[full])


def v2(n):
    """2-adic valuation of n."""
    if n == 0:
        return float('inf')
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def s2(n):
    """Binary digit sum of n."""
    return bin(n).count('1')


def main():
    print("=== 2-adic Tower and THM-J Verification ===\n")

    # Part 1: Verify THM-J criterion
    print("--- Part 1: THM-J Criterion ---")
    print("S mod 2^{n-1} universal iff s_2(n-3) <= 1\n")
    for n in range(3, 40, 2):
        criterion = s2(n - 3) <= 1
        n_minus_3 = n - 3
        print(f"  n={n:2d}: n-3={n_minus_3:2d} = 0b{n_minus_3:06b}, "
              f"s_2(n-3)={s2(n_minus_3)}, universal={'YES' if criterion else 'NO'}")

    # Part 2: H(T) mod 2^k distribution
    print("\n--- Part 2: H(T) mod 2^k Distribution ---")
    for n in [4, 5, 6, 7]:
        print(f"\n  n={n}:")
        num_edges = n * (n - 1) // 2
        max_bits = 1 << num_edges

        h_values = []
        for bits in range(max_bits):
            adj = tournament_adj(n, bits)
            H = count_ham_paths(adj, n)
            h_values.append(H)

        # Distribution of H mod 2^k
        for k in range(1, 6):
            mod = 1 << k
            dist = {}
            for H in h_values:
                r = H % mod
                dist[r] = dist.get(r, 0) + 1
            print(f"    H mod {mod:2d}: {dict(sorted(dist.items()))}")

    # Part 3: v_2(H-1) distribution
    print("\n--- Part 3: v_2(H-1) Distribution ---")
    for n in [5, 6, 7]:
        print(f"\n  n={n}:")
        num_edges = n * (n - 1) // 2
        max_bits = 1 << num_edges

        v2_dist = {}
        for bits in range(max_bits):
            adj = tournament_adj(n, bits)
            H = count_ham_paths(adj, n)
            val = v2(H - 1)
            if val == float('inf'):
                val = 'inf'
            v2_dist[val] = v2_dist.get(val, 0) + 1

        for val in sorted(k for k in v2_dist if isinstance(k, int)):
            print(f"    v_2(H-1) = {val}: {v2_dist[val]} tournaments")
        if 'inf' in v2_dist:
            print(f"    v_2(H-1) = inf (H=1): {v2_dist['inf']} tournaments")

    # Part 4: 2-adic analysis of gap values
    print("\n--- Part 4: 2-adic Analysis of Gap Values ---")
    gaps = [7, 21, 23]  # known or candidate gaps
    for g in gaps:
        print(f"  H={g}: binary={g:b}, v_2(H)={v2(g)}, v_2(H-1)={v2(g-1)}, "
              f"v_2(H+1)={v2(g+1)}")

    # Achievable H values at each n
    for n in [5, 6, 7, 8]:
        num_edges = n * (n - 1) // 2
        max_bits = 1 << num_edges
        achieved = set()
        for bits in range(max_bits):
            adj = tournament_adj(n, bits)
            H = count_ham_paths(adj, n)
            achieved.add(H)

        # Look at gaps
        max_h = max(achieved)
        all_odd = sorted(h for h in range(1, max_h + 1, 2))
        missing_odd = [h for h in all_odd if h not in achieved]
        print(f"\n  n={n}: achieved {len(achieved)} values, max={max_h}")
        print(f"    Missing odd values up to max: {missing_odd[:20]}")

        # 2-adic pattern of missing values
        if missing_odd:
            for h in missing_odd[:10]:
                print(f"      H={h}: v_2(H-1)={v2(h-1)}, "
                      f"H mod 4={h%4}, H mod 8={h%8}")

    # Part 5: S(T) computation at n=5 with t3 correlation
    print("\n--- Part 5: S(T) at n=5 ---")
    n = 5
    num_edges = n * (n - 1) // 2
    max_bits = 1 << num_edges

    s_by_t3 = {}
    for bits in range(max_bits):
        adj = tournament_adj(n, bits)
        H = count_ham_paths(adj, n)
        t3 = count_3cycles(adj, n)
        S = signed_hp_permanent(adj, n)

        if t3 not in s_by_t3:
            s_by_t3[t3] = set()
        s_by_t3[t3].add(S)

    for t3 in sorted(s_by_t3):
        vals = sorted(s_by_t3[t3])
        print(f"  t3={t3}: S values = {vals}, v_2 = {[v2(abs(s)) if s != 0 else 'inf' for s in vals]}")

    # Part 6: S(T) at n=7 (sample)
    print("\n--- Part 6: S(T) at n=7 (exhaustive) ---")
    n = 7
    num_edges = n * (n - 1) // 2
    max_bits = 1 << num_edges

    import random
    random.seed(42)

    s_mod64 = {}
    s_values = set()
    count = 0
    # Exhaustive at n=7 is 2^21 = 2M, feasible but slow. Sample instead.
    for trial in range(10000):
        bits = random.randint(0, max_bits - 1)
        adj = tournament_adj(n, bits)
        S = signed_hp_permanent(adj, n)
        s_values.add(S)
        r = S % 64
        s_mod64[r] = s_mod64.get(r, 0) + 1
        count += 1

    print(f"  {count} samples, {len(s_values)} distinct S values")
    print(f"  S mod 64: {dict(sorted(s_mod64.items()))}")
    print(f"  min(v_2(S)): {min(v2(abs(s)) if s != 0 else float('inf') for s in s_values)}")

    # Check S mod 2^6 = 48 always (THM-H prediction)
    all_48 = all(s % 64 == 48 or s % 64 == -16 % 64 for s in s_values)
    print(f"  All S ≡ 48 mod 64: {all_48}")
    for s in sorted(s_values)[:5]:
        print(f"    S={s}, S mod 64={s%64}, v_2(S)={v2(abs(s)) if s != 0 else 'inf'}")
    for s in sorted(s_values)[-5:]:
        print(f"    S={s}, S mod 64={s%64}, v_2(S)={v2(abs(s)) if s != 0 else 'inf'}")


if __name__ == "__main__":
    main()
