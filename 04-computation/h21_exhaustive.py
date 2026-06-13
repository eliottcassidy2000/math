#!/usr/bin/env python3
"""
Exhaustive check: does H(T)=21 occur for any tournament at n=3,...,8?

Uses fast Held-Karp DP for H(T) computation.

For n=7: 2^21 = 2,097,152 tournaments (feasible).
For n=8: 2^28 = 268,435,456 tournaments (need bitwise tricks or sampling).

opus-2026-03-07-S38
"""
import sys
from collections import Counter

def held_karp_count(n, adj):
    """Count Hamiltonian paths using Held-Karp DP. O(2^n * n^2)."""
    # dp[S][v] = number of Hamiltonian paths ending at v using vertex set S
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1

    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if adj[v] & (1 << u):  # v -> u
                    dp[S | (1 << u)][u] += dp[S][v]

    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def enumerate_tournaments(n):
    """Enumerate all tournaments on n vertices via bit encoding."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for bits in range(1 << m):
        # Build adjacency bitmask
        adj = [0] * n
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                adj[j] |= (1 << i)  # j -> i
            else:
                adj[i] |= (1 << j)  # i -> j
        yield adj


def count_cycles_fast(n, adj):
    """Count directed odd cycles (3-cycles and 5-cycles for small n)."""
    # Count directed 3-cycles
    t3 = 0
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                # Check a->b->c->a
                if (adj[a] >> b & 1) and (adj[b] >> c & 1) and (adj[c] >> a & 1):
                    t3 += 1
                elif (adj[a] >> c & 1) and (adj[c] >> b & 1) and (adj[b] >> a & 1):
                    t3 += 1
    return t3


def main():
    target = 21
    print(f"=== Searching for H(T) = {target} ===\n")

    for n in range(3, 8):
        edges = [(i, j) for i in range(n) for j in range(i+1, n)]
        m = len(edges)
        total = 1 << m

        if total > 300_000_000:
            print(f"n={n}: {total} tournaments — too many for exhaustive, sampling instead")
            import random
            found = False
            H_counter = Counter()
            num_samples = 2_000_000
            for seed in range(num_samples):
                rng = random.Random(seed)
                adj = [0] * n
                for i in range(n):
                    for j in range(i+1, n):
                        if rng.random() < 0.5:
                            adj[i] |= (1 << j)
                        else:
                            adj[j] |= (1 << i)
                H = held_karp_count(n, adj)
                H_counter[H] += 1
                if H == target:
                    found = True
                    print(f"  FOUND H={target} at n={n}, seed={seed}!")

            if not found:
                print(f"  H={target} NOT found in {num_samples} samples")
            # Print nearby values
            for h in sorted(H_counter.keys()):
                if abs(h - target) <= 4:
                    print(f"    H={h}: {H_counter[h]} occurrences")
            continue

        found = False
        H_counter = Counter()
        for bits in range(total):
            adj = [0] * n
            for k, (i, j) in enumerate(edges):
                if bits & (1 << k):
                    adj[j] |= (1 << i)
                else:
                    adj[i] |= (1 << j)
            H = held_karp_count(n, adj)
            H_counter[H] += 1
            if H == target:
                found = True
                # Only print first few
                if H_counter[target] <= 3:
                    print(f"  FOUND H={target} at n={n}, bits={bits}")

        if found:
            print(f"n={n}: H={target} occurs {H_counter[target]} times out of {total}")
        else:
            print(f"n={n}: H={target} NEVER occurs (out of {total} tournaments)")

        # Print all achievable H values
        achievable = sorted(H_counter.keys())
        print(f"  Achievable: {achievable}")

    # Special focus on near-21 values at n=7
    print(f"\n=== Detailed n=7 analysis: H values near 21 ===")
    n = 7
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    target_range = range(17, 26)  # 17 to 25
    examples = {h: None for h in target_range}
    H_counter = Counter()

    for bits in range(1 << m):
        adj = [0] * n
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                adj[j] |= (1 << i)
            else:
                adj[i] |= (1 << j)
        H = held_karp_count(n, adj)
        H_counter[H] += 1
        if H in target_range and examples[H] is None:
            examples[H] = bits

    print(f"H values in [17,25]:")
    for h in target_range:
        count = H_counter.get(h, 0)
        has_example = examples.get(h) is not None
        print(f"  H={h}: {count} tournaments" +
              (f" (example bits={examples[h]})" if has_example else ""))

    # Print ALL achievable H values at n=7
    print(f"\nAll achievable H values at n=7:")
    for h in sorted(H_counter.keys()):
        print(f"  H={h}: {H_counter[h]}")


if __name__ == "__main__":
    main()
