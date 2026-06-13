#!/usr/bin/env python3
"""
p23_interval_H.py -- Compute H(Interval_23) to test HYP-480

HYP-480: Cyclic interval maximizes H among all circulant tournaments for p >= 13.
Known: H(Paley_23) = 15,760,206,976,379,349 (opus-S10)
Goal: Compute H(Interval_23) and compare.

Uses Held-Karp DP with fixed start vertex (circulant symmetry: H = p * h_0).

Author: kind-pasteur-2026-03-12-S57
"""

import time
import sys


def ham_paths_fixed_start(n, adj_sets):
    """Count Hamiltonian paths in tournament starting at vertex 0.

    adj_sets[v] = frozenset of vertices v can reach.
    Returns h_0 = number of Ham paths starting at 0.
    H(T) = n * h_0 for vertex-transitive tournaments.

    Uses bitmask DP on vertices {1,...,n-1}. Mask tracks which of
    {1,...,n-1} have been visited. Vertex 0 is always visited.
    """
    N = 1 << (n - 1)  # 2^(n-1) masks

    # dp[mask * n + v] = # paths from 0 visiting {0} U support(mask), ending at v
    # Use flat list for speed. Vertex indices 0..n-1.
    dp = [0] * (N * n)
    dp[0 * n + 0] = 1  # mask=0, vertex=0: one empty path

    # Precompute neighbor lists as tuples for speed
    # adj_out[v] = list of (bit_index, vertex) for neighbors of v in {1,...,n-1}
    adj_out = [[] for _ in range(n)]
    for v in range(n):
        for u in range(1, n):
            if u in adj_sets[v]:
                adj_out[v].append((u - 1, u))  # (bit_index, vertex)
    adj_out = [tuple(lst) for lst in adj_out]

    t0 = time.time()
    last_report = t0

    for mask in range(N):
        base = mask * n

        # Progress reporting every 30 seconds
        now = time.time()
        if now - last_report > 30:
            pct = 100.0 * mask / N
            elapsed = now - t0
            rate = mask / elapsed if elapsed > 0 else 0
            eta = (N - mask) / rate if rate > 0 else 0
            print(f"    Progress: {pct:.1f}% ({mask}/{N}), "
                  f"elapsed={elapsed:.0f}s, ETA={eta:.0f}s", flush=True)
            last_report = now

        for v in range(n):
            c = dp[base + v]
            if c == 0:
                continue
            for bit_idx, u in adj_out[v]:
                if mask & (1 << bit_idx):
                    continue
                new_mask = mask | (1 << bit_idx)
                dp[new_mask * n + u] += c

    full_mask = N - 1
    base = full_mask * n
    total = sum(dp[base + v] for v in range(n))

    elapsed = time.time() - t0
    print(f"    DP completed in {elapsed:.1f}s", flush=True)

    return total


def main():
    print("=" * 70)
    print("HYP-480 TEST: H(Interval_23) vs H(Paley_23)")
    print("=" * 70)

    p = 23
    m = (p - 1) // 2  # = 11

    # Known value
    H_paley_23 = 15_760_206_976_379_349
    print(f"\n  H(Paley_23) = {H_paley_23:,} (known from opus-S10)")

    # Interval tournament: S = {1, 2, ..., 11}
    S_interval = frozenset(range(1, m + 1))
    print(f"  Interval S = {{{', '.join(map(str, sorted(S_interval)))}}}")

    # Build adjacency sets
    adj_interval = [frozenset() for _ in range(p)]
    for i in range(p):
        neighbors = frozenset(j for j in range(p) if j != i and (j - i) % p in S_interval)
        adj_interval[i] = neighbors

    print(f"\n  Computing H(Interval_23) via Held-Karp DP...")
    print(f"  State space: 2^{p-1} x {p} = {(1 << (p-1)) * p:,} entries")
    print(f"  Memory: ~{(1 << (p-1)) * p * 8 / 1e9:.1f} GB (int64 flat list)")

    t0 = time.time()
    h0_interval = ham_paths_fixed_start(p, adj_interval)
    H_interval = p * h0_interval
    elapsed = time.time() - t0

    print(f"\n  h_0(Interval_23) = {h0_interval:,}")
    print(f"  H(Interval_23) = {p} * {h0_interval:,} = {H_interval:,}")
    print(f"  Total time: {elapsed:.1f}s")

    # Comparison
    print(f"\n  COMPARISON:")
    print(f"  H(Paley_23)    = {H_paley_23:,}")
    print(f"  H(Interval_23) = {H_interval:,}")
    diff = H_interval - H_paley_23
    pct = 100.0 * diff / H_paley_23 if H_paley_23 > 0 else 0

    if diff > 0:
        print(f"  Interval WINS by {diff:,} ({pct:+.4f}%)")
        print(f"\n  ==> HYP-480 SUPPORTED at p=23: Interval > Paley")
    elif diff < 0:
        print(f"  Paley WINS by {-diff:,} ({pct:+.4f}%)")
        print(f"\n  ==> HYP-480 UNCERTAIN at p=23: Paley still wins")
    else:
        print(f"  TIE!")

    # Also verify Paley
    print(f"\n  Verifying H(Paley_23)...")
    S_paley = frozenset(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1)
    print(f"  Paley QR_23 = {{{', '.join(map(str, sorted(S_paley)))}}}")

    adj_paley = [frozenset() for _ in range(p)]
    for i in range(p):
        neighbors = frozenset(j for j in range(p) if j != i and (j - i) % p in S_paley)
        adj_paley[i] = neighbors

    h0_paley = ham_paths_fixed_start(p, adj_paley)
    H_paley_computed = p * h0_paley
    print(f"  H(Paley_23) computed = {H_paley_computed:,}")
    print(f"  Match with known value: {H_paley_computed == H_paley_23}")

    print(f"\nDONE.")


if __name__ == '__main__':
    main()
