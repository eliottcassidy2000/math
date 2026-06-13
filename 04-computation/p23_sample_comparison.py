#!/usr/bin/env python3
"""
p23_sample_comparison.py -- Sample circulant tournaments at p=23 to check HYP-480

Known: H(Interval_23) = 16,011,537,490,557,279 (max among {Paley, Interval})
Question: Is interval the GLOBAL max among ALL 2048 circulant tournaments?

Sample 30 random tournaments + include the "runner-up patterns" from smaller p.

Author: kind-pasteur-2026-03-12-S57
"""

import time
import random

H_INTERVAL = 16_011_537_490_557_279
H_PALEY = 15_760_206_976_379_349


def ham_paths_fixed_start(n, adj_flat):
    """Count Hamiltonian paths starting at vertex 0.
    adj_flat[v*n + u] = True if v -> u."""
    N = 1 << (n - 1)
    dp = [0] * (N * n)
    dp[0] = 1

    # Precompute adjacency lists
    adj_out = [[] for _ in range(n)]
    for v in range(n):
        for u in range(1, n):
            if adj_flat[v * n + u]:
                adj_out[v].append((u - 1, u))
    adj_out = [tuple(lst) for lst in adj_out]

    for mask in range(N):
        base = mask * n
        for v in range(n):
            c = dp[base + v]
            if c == 0:
                continue
            for bit_idx, u in adj_out[v]:
                if mask & (1 << bit_idx):
                    continue
                dp[(mask | (1 << bit_idx)) * n + u] += c

    full = N - 1
    return sum(dp[full * n + v] for v in range(n))


def build_adj(p, S):
    """Build flat adjacency array for circulant tournament."""
    S_set = set(S)
    adj = [False] * (p * p)
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in S_set:
                adj[i * p + j] = True
    return adj


def all_circulant_sets(p):
    """Generate all valid connection sets for circulant tournament on Z_p."""
    pairs = []
    used = set()
    for a in range(1, p):
        if a not in used:
            b = p - a
            pairs.append((a, b))
            used.add(a)
            used.add(b)
    results = []
    for bits in range(2 ** len(pairs)):
        S = []
        for i, (a, b) in enumerate(pairs):
            S.append(a if (bits >> i) & 1 else b)
        results.append(tuple(sorted(S)))
    return results


def main():
    p = 23
    m = (p - 1) // 2
    print("=" * 70)
    print(f"HYP-480 UNIVERSALITY TEST: Is Interval max at p={p}?")
    print(f"=" * 70)
    print(f"  H(Interval_{p}) = {H_INTERVAL:,}")
    print(f"  H(Paley_{p})    = {H_PALEY:,}")

    all_S = all_circulant_sets(p)
    n_total = len(all_S)
    print(f"  Total circulant tournaments: {n_total}")

    # Identify key tournaments
    S_interval = tuple(sorted(range(1, m + 1)))
    S_paley = tuple(sorted(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1))

    # Sample strategy: random + structured candidates
    random.seed(42)
    sample_indices = random.sample(range(n_total), min(25, n_total))
    sample = [all_S[i] for i in sample_indices]

    # Add tournaments "close to" interval (swap one element)
    # These are the most likely competitors
    S_int_set = set(S_interval)
    for s_in in S_interval:
        s_out = p - s_in  # the element not in S
        new_S = tuple(sorted((S_int_set - {s_in}) | {s_out}))
        if new_S not in sample:
            sample.append(new_S)

    # Remove Paley and Interval from sample (already known)
    sample = [S for S in sample if S != S_interval and S != S_paley]

    print(f"  Sampling {len(sample)} tournaments (25 random + {len(sample) - 25} near-interval)")
    print()

    max_H = H_INTERVAL
    max_S = S_interval
    beats_interval = 0

    for i, S in enumerate(sample):
        adj = build_adj(p, S)
        t0 = time.time()
        h0 = ham_paths_fixed_start(p, adj)
        H = p * h0
        elapsed = time.time() - t0

        marker = ""
        if H > H_INTERVAL:
            marker = " *** BEATS INTERVAL ***"
            beats_interval += 1
            if H > max_H:
                max_H = H
                max_S = S

        pct = 100.0 * (H - H_INTERVAL) / H_INTERVAL
        print(f"  [{i + 1:2d}/{len(sample)}] S={list(S)}: H={H:,} ({pct:+.4f}%) [{elapsed:.1f}s]{marker}",
              flush=True)

    print()
    print("=" * 70)
    print(f"RESULT: {beats_interval} out of {len(sample)} sampled tournaments beat interval")
    if beats_interval == 0:
        print(f"  ==> HYP-480 STRONGLY SUPPORTED at p={p}")
        print(f"  Interval appears to be universal max among circulants")
    else:
        print(f"  ==> HYP-480 REFUTED at p={p}!")
        print(f"  Best: S={list(max_S)}, H={max_H:,}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
