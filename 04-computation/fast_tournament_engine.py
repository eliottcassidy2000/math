#!/usr/bin/env python3
"""
fast_tournament_engine.py — kind-pasteur-2026-03-22-S20m

A FAST tournament computation engine using all our discoveries:
  1. c3 from scores in O(n) [THM proved]
  2. H = C_n - S_2 exact at n<=4 [THM proved]
  3. Gray code traversal with incremental updates
  4. Tr(A^3) for c3 verification
  5. Batch processing of all tournaments at given n

The goal: compute H for ALL tournaments at n=6 (32768) as fast as possible,
collecting statistics that reveal new patterns.

Author: kind-pasteur-2026-03-22-S20m
"""

import sys
import numpy as np
from math import comb, log
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

# ========================================================================
# FAST H COMPUTATION
# ========================================================================

def H_dp(A, n):
    """Held-Karp DP for exact H. O(n^2 * 2^n)."""
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[mask, v] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]:
                    dp[mask | (1 << w), w] += dp[mask, v]
    return int(dp[full].sum())

def bits_to_adj(bits, n, pairs):
    """Convert bit pattern to adjacency matrix."""
    A = np.zeros((n, n), dtype=np.int8)
    for k in range(len(pairs)):
        i, j = pairs[k]
        if (bits >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    return A

def scores_from_bits(bits, n, pairs):
    """Compute scores directly from bits without building full matrix."""
    s = np.zeros(n, dtype=int)
    for k in range(len(pairs)):
        i, j = pairs[k]
        if (bits >> k) & 1:
            s[i] += 1
        else:
            s[j] += 1
    return s

def c3_fast(n, s):
    """c3 from scores in O(n)."""
    S2 = int(sum(s * s))
    return comb(n, 3) - (S2 - comb(n, 2)) // 2

def H_from_scores(n, s):
    """H from scores (exact at n<=4)."""
    S2 = int(sum(s * s))
    return 1 + n*(n-1)*(2*n-1)//6 - S2

# ========================================================================
# EXHAUSTIVE COMPUTATION WITH STATISTICS
# ========================================================================

def exhaustive_analysis(n, max_n_for_H=6):
    """Compute all tournament invariants at order n."""
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    total = 2**m

    print(f"\n{'='*60}")
    print(f"  EXHAUSTIVE ANALYSIS AT n={n} ({total} tournaments)")
    print(f"{'='*60}")

    start = time.time()

    # Collect data
    H_dist = defaultdict(int)
    c3_dist = defaultdict(int)
    score_to_H = defaultdict(set)
    H_to_scores = defaultdict(set)
    delta_dist = defaultdict(int)

    # Angle counts
    N120_to_H = defaultdict(set)

    # The H sequence along bit-order (not Gray code, for speed)
    prev_H = None
    H_sum = 0
    H_sq_sum = 0

    for bits in range(total):
        # Fast score computation
        s = scores_from_bits(bits, n, pairs)
        S2 = int(sum(s * s))
        c3 = comb(n, 3) - (S2 - comb(n, 2)) // 2
        score_key = tuple(sorted(s))

        # H computation
        if n <= 4:
            H = 1 + n*(n-1)*(2*n-1)//6 - S2
        elif n <= max_n_for_H:
            A = bits_to_adj(bits, n, pairs)
            H = H_dp(A, n)
        else:
            H = None

        if H is not None:
            H_dist[H] += 1
            H_sum += H
            H_sq_sum += H * H
            score_to_H[score_key].add(H)
            H_to_scores[H].add(score_key)

            # Angle count
            N120 = n*(n-1)*(n-2)//2 - S2
            N120_to_H[N120].add(H)

        c3_dist[c3] += 1

        if bits % 10000 == 0 and bits > 0:
            elapsed = time.time() - start
            rate = bits / elapsed
            eta = (total - bits) / rate
            print(f"  ... {bits}/{total} ({100*bits/total:.0f}%) "
                  f"rate={rate:.0f}/s ETA={eta:.0f}s")

    elapsed = time.time() - start
    print(f"\n  Completed in {elapsed:.1f}s ({total/elapsed:.0f} tournaments/s)")

    # ================================================================
    # RESULTS
    # ================================================================

    if H_dist:
        H_vals = sorted(H_dist.keys())
        print(f"\n  H DISTRIBUTION:")
        for H in H_vals:
            bar = "#" * min(50, H_dist[H] * 50 // max(H_dist.values()))
            print(f"    H={H:>5d}: {H_dist[H]:>6d}  {bar}")

        # H statistics
        mean_H = H_sum / total
        var_H = H_sq_sum / total - mean_H**2
        print(f"\n  H STATISTICS:")
        print(f"    Mean H = {mean_H:.4f}")
        print(f"    Var H = {var_H:.4f}")
        print(f"    E[H] = n!/2^(n-1) = {int(np.prod(np.arange(1, n+1))) / 2**(n-1):.4f}")
        print(f"    Min H = {min(H_vals)}, Max H = {max(H_vals)}")

        # Forbidden values
        all_odd = set(range(1, max(H_vals)+1, 2))
        gaps = sorted(all_odd - set(H_vals))
        if gaps:
            print(f"\n  FORBIDDEN H VALUES (odd gaps): {gaps[:20]}{'...' if len(gaps)>20 else ''}")

        # Score determination
        multi_H = sum(1 for v in score_to_H.values() if len(v) > 1)
        total_scores = len(score_to_H)
        print(f"\n  SCORE DETERMINATION:")
        print(f"    Score classes: {total_scores}")
        print(f"    Score classes with multiple H: {multi_H}")
        if multi_H > 0:
            for score_key, Hs in sorted(score_to_H.items()):
                if len(Hs) > 1:
                    print(f"      scores={list(score_key)}: H in {sorted(Hs)}")

        # OCR computation
        score_means = {}
        score_counts = {}
        for bits in range(total):
            s = scores_from_bits(bits, n, pairs)
            S2 = int(sum(s * s))
            score_key = tuple(sorted(s))
            if score_key not in score_means:
                score_means[score_key] = 0
                score_counts[score_key] = 0

        # Recompute with H
        if n <= max_n_for_H:
            score_H_sums = defaultdict(float)
            score_H_counts = defaultdict(int)
            for bits in range(min(total, 50000)):  # sample for large n
                s = scores_from_bits(bits, n, pairs)
                score_key = tuple(sorted(s))
                if n <= 4:
                    S2 = int(sum(s * s))
                    H = 1 + n*(n-1)*(2*n-1)//6 - S2
                else:
                    A = bits_to_adj(bits, n, pairs)
                    H = H_dp(A, n)
                score_H_sums[score_key] += H
                score_H_counts[score_key] += 1

            # Var(E[H|score]) / Var(H)
            cond_means = {k: score_H_sums[k]/score_H_counts[k] for k in score_H_sums}
            overall_mean = sum(score_H_sums.values()) / sum(score_H_counts.values())
            var_cond = sum(score_H_counts[k] * (cond_means[k] - overall_mean)**2
                         for k in cond_means) / sum(score_H_counts.values())
            if var_H > 0:
                OCR = var_cond / var_H
                print(f"\n  OCR (score determines H):")
                print(f"    OCR = {OCR:.6f} ({100*OCR:.2f}%)")

    # c3 distribution
    print(f"\n  c3 DISTRIBUTION:")
    for c3_val in sorted(c3_dist.keys()):
        print(f"    c3={c3_val:>3d}: {c3_dist[c3_val]:>6d} tournaments")

    # N120 -> H relationship
    if N120_to_H:
        print(f"\n  N_120 -> H:")
        for N120 in sorted(N120_to_H.keys()):
            Hs = sorted(N120_to_H[N120])
            det = "EXACT" if len(Hs) == 1 else "multiple"
            print(f"    N_120={N120:>4d}: H in {Hs} ({det})")

    # Delta analysis (Gray code)
    if n <= 5 and H_dist:
        print(f"\n  GRAY CODE DELTA ANALYSIS:")
        prev_H = None
        delta_counts = defaultdict(int)

        def gray_gen(m):
            for i in range(2**m):
                yield i ^ (i >> 1)

        for gray in gray_gen(m):
            s = scores_from_bits(gray, n, pairs)
            S2 = int(sum(s * s))
            if n <= 4:
                H = 1 + n*(n-1)*(2*n-1)//6 - S2
            else:
                A = bits_to_adj(gray, n, pairs)
                H = H_dp(A, n)

            if prev_H is not None:
                delta = H - prev_H
                delta_counts[delta] += 1
            prev_H = H

        for d in sorted(delta_counts.keys()):
            print(f"    delta={d:>+4d}: {delta_counts[d]:>5d}")

        # Check for missing deltas
        min_d = min(delta_counts.keys())
        max_d = max(delta_counts.keys())
        all_even = set(range(min_d, max_d+1, 2))
        missing = sorted(all_even - set(delta_counts.keys()))
        if missing:
            print(f"    MISSING deltas: {missing}")

    return H_dist

# ========================================================================
# MAIN
# ========================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("  FAST TOURNAMENT ENGINE")
    print("=" * 60)

    # Run at n=3,4,5 (fast)
    for n in [3, 4, 5]:
        exhaustive_analysis(n)

    # n=6: this is the big one
    print("\n  Starting n=6 exhaustive computation...")
    exhaustive_analysis(6, max_n_for_H=6)
