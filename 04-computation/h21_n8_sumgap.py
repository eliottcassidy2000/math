"""
h21_n8_sumgap.py — Investigate whether H=21 occurs at n=8.

H(T) = number of Hamiltonian paths in tournament T.
Uses DP on bitmasks: dp[S][v] = # HP using vertex set S ending at v.
Total work per tournament: O(2^n * n^2) ~ 2^8 * 64 = 16384 ops.

We generate 500k random tournaments on n=8 and check:
1. Does H=21 ever occur?
2. Distribution of H values near 21
3. The "sum gap" alpha_1 + 2*alpha_2 = 10 question
   (H=21 means I(Omega,2)=21, i.e., 1+2*alpha_1+4*alpha_2+...=21,
    so alpha_1+2*alpha_2+... = 10 with higher terms contributing)

Since H = # Hamiltonian paths is much easier to compute than I(Omega,2),
we just check whether H=21 is achievable.
"""

import numpy as np
from collections import Counter
import time
import sys

def count_hp_dp(adj, n):
    """Count Hamiltonian paths using bitmask DP.

    dp[mask][v] = number of Hamiltonian paths using vertices in 'mask' ending at v.
    mask is a bitmask of n bits.
    """
    full = (1 << n) - 1
    # dp[mask][v]: use arrays for speed
    # Initialize: single vertex paths
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1

    for mask in range(1, 1 << n):
        popcount = bin(mask).count('1')
        if popcount < 2:
            continue
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            # v is the last vertex; sum over predecessors u
            prev_mask = mask ^ (1 << v)
            for u in range(n):
                if (prev_mask & (1 << u)) and adj[u][v]:
                    dp[mask, v] += dp[prev_mask, u]

    return int(dp[full].sum())


def count_hp_dp_fast(adj, n):
    """Faster DP using Python lists instead of numpy for small n."""
    full = (1 << n) - 1
    # dp[mask] is a list of n values
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1

    for mask in range(3, 1 << n):  # skip 0 and single-bit masks
        # For each set bit v in mask
        m = mask
        while m:
            v = (m & -m).bit_length() - 1  # lowest set bit
            m &= m - 1
            prev_mask = mask ^ (1 << v)
            if prev_mask == 0:
                continue
            total = 0
            # Sum over predecessors u in prev_mask that have edge u->v
            pm = prev_mask
            while pm:
                u = (pm & -pm).bit_length() - 1
                pm &= pm - 1
                if adj[u][v]:
                    total += dp[prev_mask][u]
            dp[mask][v] = total

    return sum(dp[full])


def main():
    n = 8
    num_samples = 500_000
    batch_size = 1000

    print(f"Investigating H=21 at n={n}")
    print(f"Generating {num_samples} random tournaments")
    print(f"Using bitmask DP: O(2^{n} * {n}^2) = {(1<<n)*n*n} ops per tournament")
    print()

    rng = np.random.default_rng(42)

    h_counter = Counter()
    t0 = time.time()

    # Pre-build adjacency for speed
    for batch_start in range(0, num_samples, batch_size):
        batch_end = min(batch_start + batch_size, num_samples)
        actual_batch = batch_end - batch_start

        for _ in range(actual_batch):
            # Generate random tournament: for each pair i<j, flip a coin
            adj = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1, n):
                    if rng.random() < 0.5:
                        adj[i][j] = 1
                    else:
                        adj[j][i] = 1

            h = count_hp_dp_fast(adj, n)
            h_counter[h] += 1

        elapsed = time.time() - t0
        done = batch_end
        rate = done / elapsed if elapsed > 0 else 0
        eta = (num_samples - done) / rate if rate > 0 else 0

        if done % 10000 == 0 or done == num_samples:
            print(f"  {done}/{num_samples} done  ({rate:.0f}/s, ETA {eta:.0f}s)", flush=True)

    total_time = time.time() - t0
    print(f"\nCompleted in {total_time:.1f}s ({num_samples/total_time:.0f} tournaments/s)")

    # Results
    print("\n" + "="*60)
    print("FULL H DISTRIBUTION (sorted by H value)")
    print("="*60)
    for h_val in sorted(h_counter.keys()):
        count = h_counter[h_val]
        pct = 100.0 * count / num_samples
        print(f"  H={h_val:6d}: {count:8d} ({pct:6.3f}%)")

    print("\n" + "="*60)
    print("H VALUES NEAR 21")
    print("="*60)
    for h_val in range(10, 32):
        count = h_counter.get(h_val, 0)
        marker = " <--- TARGET" if h_val == 21 else ""
        print(f"  H={h_val:4d}: {count:8d}{marker}")

    print("\n" + "="*60)
    print("KEY QUESTION: Does H=21 occur?")
    print("="*60)
    h21_count = h_counter.get(21, 0)
    if h21_count > 0:
        print(f"  YES! H=21 occurred {h21_count} times out of {num_samples}")
    else:
        print(f"  NO. H=21 never occurred in {num_samples} samples.")

    # Check odd values near 21
    print("\n" + "="*60)
    print("ODD H VALUES (H must be odd for odd n by Redei)")
    print("="*60)
    print("  Note: n=8 is EVEN, so H can be even or odd.")
    odd_h = {h: c for h, c in h_counter.items() if h % 2 == 1}
    even_h = {h: c for h, c in h_counter.items() if h % 2 == 0}
    print(f"  Odd H values found:  {len(odd_h)} distinct values, {sum(odd_h.values())} total")
    print(f"  Even H values found: {len(even_h)} distinct values, {sum(even_h.values())} total")

    # Summary statistics
    all_h = []
    for h_val, count in h_counter.items():
        all_h.extend([h_val] * count)
    all_h = np.array(all_h)
    print(f"\n  min(H) = {all_h.min()}")
    print(f"  max(H) = {all_h.max()}")
    print(f"  mean(H) = {all_h.mean():.2f}")
    print(f"  median(H) = {np.median(all_h):.1f}")
    print(f"  std(H) = {all_h.std():.2f}")

    # Check: what H values are MISSING in a range?
    h_min, h_max = int(all_h.min()), int(all_h.max())
    missing = [h for h in range(h_min, h_max+1) if h not in h_counter]
    if missing:
        print(f"\n  Missing H values in [{h_min}, {h_max}]: {missing[:50]}")
        if len(missing) > 50:
            print(f"    ... and {len(missing)-50} more")

    # Check H mod 2 pattern
    print("\n" + "="*60)
    print("H PARITY ANALYSIS")
    print("="*60)
    # At n=8 (even), what's the H parity pattern?
    for parity_name, subset in [("even", even_h), ("odd", odd_h)]:
        if subset:
            vals = sorted(subset.keys())
            print(f"  {parity_name} H values: {vals[:30]}{'...' if len(vals)>30 else ''}")


if __name__ == "__main__":
    main()
