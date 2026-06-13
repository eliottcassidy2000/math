#!/usr/bin/env python3
"""
opus-2026-04-05-S25: Three practical speedups for A000568 computation.

STRATEGY 1: DP Memoization on edge-count
  Group partitions by their total edge count e(λ). Accumulate 1/z_λ
  per edge count. Final answer = Σ_e bucket[e] × 2^e / n!.
  This is what the GMP code already does (bucket accumulation).
  But we can go further: memoize partial states in the partition DP.

STRATEGY 2: Partition Pruning
  The contribution of a partition λ is (n!/z_λ) × 2^{e(λ)}.
  For large partitions dominated by 1s: z_λ ≈ n! (large), e(λ) ≈ 0 (small).
  Prune partitions where 2^{e(λ)} / z_λ < ε / (number of remaining partitions).

  Key insight: the edge count e(λ) grows with the number of large parts.
  Partitions with many 1s have e(λ) ≈ 0 and contribute ≈ 1 to the sum.
  Partitions with few large parts have e(λ) ≈ C(n,2)/2 and z_λ small.

STRATEGY 3: Asymptotic Recursion
  a(n) = a(n-1) × 2^{n-1}/n × (1 - correction)
  The correction comes from the 3-cycle partition type.
  For APPROXIMATE values, this gives log-n digits of precision.
  For EXACT computation, use as a CHECKSUM.

Let's benchmark all three and compare with the GMP reference.
"""

import time
from math import factorial, gcd, log2
from functools import lru_cache
from collections import defaultdict

# ─── Reference values ───
KNOWN = {
    1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
    9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168,
    13: 48542114686912, 14: 28401423719122304, 15: 31021002160355166848,
}

# ═══════════════════════════════════════════════════════
# METHOD 1: Bucket DP with memoization
# ═══════════════════════════════════════════════════════

def a000568_bucket_dp(n):
    """
    Compute a(n) using the standard Burnside formula with bucket accumulation.
    This is the Python equivalent of the GMP C code's approach.
    """
    if n <= 2:
        return 1

    nfact = factorial(n)
    odd_parts = list(range(n if n % 2 else n-1, 0, -2))

    # Precompute gcd table
    gcd_tab = {}
    for a in odd_parts:
        for b in odd_parts:
            gcd_tab[(a,b)] = gcd(a,b)

    # bucket[e] = sum of n!/z_lambda for all partitions with edge count e
    buckets = defaultdict(int)
    count = [0]

    def recurse(idx, remaining, parts_placed, edges):
        """
        idx: current index in odd_parts
        remaining: sum left to allocate
        parts_placed: list of (k, m) pairs placed so far
        edges: accumulated edge count
        """
        if remaining == 0:
            # Compute n!/z_lambda
            z = 1
            for k, m in parts_placed:
                z *= k**m * factorial(m)
            buckets[edges] += nfact // z
            count[0] += 1
            return

        if idx >= len(odd_parts):
            return

        k = odd_parts[idx]
        if k > remaining:
            recurse(idx + 1, remaining, parts_placed, edges)
            return

        # m = 0
        recurse(idx + 1, remaining, parts_placed, edges)

        # m = 1, 2, ...
        for m in range(1, remaining // k + 1):
            # Self-edges: m*(k-1)/2 + C(m,2)*k
            e_self = m * (k - 1) // 2 + m * (m - 1) // 2 * k

            # Cross-edges with previously placed parts
            e_cross = sum(m * pm * gcd_tab[(k, pk)] for pk, pm in parts_placed)

            new_edges = edges + e_self + e_cross
            recurse(idx + 1, remaining - m * k,
                    parts_placed + [(k, m)], new_edges)

    recurse(0, n, [], 0)

    # Final: a(n) = (1/n!) * Σ_e bucket[e] * 2^e
    total = sum(coeff * (2 ** e) for e, coeff in buckets.items())
    return total // nfact

# ═══════════════════════════════════════════════════════
# METHOD 2: Pruned enumeration with bounds
# ═══════════════════════════════════════════════════════

def a000568_pruned(n, prune_threshold=1e-15):
    """
    Like bucket DP but prune partitions with negligible contribution.

    Bound: contribution of partition λ ≤ 2^{e_max} / z_λ
    where e_max = C(n,2)/2 (max possible edge count).

    If 2^{e(λ)} / z_λ is tiny relative to the running sum, skip.
    Actually the tighter bound uses the ACTUAL edge count.

    For exact computation, we DON'T prune (would lose precision).
    For approximate computation, prune aggressively.
    """
    if n <= 2:
        return 1

    nfact = factorial(n)
    max_possible_edges = n * (n - 1) // 4  # approximate max for tournaments
    odd_parts = list(range(n if n % 2 else n-1, 0, -2))

    gcd_tab = {}
    for a in odd_parts:
        for b in odd_parts:
            gcd_tab[(a,b)] = gcd(a,b)

    total = [0]
    pruned_count = [0]
    total_count = [0]

    def recurse(idx, remaining, parts_placed, edges, z_so_far):
        if remaining == 0:
            contrib = (nfact // z_so_far) * (2 ** edges)
            total[0] += contrib
            total_count[0] += 1
            return

        if idx >= len(odd_parts):
            return

        k = odd_parts[idx]
        if k > remaining:
            recurse(idx + 1, remaining, parts_placed, edges, z_so_far)
            return

        # m = 0
        recurse(idx + 1, remaining, parts_placed, edges, z_so_far)

        # Pruning: even with maximum possible remaining edges,
        # can this partition's contribution matter?
        # Remaining parts all being 1 gives edges ≈ 0 extra.
        # Remaining parts being one big part gives most edges.
        max_remaining_edges = remaining * (remaining - 1) // 4

        for m in range(1, remaining // k + 1):
            new_z = z_so_far * k**m * factorial(m)
            e_self = m * (k - 1) // 2 + m * (m - 1) // 2 * k
            e_cross = sum(m * pm * gcd_tab[(k, pk)] for pk, pm in parts_placed)
            new_edges = edges + e_self + e_cross

            # Pruning check: if 2^{new_edges + max_remaining_edges} / new_z < threshold * total
            if total[0] > 0 and prune_threshold > 0:
                max_contrib = (nfact / new_z) * 2**(new_edges + max_remaining_edges)
                if max_contrib < prune_threshold * total[0]:
                    pruned_count[0] += 1
                    continue

            recurse(idx + 1, remaining - m * k,
                    parts_placed + [(k, m)], new_edges, new_z)

    recurse(0, n, [], 0, 1)

    result = total[0] // nfact
    return result, total_count[0], pruned_count[0]

# ═══════════════════════════════════════════════════════
# METHOD 3: Asymptotic recursion
# ═══════════════════════════════════════════════════════

def a000568_asymptotic(n, prev_value=None):
    """
    a(n) ≈ a(n-1) * 2^{n-1} / n * (1 - correction)

    The correction comes from the 3-cycle partition contribution.
    Main term: identity partition (1^n) → 2^{C(n,2)}/n!
    3-cycle correction: partitions with one 3-cycle → C(n,3) * 2^{e}/z

    From S14: correction ≈ (n-1)(n-2)(n-4) / 4^{n-2}
    """
    if n <= 2:
        return 1
    if prev_value is None:
        # Bootstrap from known or recursive
        if n in KNOWN:
            return KNOWN[n]
        prev_value = a000568_asymptotic(n - 1)

    # Main term: a(n-1) * 2^{n-1} / n
    main = prev_value * 2**(n-1) / n

    # 3-cycle correction
    if n >= 5:
        correction = (n-1) * (n-2) * (n-4) / 4**(n-2)
    else:
        correction = 0

    return main * (1 - correction)

def a000568_asymptotic_chain(max_n):
    """Compute the asymptotic approximation for all n up to max_n."""
    results = {1: 1, 2: 1}
    for n in range(3, max_n + 1):
        results[n] = a000568_asymptotic(n, results[n-1])
    return results

# ═══════════════════════════════════════════════════════
# METHOD 4: HYBRID — exact small terms + asymptotic tail
# ═══════════════════════════════════════════════════════

def partition_contribution_sorted(n):
    """
    Enumerate partitions sorted by contribution magnitude (descending).
    The largest contributions come from: identity (1^n), single 3-cycle + 1s,
    two 3-cycles + 1s, single 5-cycle + 1s, etc.

    Returns: list of (edge_count, z_lambda, contribution_fraction) tuples.
    """
    nfact = factorial(n)
    odd_parts = list(range(n if n % 2 else n-1, 0, -2))

    gcd_tab = {}
    for a in odd_parts:
        for b in odd_parts:
            gcd_tab[(a,b)] = gcd(a,b)

    partitions = []

    def recurse(idx, remaining, parts, edges, z):
        if remaining == 0:
            contrib = (nfact // z) * (2 ** edges)
            partitions.append((contrib, edges, z, list(parts)))
            return
        if idx >= len(odd_parts):
            return
        k = odd_parts[idx]
        if k > remaining:
            recurse(idx + 1, remaining, parts, edges, z)
            return
        recurse(idx + 1, remaining, parts, edges, z)
        for m in range(1, remaining // k + 1):
            e_self = m * (k-1)//2 + m*(m-1)//2*k
            e_cross = sum(m * pm * gcd_tab[(k, pk)] for pk, pm in parts)
            recurse(idx + 1, remaining - m*k,
                    parts + [(k, m)], edges + e_self + e_cross,
                    z * k**m * factorial(m))

    recurse(0, n, [], 0, 1)
    partitions.sort(reverse=True)
    return partitions

# ═══════════════════════════════════════════════════════
# BENCHMARKS
# ═══════════════════════════════════════════════════════

if __name__ == "__main__":
    print("="*70)
    print("A000568 SPEEDUP BENCHMARKS")
    print("="*70)

    # Method 1: Bucket DP (exact)
    print("\n── Method 1: Bucket DP (exact) ──")
    for n in range(3, 14):
        t0 = time.time()
        result = a000568_bucket_dp(n)
        elapsed = time.time() - t0
        expected = KNOWN.get(n)
        match = "✓" if result == expected else f"✗ (got {result})"
        print(f"  a({n:2d}) = {result:>20d}  {elapsed:8.4f}s  {match}")
        if elapsed > 10:
            break

    # Method 2: Pruned (exact with threshold=0)
    print("\n── Method 2: Pruned enumeration (threshold=0 = exact) ──")
    for n in [10, 12, 14]:
        if n > 14: break
        t0 = time.time()
        result, total, pruned = a000568_pruned(n, prune_threshold=0)
        elapsed = time.time() - t0
        match = "✓" if result == KNOWN.get(n) else "✗"
        print(f"  a({n:2d}) = {result:>20d}  {elapsed:8.4f}s  partitions={total}  {match}")

    # Method 2b: Aggressive pruning (approximate)
    print("\n── Method 2b: Aggressive pruning (approximate) ──")
    for n in [15, 20, 25]:
        t0 = time.time()
        result, total, pruned = a000568_pruned(n, prune_threshold=1e-10)
        elapsed = time.time() - t0
        if n in KNOWN:
            rel_err = abs(result - KNOWN[n]) / KNOWN[n]
            print(f"  a({n:2d}) ≈ {result:>20d}  {elapsed:8.4f}s  "
                  f"counted={total} pruned={pruned}  err={rel_err:.2e}")
        else:
            print(f"  a({n:2d}) ≈ {result:>20d}  {elapsed:8.4f}s  "
                  f"counted={total} pruned={pruned}")

    # Method 3: Asymptotic recursion
    print("\n── Method 3: Asymptotic recursion ──")
    approx = a000568_asymptotic_chain(20)
    for n in range(3, 16):
        exact = KNOWN.get(n)
        if exact:
            rel_err = abs(approx[n] - exact) / exact
            print(f"  a({n:2d}): approx={approx[n]:>20.1f}  exact={exact:>20d}  err={rel_err:.6f}")

    # Method 4: Partition contribution ranking
    print("\n── Method 4: Partition contribution ranking (n=10) ──")
    parts = partition_contribution_sorted(10)
    total = sum(c for c, _, _, _ in parts)
    cumulative = 0
    print(f"  Total partitions: {len(parts)}, total sum = {total}")
    for i, (contrib, e, z, p) in enumerate(parts[:10]):
        cumulative += contrib
        frac = cumulative / total
        print(f"  #{i+1}: partition={p}, e={e}, contrib/total={contrib/total:.6f}, cumulative={frac:.6f}")

    # How many partitions needed for 99%, 99.9%, 99.99% of the sum?
    cumulative = 0
    thresholds = [0.99, 0.999, 0.9999, 0.99999]
    threshold_idx = 0
    for i, (contrib, e, z, p) in enumerate(parts):
        cumulative += contrib
        while threshold_idx < len(thresholds) and cumulative / total >= thresholds[threshold_idx]:
            print(f"  {thresholds[threshold_idx]*100:.3f}% reached at partition #{i+1}/{len(parts)}")
            threshold_idx += 1

    print("\n\nDONE.")
