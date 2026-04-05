#!/usr/bin/env python3
"""
a000568_truncated_s24.py — opus-2026-04-05-S24

WILD NEW APPROACH: Truncated partition enumeration for A000568.

KEY INSIGHT: In the Burnside formula
  A000568(n) = Σ_{λ odd-part partition of n} (1/z_λ) × 2^{e(λ)}

the identity partition (1^n) dominates with contribution 2^{C(n,2)}/n!.
Every other partition is EXPONENTIALLY suppressed: a single 3-cycle costs
a factor of 2^{-(2n-4)}, a 5-cycle costs 2^{-(4n-12)}, etc.

The total is an exact sum of rational numbers. To get the exact integer,
we only need partitions whose cumulative contribution exceeds 0.5 in the
last decimal place. This requires a "departure budget" of only
B_max ≈ 0.35n, giving p(B_max/2) ≈ p(n/6) partitions instead of p(n).

SPEEDUP at n=200: p(200) ≈ 4×10¹² → p(33) ≈ 10,000. Factor ~10⁸.
SPEEDUP at n=500: p(500) ≈ 10²² → p(88) ≈ 3×10⁶. Factor ~10¹⁶.

Method:
1. Enumerate "departure configurations" — multisets of odd parts ≥ 3
   with total departure Σ(l_i - 1) ≤ B_max
2. For each configuration, compute the Burnside contribution using mpmath
3. Sum all contributions with sufficient precision
4. Round to nearest integer

This is a POLYNOMIAL-TIME algorithm (per digit of output) when B_max grows
linearly with n. The current GMP approach is exponential-time (partition count).
"""

from mpmath import mp, mpf, fac, log, ceil, power, floor as mpfloor
from math import gcd, comb, factorial
from collections import Counter
from itertools import combinations_with_replacement
import time
import sys

sys.stdout.reconfigure(encoding='utf-8')

def compute_e(parts_counter, n):
    """Compute edge exponent e(λ) for a partition given as {size: multiplicity}.

    e(λ) = Σ_k m_k * floor(k/2)
          + Σ_k C(m_k, 2) * k
          + Σ_{k<l} m_k * m_l * gcd(k, l)

    where parts_counter = {k: m_k} including 1-parts.
    """
    e = 0
    sizes = sorted(parts_counter.keys())
    for k in sizes:
        mk = parts_counter[k]
        e += mk * (k // 2)             # Σ m_k * floor(k/2)
        e += mk * (mk - 1) // 2 * k   # Σ C(m_k, 2) * k
    for i, k in enumerate(sizes):
        for l in sizes[i+1:]:
            e += parts_counter[k] * parts_counter[l] * gcd(k, l)
    return e

def compute_z(parts_counter):
    """Compute z_λ = Π_k k^m_k * m_k!"""
    z = 1
    for k, mk in parts_counter.items():
        z *= k**mk * factorial(mk)
    return z

def deficit(parts_counter, n):
    """Compute deficit = C(n,2) - e(λ)."""
    return comb(n, 2) - compute_e(parts_counter, n)

def enumerate_departure_configs(budget_max, n):
    """Enumerate all multisets of odd parts ≥ 3 with total departure ≤ budget_max
    and total sum ≤ n.

    Each odd part k ≥ 3 has departure k-1 (which is even, ≥ 2).
    A configuration is a Counter {k: multiplicity}.

    Yields (config, total_used) where total_used = Σ k*m_k ≤ n.
    """
    # Recursive enumeration
    odd_sizes = list(range(3, n+1, 2))  # 3, 5, 7, ..., n (or n-1)

    def recurse(idx, remaining_budget, remaining_n):
        """Yield configurations using sizes[idx:]."""
        yield Counter()  # empty configuration (departure 0)

        for i in range(idx, len(odd_sizes)):
            k = odd_sizes[i]
            dep = k - 1  # departure per copy

            if dep > remaining_budget:
                break  # larger sizes have even more departure

            # Add 1, 2, 3, ... copies of size k
            for m in range(1, remaining_budget // dep + 1):
                total_dep = m * dep
                total_used = m * k

                if total_dep > remaining_budget or total_used > remaining_n:
                    break

                # This configuration: m copies of k, plus whatever comes after
                base = Counter({k: m})
                for sub_config in recurse(i + 1, remaining_budget - total_dep,
                                          remaining_n - total_used):
                    config = base + sub_config
                    yield config

    # Remove the empty config (it's the identity, handled separately)
    first = True
    for config in recurse(0, budget_max, n):
        if first:
            first = False
            yield config  # identity (empty)
        else:
            yield config

def a000568_truncated(n, safety_digits=50):
    """Compute A000568(n) using truncated partition enumeration."""

    cn2 = n * (n - 1) // 2  # C(n, 2)

    # Estimate number of digits in result
    # A000568(n) ≈ 2^{C(n,2)} / n!
    # log10 ≈ C(n,2) * log10(2) - log10(n!)
    import math
    log10_result = cn2 * math.log10(2) - sum(math.log10(k) for k in range(1, n+1))
    num_digits = int(log10_result) + 1

    # Set mpmath precision
    mp.dps = num_digits + safety_digits + 50

    print(f"  n={n}: est {num_digits} digits, precision={mp.dps} dps")

    # Estimate budget needed
    # Each 3-cycle reduces contribution by factor ~2^{-(2n-4)} ≈ 10^{-(0.6n)}
    # Need levels until error < 0.5
    # A very conservative estimate:
    deficit_per_3cycle = 2 * (n - 2)  # deficit of one 3-cycle
    if deficit_per_3cycle > 0:
        levels_needed = int(math.ceil(num_digits * math.log(10) / math.log(2) / deficit_per_3cycle)) + 5
    else:
        levels_needed = n

    budget_max = levels_needed * 2  # departure budget (each 3-cycle has departure 2)
    budget_max = min(budget_max, n)  # can't exceed n

    print(f"  Deficit per 3-cycle: {deficit_per_3cycle}")
    print(f"  Levels needed: {levels_needed}")
    print(f"  Departure budget: {budget_max}")

    # Count and enumerate configurations
    t0 = time.time()
    total = mpf(0)
    config_count = 0
    max_deficit_seen = 0

    # Precompute n! as mpf
    nfact = fac(n)

    for config in enumerate_departure_configs(budget_max, n):
        # Build full partition: config + enough 1-parts
        used = sum(k * m for k, m in config.items())
        ones = n - used
        if ones < 0:
            continue

        # Full partition counter
        parts = Counter(config)
        if ones > 0:
            parts[1] = ones

        # Compute e(λ)
        e = compute_e(parts, n)

        # Compute z_λ
        z = compute_z(parts)

        # Contribution = (1/z_λ) × 2^{e(λ)}
        # Use mpmath: contribution = power(2, e) / z
        contrib = power(2, e) / mpf(z)
        total += contrib

        config_count += 1
        d = cn2 - e
        max_deficit_seen = max(max_deficit_seen, d)

        if config_count <= 10 or config_count % 10000 == 0:
            ratio_to_identity = float(contrib / (power(2, cn2) / nfact)) if cn2 > 0 else 1
            if config_count <= 10:
                parts_str = dict(sorted(((k, v) for k, v in config.items()), key=lambda x: x[0]))
                print(f"    #{config_count}: parts={parts_str if parts_str else '{}'}, "
                      f"e={e}, deficit={d}, ratio_to_id={ratio_to_identity:.3e}")
            elif config_count % 10000 == 0:
                print(f"    #{config_count} configs processed, max_deficit={max_deficit_seen}")

    elapsed = time.time() - t0

    # A000568 = Σ (1/z_λ) × 2^{e(λ)} — no extra n! division needed
    # (The Burnside formula already includes the 1/n! via z_λ)
    result = total

    # Round to nearest integer
    result_int = int(mpfloor(result + mpf('0.5')))

    # Verify: check the fractional part is close to 0 or 1
    frac_part = float(result - mpfloor(result))
    if frac_part > 0.5:
        frac_part = 1 - frac_part

    print(f"\n  Configs enumerated: {config_count}")
    print(f"  Max deficit seen: {max_deficit_seen}")
    print(f"  Time: {elapsed:.3f}s")
    print(f"  Fractional part distance from integer: {frac_part:.2e}")
    print(f"  Result digits: {len(str(result_int))}")

    return result_int, elapsed, config_count


# ================================================================
# Verification and benchmarking
# ================================================================

# Known values for verification
KNOWN = {
    3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880, 9: 191536,
    10: 9733056, 11: 903753248, 12: 154108311168,
    13: 48542114686912, 14: 28401423719122304,
    15: 30995890806033613824,
    20: 4951364705546250940976699872,
}

print("=" * 70)
print("A000568 via TRUNCATED PARTITION ENUMERATION")
print("=" * 70)

results = {}

for n in [5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 150, 200]:
    print(f"\n{'='*50}")
    print(f"n = {n}")
    print(f"{'='*50}")

    result, elapsed, num_configs = a000568_truncated(n)

    # Verify against known values
    if n in KNOWN:
        correct = result == KNOWN[n]
        print(f"  Known value: {KNOWN[n]}")
        print(f"  Our result:  {result}")
        print(f"  CORRECT: {correct}")
        if not correct:
            print(f"  *** MISMATCH! Difference: {result - KNOWN[n]} ***")
    else:
        print(f"  Result: {str(result)[:80]}...{str(result)[-20:]}" if len(str(result)) > 100
              else f"  Result: {result}")

    results[n] = {
        'result': result,
        'time': elapsed,
        'configs': num_configs,
        'digits': len(str(result)),
    }

# Summary
print(f"\n\n{'='*70}")
print("BENCHMARK SUMMARY")
print(f"{'='*70}")
print(f"  {'n':>5} | {'configs':>10} | {'time (s)':>10} | {'digits':>7} | {'configs/s':>10} | vs p(n)")
print(f"  {'-'*5}-|-{'-'*10}-|-{'-'*10}-|-{'-'*7}-|-{'-'*10}-|-------")
for n in sorted(results.keys()):
    r = results[n]
    rate = r['configs'] / r['time'] if r['time'] > 0 else float('inf')
    # Rough estimate of p(n)
    import math
    pn_est = math.exp(math.pi * math.sqrt(2*n/3)) / (4*n*math.sqrt(3)) if n > 1 else 1
    speedup = pn_est / r['configs'] if r['configs'] > 0 else 0
    print(f"  {n:>5} | {r['configs']:>10} | {r['time']:>10.3f} | {r['digits']:>7} | {rate:>10.0f} | {speedup:>.1e}")
