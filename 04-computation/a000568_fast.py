#!/usr/bin/env python3
"""
a000568_fast.py — Optimized computation of A000568 (nonisomorphic tournaments).

Key insight: DP over odd part sizes, tracking only the divisor-sum profile
needed for future gcd cross-terms. As we process larger parts, fewer
divisors remain relevant, pruning the state space dramatically.

Additional optimizations:
1. Only track divisors that will appear in future parts
2. Use integer arithmetic (multiply through by denominators)
3. Prune states that can't reach total n

Author: opus-2026-03-07-S46f
"""

from math import factorial, gcd
from fractions import Fraction
from time import perf_counter
from functools import lru_cache

def euler_totient(n):
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result

def a000568_fast(n):
    """Optimized DP computation of A000568(n).

    Process odd parts k = 1, 3, 5, ... in order.
    Track divisor sums S_d only for divisors d that divide some FUTURE part.
    """
    if n <= 1:
        return 1

    phi_cache = {i: euler_totient(i) for i in range(1, n + 1)}

    odd_parts = list(range(1, n + 1, 2))

    # For each part k, find which odd divisors d | k
    divisors_of = {}
    for k in odd_parts:
        divisors_of[k] = [d for d in range(1, k + 1) if k % d == 0 and d % 2 == 1]

    # When processing part k, we need to track S_d for all d that divide
    # some part k' >= k. After processing k, we can drop d values that
    # don't divide any remaining part.

    # Compute: for each divisor d, which is the largest part it divides?
    max_part_for_div = {}
    for k in odd_parts:
        for d in divisors_of[k]:
            max_part_for_div[d] = k  # last update = largest

    # For efficiency, track the set of "active" divisors at each stage
    def active_divisors_after(k):
        """Divisors still relevant after processing part k."""
        return sorted(d for d, maxk in max_part_for_div.items() if maxk > k)

    def active_divisors_at(k):
        """Divisors relevant when processing part k (including k's own)."""
        return sorted(d for d, maxk in max_part_for_div.items() if maxk >= k)

    # DP state: (total_size, tuple of S_d values for active divisors)
    # Weight: Fraction accumulator

    # Start: processing k=1
    # Active divisors at k=1: all odd divisors ≤ n that divide some part ≥ 1
    # = {1, 3, 5, 7, ...} up to n (d=1 divides everything)

    active = active_divisors_at(odd_parts[0])
    div_index = {d: i for i, d in enumerate(active)}

    init_state = tuple(0 for _ in active)
    dp = {(0, init_state): Fraction(1)}

    for ki, k in enumerate(odd_parts):
        if k > n:
            break

        divs_k = divisors_of[k]
        div_idx_k = [div_index[d] for d in divs_k if d in div_index]
        phi_k = [phi_cache[d] for d in divs_k if d in div_index]

        new_dp = {}
        for (total, state), weight in dp.items():
            max_m = (n - total) // k
            for m in range(0, max_m + 1):
                new_total = total + m * k

                # Pruning: can we still reach n with remaining parts?
                remaining = n - new_total
                if ki + 1 < len(odd_parts):
                    max_remaining_part = odd_parts[-1]  # largest available
                    # loosest check: remaining must be achievable
                    if remaining < 0:
                        continue
                else:
                    if remaining != 0:
                        continue

                # Compute Δt
                cross = sum(phi_k[i] * state[div_idx_k[i]] for i in range(len(div_idx_k)))
                delta_t = m * cross + m * (m - 1) * k // 2 + m * (k - 1) // 2

                factor = Fraction(1 << delta_t, k**m * factorial(m))

                # Update state
                new_state = list(state)
                for idx in div_idx_k:
                    new_state[idx] += m
                new_state = tuple(new_state)

                key = (new_total, new_state)
                val = weight * factor
                if key in new_dp:
                    new_dp[key] += val
                else:
                    new_dp[key] = val

        dp = new_dp

        # After processing k, compress state by dropping divisors
        # that won't appear in any future part
        if ki + 1 < len(odd_parts):
            next_k = odd_parts[ki + 1]
            new_active = active_divisors_at(next_k)

            if len(new_active) < len(active):
                # Need to project state: drop indices not in new_active
                new_div_index = {d: i for i, d in enumerate(new_active)}
                keep_indices = [div_index[d] for d in new_active if d in div_index]

                projected_dp = {}
                for (total, state), weight in dp.items():
                    new_state = tuple(state[i] for i in keep_indices)
                    key = (total, new_state)
                    if key in projected_dp:
                        projected_dp[key] += weight
                    else:
                        projected_dp[key] = weight

                dp = projected_dp
                active = new_active
                div_index = new_div_index

    # Sum over all states with total = n
    result = sum(w for (total, state), w in dp.items() if total == n)
    return int(result)


def a000568_fast_int(n):
    """Integer-arithmetic version (no Fraction overhead).

    Multiply through by the LCM of all possible z_λ denominators.
    Actually, track numerator and denominator separately, reducing periodically.
    """
    if n <= 1:
        return 1

    phi_cache = {i: euler_totient(i) for i in range(1, n + 1)}
    odd_parts = list(range(1, n + 1, 2))

    divisors_of = {}
    for k in odd_parts:
        divisors_of[k] = [d for d in range(1, k + 1) if k % d == 0 and d % 2 == 1]

    max_part_for_div = {}
    for k in odd_parts:
        for d in divisors_of[k]:
            max_part_for_div[d] = k

    def active_divisors_at(k):
        return sorted(d for d, maxk in max_part_for_div.items() if maxk >= k)

    active = active_divisors_at(odd_parts[0])
    div_index = {d: i for i, d in enumerate(active)}
    init_state = tuple(0 for _ in active)

    # dp[(total, state)] = (numerator, denominator) as integers
    dp = {(0, init_state): (1, 1)}

    for ki, k in enumerate(odd_parts):
        if k > n:
            break

        divs_k = divisors_of[k]
        div_idx_k = [div_index[d] for d in divs_k if d in div_index]
        phi_k = [phi_cache[d] for d in divs_k if d in div_index]

        new_dp = {}
        for (total, state), (num, den) in dp.items():
            max_m = (n - total) // k
            for m in range(0, max_m + 1):
                new_total = total + m * k
                if new_total > n:
                    break

                cross = sum(phi_k[i] * state[div_idx_k[i]] for i in range(len(div_idx_k)))
                delta_t = m * cross + m * (m - 1) * k // 2 + m * (k - 1) // 2

                new_num = num * (1 << delta_t)
                new_den = den * (k**m * factorial(m))

                # Reduce
                g = gcd(new_num, new_den)
                new_num //= g
                new_den //= g

                new_state = list(state)
                for idx in div_idx_k:
                    new_state[idx] += m
                new_state = tuple(new_state)

                key = (new_total, new_state)
                if key in new_dp:
                    on, od = new_dp[key]
                    # Add fractions: on/od + new_num/new_den
                    combined_num = on * new_den + new_num * od
                    combined_den = od * new_den
                    g = gcd(combined_num, combined_den)
                    new_dp[key] = (combined_num // g, combined_den // g)
                else:
                    new_dp[key] = (new_num, new_den)

        dp = new_dp

        # Project state
        if ki + 1 < len(odd_parts):
            next_k = odd_parts[ki + 1]
            new_active = active_divisors_at(next_k)
            if len(new_active) < len(active):
                new_div_index = {d: i for i, d in enumerate(new_active)}
                keep_indices = [div_index[d] for d in new_active if d in div_index]
                projected_dp = {}
                for (total, state), val in dp.items():
                    new_state = tuple(state[i] for i in keep_indices)
                    key = (total, new_state)
                    if key in projected_dp:
                        on, od = projected_dp[key]
                        nn, nd = val
                        combined_num = on * nd + nn * od
                        combined_den = od * nd
                        g = gcd(combined_num, combined_den)
                        projected_dp[key] = (combined_num // g, combined_den // g)
                    else:
                        projected_dp[key] = val
                dp = projected_dp
                active = new_active
                div_index = new_div_index

    result_num = 0
    result_den = 1
    for (total, state), (num, den) in dp.items():
        if total == n:
            result_num = result_num * den + num * result_den
            result_den = result_den * den
            g = gcd(result_num, result_den)
            result_num //= g
            result_den //= g

    assert result_den == 1, f"Result is not integer: {result_num}/{result_den}"
    return result_num


# Known values for validation
known = {0: 1, 1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
         9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168}

if __name__ == "__main__":
    print("=" * 60)
    print("A000568 FAST COMPUTATION")
    print("=" * 60)

    # Validate
    print("\nValidation (Fraction version):")
    for n in range(13):
        t0 = perf_counter()
        result = a000568_fast(n)
        dt = perf_counter() - t0
        expected = known.get(n, "?")
        match = "✓" if result == expected else f"✗ (expected {expected})"
        print(f"  a({n:2d}) = {result:>15}  [{dt:.4f}s]  {match}")

    print("\nValidation (integer version):")
    for n in range(13):
        t0 = perf_counter()
        result = a000568_fast_int(n)
        dt = perf_counter() - t0
        expected = known.get(n, "?")
        match = "✓" if result == expected else f"✗ (expected {expected})"
        print(f"  a({n:2d}) = {result:>15}  [{dt:.4f}s]  {match}")

    # Benchmark larger n
    print("\n" + "=" * 60)
    print("SCALING TEST")
    print("=" * 60)

    for n in [20, 30, 40, 50, 60, 70, 80, 90, 100]:
        t0 = perf_counter()
        result = a000568_fast(n)
        dt = perf_counter() - t0
        digits = len(str(result))
        print(f"  a({n:3d}) = {str(result)[:30]}{'...' if digits > 30 else '':3s} "
              f"({digits:4d} digits)  [{dt:8.3f}s]")
