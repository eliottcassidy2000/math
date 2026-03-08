#!/usr/bin/env python3
"""
a000568_int_v2.py — A000568 via integer DP with periodic GCD reduction.

Same as a000568_int but after each round, divides all weights by their
common GCD to keep integers small. Also tracks accumulated divisor to
recover the rational result.

Author: opus-2026-03-07-S46f
"""

import gmpy2
from gmpy2 import mpz
from time import perf_counter
from functools import reduce


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


def a000568_int_v2(n):
    """Compute A000568(n) using integer DP with periodic GCD reduction."""
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

    fact_cache = [mpz(1)]
    for i in range(1, n + 1):
        fact_cache.append(fact_cache[-1] * i)

    active = active_divisors_at(odd_parts[0])
    div_index = {d: i for i, d in enumerate(active)}
    init_state = tuple(0 for _ in active)

    dp = {(0, init_state): mpz(1)}
    D_so_far = mpz(1)

    for ki, k in enumerate(odd_parts):
        if k > n:
            break

        divs_k = divisors_of[k]
        div_idx_k = [div_index[d] for d in divs_k if d in div_index]
        phi_k = [phi_cache[d] for d in divs_k if d in div_index]

        k_mpz = mpz(k)
        max_m_k = n // k
        denom_budget = k_mpz ** max_m_k * fact_cache[max_m_k]

        new_dp = {}
        for (total, state), weight in dp.items():
            max_m_here = (n - total) // k
            for m in range(0, max_m_here + 1):
                new_total = total + m * k

                cross = sum(phi_k[i] * state[div_idx_k[i]] for i in range(len(div_idx_k)))
                delta_t = m * cross + m * (m - 1) * k // 2 + m * (k - 1) // 2

                # Integer factor: 2^Δt · k^{max_m_k - m} · max_m_k! / m!
                int_factor = (mpz(1) << delta_t) * k_mpz ** (max_m_k - m) * (fact_cache[max_m_k] // fact_cache[m])

                new_state = list(state)
                for idx in div_idx_k:
                    new_state[idx] += m
                new_state = tuple(new_state)

                key = (new_total, new_state)
                val = weight * int_factor
                if key in new_dp:
                    new_dp[key] += val
                else:
                    new_dp[key] = val

        D_so_far *= denom_budget
        dp = new_dp

        # Periodic GCD reduction: divide all weights by their common GCD
        if dp:
            vals = list(dp.values())
            g = vals[0]
            for v in vals[1:]:
                g = gmpy2.gcd(g, v)
                if g == 1:
                    break
            if g > 1:
                dp = {k: v // g for k, v in dp.items()}
                D_so_far //= g

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
                        projected_dp[key] += val
                    else:
                        projected_dp[key] = val
                dp = projected_dp
                active = new_active
                div_index = new_div_index

    result = mpz(0)
    for (total, state), val in dp.items():
        if total == n:
            result += val

    assert result % D_so_far == 0, f"Result not divisible by D"
    return int(result // D_so_far)


known = {0: 1, 1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
         9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168}

if __name__ == "__main__":
    print("=" * 70)
    print("A000568 via integer DP with GCD reduction")
    print("=" * 70, flush=True)

    print("\nValidation:", flush=True)
    for n_val in range(13):
        t0 = perf_counter()
        result = a000568_int_v2(n_val)
        dt = perf_counter() - t0
        expected = known.get(n_val, "?")
        match = "OK" if result == expected else f"FAIL (got {result}, expected {expected})"
        print(f"  a({n_val:2d}) = {result:>15}  [{dt:.4f}s]  {match}", flush=True)

    print("\nScaling:", flush=True)
    for n_test in [20, 30, 40, 50, 60, 70, 76]:
        t0 = perf_counter()
        result = a000568_int_v2(n_test)
        dt = perf_counter() - t0
        s = str(result)
        print(f"  a({n_test:3d}) ({len(s):4d} digits)  [{dt:8.3f}s]", flush=True)
