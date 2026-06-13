#!/usr/bin/env python3
"""
a000568_int_v3.py — A000568 via integer DP, optimized Python.

Optimizations over v2:
1. Precompute all int_factors for each (k, m) pair
2. Use list-based state representation
3. Batch GCD using a sampled subset for speed
4. Inline critical operations

Author: opus-2026-03-07-S46f
"""

import gmpy2
from gmpy2 import mpz
from time import perf_counter


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


def a000568_int_v3(n):
    """Compute A000568(n) using optimized integer DP."""
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

        # Precompute k^{max_m_k - m} * (max_m_k! / m!) for each m
        # These don't depend on the state, only on m
        k_pow_complement = [mpz(0)] * (max_m_k + 1)
        fact_ratio = [mpz(0)] * (max_m_k + 1)
        for m in range(max_m_k + 1):
            k_pow_complement[m] = k_mpz ** (max_m_k - m)
            fact_ratio[m] = fact_cache[max_m_k] // fact_cache[m]

        # Precompute base_factor[m] = k^{max_m_k - m} * max_m_k! / m!
        base_factor = [k_pow_complement[m] * fact_ratio[m] for m in range(max_m_k + 1)]

        # Self-contribution: m(m-1)k/2 + m(k-1)/2
        self_delta = [m * (m - 1) * k // 2 + m * (k - 1) // 2 for m in range(max_m_k + 1)]

        num_dk = len(div_idx_k)
        new_dp = {}

        for (total, state), weight in dp.items():
            max_m_here = (n - total) // k

            for m in range(max_m_here + 1):
                # Compute cross contribution
                cross = 0
                for i in range(num_dk):
                    cross += phi_k[i] * state[div_idx_k[i]]

                delta_t = m * cross + self_delta[m]

                # int_factor = 2^Δt * base_factor[m]
                int_factor = (mpz(1) << delta_t) * base_factor[m]

                # Build new state
                new_state = list(state)
                for i in range(num_dk):
                    new_state[div_idx_k[i]] += m

                key = (total + m * k, tuple(new_state))
                val = weight * int_factor
                if key in new_dp:
                    new_dp[key] += val
                else:
                    new_dp[key] = val

        D_so_far *= denom_budget
        dp = new_dp

        # GCD reduction using sampling for speed
        if dp:
            vals = list(dp.values())
            nv = len(vals)
            if nv <= 100:
                g = vals[0]
                for v in vals[1:]:
                    g = gmpy2.gcd(g, v)
                    if g == 1:
                        break
            else:
                # Sample ~50 values for GCD estimate, then verify
                import random
                sample_idx = random.sample(range(nv), min(50, nv))
                g = vals[sample_idx[0]]
                for idx in sample_idx[1:]:
                    g = gmpy2.gcd(g, vals[idx])
                    if g == 1:
                        break
                if g > 1:
                    # Verify against all
                    for v in vals:
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

    assert result % D_so_far == 0
    return int(result // D_so_far)


known = {0: 1, 1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
         9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168}

if __name__ == "__main__":
    print("=" * 70)
    print("A000568 via optimized integer DP (v3)")
    print("=" * 70, flush=True)

    print("\nValidation:", flush=True)
    for n_val in range(13):
        t0 = perf_counter()
        result = a000568_int_v3(n_val)
        dt = perf_counter() - t0
        expected = known.get(n_val, "?")
        match = "OK" if result == expected else f"FAIL (got {result}, expected {expected})"
        print(f"  a({n_val:2d}) = {result:>15}  [{dt:.4f}s]  {match}", flush=True)

    print("\nBenchmark:", flush=True)
    for n_test in [20, 30, 40, 50, 60, 70, 76, 80]:
        t0 = perf_counter()
        result = a000568_int_v3(n_test)
        dt = perf_counter() - t0
        s = str(result)
        print(f"  a({n_test:3d}) ({len(s):4d} digits)  [{dt:8.3f}s]", flush=True)
