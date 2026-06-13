#!/usr/bin/env python3
"""
a000568_gmpy2.py — A000568 via DP with gmpy2 for fast big-integer arithmetic.

Uses gmpy2.mpq (rational numbers backed by GMP) instead of Python fractions.Fraction.
Expected 3-10x speedup from GMP's optimized arbitrary-precision arithmetic.

Author: opus-2026-03-07-S46f
"""

import gmpy2
from gmpy2 import mpz, mpq
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

def a000568_gmpy2(n):
    """Compute A000568(n) using DP with gmpy2 rational arithmetic."""
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

    dp = {(0, init_state): mpq(1)}

    # Precompute factorials as mpz
    fact_cache = [mpz(1)]
    for i in range(1, n + 1):
        fact_cache.append(fact_cache[-1] * i)

    for ki, k in enumerate(odd_parts):
        if k > n:
            break

        divs_k = divisors_of[k]
        div_idx_k = [div_index[d] for d in divs_k if d in div_index]
        phi_k = [phi_cache[d] for d in divs_k if d in div_index]

        k_mpz = mpz(k)

        new_dp = {}
        for (total, state), weight in dp.items():
            max_m = (n - total) // k
            for m in range(0, max_m + 1):
                new_total = total + m * k

                cross = sum(phi_k[i] * state[div_idx_k[i]] for i in range(len(div_idx_k)))
                delta_t = m * cross + m * (m - 1) * k // 2 + m * (k - 1) // 2

                # factor = 2^{delta_t} / (k^m * m!)
                num = mpz(1) << delta_t
                den = k_mpz ** m * fact_cache[m]
                factor = mpq(num, den)

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

    result = mpq(0)
    for (total, state), val in dp.items():
        if total == n:
            result += val

    return int(result)


known = {0: 1, 1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
         9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168}

if __name__ == "__main__":
    print("=" * 70)
    print("A000568 via gmpy2 DP")
    print("=" * 70, flush=True)

    print("\nValidation:", flush=True)
    for n in range(13):
        t0 = perf_counter()
        result = a000568_gmpy2(n)
        dt = perf_counter() - t0
        expected = known.get(n, "?")
        match = "OK" if result == expected else f"FAIL (got {result}, expected {expected})"
        print(f"  a({n:2d}) = {result:>15}  [{dt:.4f}s]  {match}", flush=True)

    oeis_76 = 457153791459731763873714717642097124535855865100705012123154296782586066392599174177717111708812242945048119323356894628050571858925201215481427701704774498054845450499257140261101951670038378795297655176726665694457285914072690084154974004045609532979011190165711892081116449736469812420129338461848172656662018517999846546136461519522876668617723623760713964531731492213163219661503132797160477050776498667637823035846487741247681442428085385480808545797005707547708306606906039466675985395538562415540138435374979216235780263415436443048770796035424757246769390586508396300108786258541736828941774383334172760307820994535160154251942298642998619446487984483719749964075470321740294755571922446326714092505683852691999256324957470064758984015872

    print("\nScaling:", flush=True)
    for n_test in [20, 30, 40, 50, 60, 70, 76, 80, 90, 100]:
        t0 = perf_counter()
        result = a000568_gmpy2(n_test)
        dt = perf_counter() - t0
        s = str(result)
        extra = ""
        if n_test == 76:
            extra = f"  OEIS={'OK' if result == oeis_76 else 'FAIL'}"
        print(f"  a({n_test:3d}) ({len(s):4d} digits)  [{dt:8.3f}s]{extra}", flush=True)
