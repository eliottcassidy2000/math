#!/usr/bin/env python3
"""
a000568_batch_crt_v2.py — Optimized batched CRT for A000568.

Improvements over v1:
- Precompute pow(2, e, p) using discrete log tables for small-order groups
- Cache pow2 results more aggressively
- Use ctypes or manual batching for pow operations

Author: opus-2026-03-07-S46f
"""

import numpy as np
from math import factorial, gcd
from time import perf_counter
import sys

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

def primes_above(start, count):
    """Generate `count` primes starting above `start`."""
    primes = []
    n = start + 1
    if n % 2 == 0:
        n += 1
    while len(primes) < count:
        is_p = True
        i = 3
        while i * i <= n:
            if n % i == 0:
                is_p = False
                break
            i += 2
        if n > 1 and is_p:
            primes.append(n)
        n += 2
    return primes


def a000568_batch_crt_v2(n, verbose=True):
    """Compute A000568(n) using batched CRT with aggressive caching."""
    if n <= 1:
        return 1

    import math
    est_bits = int(n * (n - 1) / 2 - n * math.log2(max(n, 2)) + n) + 100
    if est_bits < 100:
        est_bits = 100

    prime_start = max(n, 2**30)
    prime_bits = 30
    num_primes = est_bits // prime_bits + 5

    primes = primes_above(prime_start, num_primes)
    P = np.array(primes, dtype=np.int64)

    if verbose:
        print(f"  n={n}: {num_primes} primes, est {est_bits} bits", flush=True)

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

    NP = num_primes
    ones = np.ones(NP, dtype=np.int64)

    # Precompute 2 mod P
    two_mod_P = np.array([2 % p for p in primes], dtype=np.int64)

    # Build pow(2, e, P) cache using repeated squaring in numpy
    # For each prime p, ord_p(2) | (p-1). Max exponent we'll see is ~n^2.
    # Strategy: precompute 2^{2^i} mod P for i=0..60, then combine bits.
    max_exp = n * n  # rough upper bound on delta_t

    # Precompute 2^{2^i} mod P for i = 0, 1, ..., 60
    two_pow2i = [two_mod_P.copy()]  # 2^{2^0} = 2
    for i in range(60):
        two_pow2i.append((two_pow2i[-1] * two_pow2i[-1]) % P)

    def vec_pow2_fast(exp):
        """Compute 2^exp mod P using precomputed squares."""
        result = ones.copy()
        e = exp
        i = 0
        while e > 0:
            if e & 1:
                result = (result * two_pow2i[i]) % P
            e >>= 1
            i += 1
        return result

    # Cache for vec_pow2_fast results
    pow2_cache = {}
    def vec_pow2_cached(exp):
        if exp in pow2_cache:
            return pow2_cache[exp]
        r = vec_pow2_fast(exp)
        if len(pow2_cache) < 100000:  # limit cache size
            pow2_cache[exp] = r
        return r

    def vec_mod_inv(x):
        """Compute x^{-1} mod each p in P."""
        result = np.empty(NP, dtype=np.int64)
        for j, p in enumerate(primes):
            result[j] = pow(x % p, p - 2, p)
        return result

    def vec_mul(a, b):
        return (a * b) % P

    def vec_add(a, b):
        return (a + b) % P

    # DP
    active = active_divisors_at(odd_parts[0])
    div_index = {d: i for i, d in enumerate(active)}
    init_state = tuple(0 for _ in active)

    dp = {(0, init_state): ones.copy()}

    for ki, k in enumerate(odd_parts):
        if k > n:
            break

        divs_k = divisors_of[k]
        div_idx_k = [div_index[d] for d in divs_k if d in div_index]
        phi_k = [phi_cache[d] for d in divs_k if d in div_index]

        max_m = n // k

        # Precompute self_t coefficients for each m
        # self_t(m) = m*(m-1)*k//2 + m*(k-1)//2
        self_t_vecs = []
        for m in range(max_m + 1):
            st = m * (m - 1) * k // 2 + m * (k - 1) // 2
            self_t_vecs.append(vec_pow2_cached(st))

        # Precompute (k^m * m!)^{-1} mod P
        k_inv = vec_mod_inv(k)
        k_pow_inv = [ones.copy()]
        for m in range(1, max_m + 1):
            k_pow_inv.append(vec_mul(k_pow_inv[-1], k_inv))

        fact_inv = [ones.copy()]
        fact = 1
        for m in range(1, max_m + 1):
            fact *= m
            fact_inv.append(vec_mod_inv(fact))

        # Combined: 2^{self_t} * (k^m)^{-1} * (m!)^{-1}
        coeff = []
        for m in range(max_m + 1):
            coeff.append(vec_mul(vec_mul(self_t_vecs[m], k_pow_inv[m]), fact_inv[m]))

        new_dp = {}
        for (total, state), weight in dp.items():
            max_m_here = (n - total) // k

            cross = sum(phi_k[i] * state[div_idx_k[i]] for i in range(len(div_idx_k)))

            # 2^cross mod P (cached)
            two_c = vec_pow2_cached(cross)
            # (2^cross)^m via repeated multiply
            two_c_pow = ones.copy()

            for m in range(max_m_here + 1):
                new_total = total + m * k

                # factor = (2^cross)^m * coeff[m]
                factor = vec_mul(two_c_pow, coeff[m])

                new_state = list(state)
                for idx in div_idx_k:
                    new_state[idx] += m
                new_state = tuple(new_state)

                key = (new_total, new_state)
                val = vec_mul(weight, factor)
                if key in new_dp:
                    new_dp[key] = vec_add(new_dp[key], val)
                else:
                    new_dp[key] = val

                two_c_pow = vec_mul(two_c_pow, two_c)

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
                        projected_dp[key] = vec_add(projected_dp[key], val)
                    else:
                        projected_dp[key] = val
                dp = projected_dp
                active = new_active
                div_index = new_div_index

        if verbose and k <= 21:
            print(f"    k={k:3d}: {len(dp)} states", flush=True)

    # Sum
    result_vec = np.zeros(NP, dtype=np.int64)
    for (total, state), val in dp.items():
        if total == n:
            result_vec = vec_add(result_vec, val)

    # CRT
    residues = [int(result_vec[i]) for i in range(NP)]
    result = residues[0]
    modulus = primes[0]
    for i in range(1, len(residues)):
        r = residues[i]
        p = primes[i]
        diff = (r - result % p) % p
        t = diff * pow(modulus, p - 2, p) % p
        result = result + modulus * t
        modulus = modulus * p

    if result < 0:
        result += modulus

    return result


known = {0: 1, 1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
         9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168}

if __name__ == "__main__":
    print("=" * 70)
    print("A000568 via BATCHED CRT v2 (fast pow2)")
    print("=" * 70, flush=True)

    print("\nValidation:", flush=True)
    for n in range(13):
        t0 = perf_counter()
        result = a000568_batch_crt_v2(n, verbose=False)
        dt = perf_counter() - t0
        expected = known.get(n, "?")
        match = "OK" if result == expected else f"FAIL (got {result}, expected {expected})"
        print(f"  a({n:2d}) = {result:>15}  [{dt:.4f}s]  {match}", flush=True)

    print("\n" + "=" * 70)
    print("SCALING TEST")
    print("=" * 70, flush=True)

    oeis_76 = 457153791459731763873714717642097124535855865100705012123154296782586066392599174177717111708812242945048119323356894628050571858925201215481427701704774498054845450499257140261101951670038378795297655176726665694457285914072690084154974004045609532979011190165711892081116449736469812420129338461848172656662018517999846546136461519522876668617723623760713964531731492213163219661503132797160477050776498667637823035846487741247681442428085385480808545797005707547708306606906039466675985395538562415540138435374979216235780263415436443048770796035424757246769390586508396300108786258541736828941774383334172760307820994535160154251942298642998619446487984483719749964075470321740294755571922446326714092505683852691999256324957470064758984015872

    for n_test in [20, 30, 40, 50, 60, 70, 76, 80, 90, 100]:
        t0 = perf_counter()
        result = a000568_batch_crt_v2(n_test, verbose=True)
        dt = perf_counter() - t0
        s = str(result)
        extra = ""
        if n_test == 76:
            extra = f"  OEIS={'OK' if result == oeis_76 else 'FAIL'}"
        print(f"  a({n_test:3d}) = {s[:40]}{'...' if len(s) > 40 else '':3s} "
              f"({len(s):4d} digits)  [{dt:8.3f}s]{extra}", flush=True)
