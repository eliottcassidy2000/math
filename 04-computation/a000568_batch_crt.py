#!/usr/bin/env python3
"""
a000568_batch_crt.py — Vectorized CRT for A000568.

Key idea: Run the DP ONCE, but at each state compute residues mod ALL primes
simultaneously using numpy arrays. This avoids running the DP num_primes times.

Each dp value is a numpy array of shape (num_primes,) holding residues mod p_i.
All modular arithmetic is vectorized.

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

def a000568_batch_crt(n, verbose=True):
    """Compute A000568(n) using batched CRT: one DP pass, vectorized over primes."""
    if n <= 1:
        return 1

    import math
    est_bits = int(n * (n - 1) / 2 - n * math.log2(max(n, 2)) + n) + 100
    if est_bits < 100:
        est_bits = 100

    # Use large primes (~30 bits) for fewer CRT steps
    prime_start = max(n, 2**30)  # primes near 2^30 ≈ 1B
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

    # Precompute vectorized modular inverses
    # pow(base, exp, mod) doesn't vectorize in numpy, but we can batch it
    def vec_mod_inv(x, P):
        """Compute x^{-1} mod each p in P. x is a Python int."""
        result = np.empty(len(P), dtype=np.int64)
        for j, p in enumerate(primes):
            result[j] = pow(x % p, p - 2, p)
        return result

    def vec_pow2(exp, P):
        """Compute 2^exp mod each p in P. exp is a Python int."""
        result = np.empty(len(P), dtype=np.int64)
        for j, p in enumerate(primes):
            result[j] = pow(2, exp, p)
        return result

    def vec_mul(a, b, P):
        """Element-wise (a * b) mod P using 128-bit intermediates."""
        # Use Python int to avoid overflow, then back to array
        # Actually, for 30-bit primes, 64-bit product fits: 30+30=60 < 63
        return (a * b) % P

    def vec_add(a, b, P):
        return (a + b) % P

    # Precompute inverses for k^m and m! that we'll need
    inv_k_cache = {}  # k -> vec of k^{-1} mod P
    inv_fact_cache = {}  # m -> vec of (m!)^{-1} mod P

    # DP: states -> numpy arrays of residues
    active = active_divisors_at(odd_parts[0])
    div_index = {d: i for i, d in enumerate(active)}
    init_state = tuple(0 for _ in active)
    ones = np.ones(num_primes, dtype=np.int64)

    dp = {(0, init_state): ones.copy()}

    for ki, k in enumerate(odd_parts):
        if k > n:
            break

        divs_k = divisors_of[k]
        div_idx_k = [div_index[d] for d in divs_k if d in div_index]
        phi_k = [phi_cache[d] for d in divs_k if d in div_index]

        max_m = n // k

        # Precompute k^m mod P and (m!)^{-1} mod P
        k_pow_vecs = [ones.copy()]  # k^0 = 1
        k_mod_P = np.array([k % p for p in primes], dtype=np.int64)
        for m in range(1, max_m + 1):
            k_pow_vecs.append(vec_mul(k_pow_vecs[-1], k_mod_P, P))

        fact_inv_vecs = [ones.copy()]  # (0!)^{-1} = 1
        fact = 1
        for m in range(1, max_m + 1):
            fact = fact * m
            fact_inv_vecs.append(vec_mod_inv(fact, P))

        # Precompute k_pow_inv: (k^m)^{-1} mod P
        k_pow_inv_vecs = [ones.copy()]
        k_inv = vec_mod_inv(k, P)
        for m in range(1, max_m + 1):
            k_pow_inv_vecs.append(vec_mul(k_pow_inv_vecs[-1], k_inv, P))

        # Precompute 2^{self_t(m)} for each m, where self_t = m(m-1)k/2 + m(k-1)/2
        self_t_vecs = []
        for m in range(0, max_m + 1):
            st = m * (m - 1) * k // 2 + m * (k - 1) // 2
            self_t_vecs.append(vec_pow2(st, P))

        # Combined coefficient: 2^{self_t} * (k^m)^{-1} * (m!)^{-1} mod P
        coeff_vecs = []
        for m in range(0, max_m + 1):
            coeff_vecs.append(vec_mul(vec_mul(self_t_vecs[m], k_pow_inv_vecs[m], P), fact_inv_vecs[m], P))

        # Cache 2^cross vectors: cross -> vec
        pow2_cross_cache = {}

        new_dp = {}
        for (total, state), weight in dp.items():
            max_m_here = (n - total) // k

            cross = sum(phi_k[i] * state[div_idx_k[i]] for i in range(len(div_idx_k)))

            # Precompute (2^cross)^m for m=0..max_m_here
            if cross not in pow2_cross_cache:
                pow2_cross_cache[cross] = vec_pow2(cross, P)
            two_c = pow2_cross_cache[cross]
            # (2^cross)^m via repeated multiply
            two_c_pow = ones.copy()  # m=0

            for m in range(0, max_m_here + 1):
                new_total = total + m * k

                # factor = (2^cross)^m * coeff[m]
                factor = vec_mul(two_c_pow, coeff_vecs[m], P)

                new_state = list(state)
                for idx in div_idx_k:
                    new_state[idx] += m
                new_state = tuple(new_state)

                key = (new_total, new_state)
                val = vec_mul(weight, factor, P)
                if key in new_dp:
                    new_dp[key] = vec_add(new_dp[key], val, P)
                else:
                    new_dp[key] = val

                # Update for next m: (2^cross)^{m+1} = (2^cross)^m * 2^cross
                two_c_pow = vec_mul(two_c_pow, two_c, P)

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
                        projected_dp[key] = vec_add(projected_dp[key], val, P)
                    else:
                        projected_dp[key] = val
                dp = projected_dp
                active = new_active
                div_index = new_div_index

        if verbose and k <= 21:
            print(f"    k={k:3d}: {len(dp)} states", flush=True)

    # Sum over states with total = n
    result_vec = np.zeros(num_primes, dtype=np.int64)
    for (total, state), val in dp.items():
        if total == n:
            result_vec = vec_add(result_vec, val, P)

    # CRT reconstruction
    residues = [int(result_vec[i]) for i in range(num_primes)]

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


# Known values
known = {0: 1, 1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
         9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168}

if __name__ == "__main__":
    print("=" * 70)
    print("A000568 via BATCHED CRT (vectorized over primes)")
    print("=" * 70)

    print("\nValidation:", flush=True)
    for n in range(13):
        t0 = perf_counter()
        result = a000568_batch_crt(n, verbose=False)
        dt = perf_counter() - t0
        expected = known.get(n, "?")
        match = "OK" if result == expected else f"FAIL (got {result}, expected {expected})"
        print(f"  a({n:2d}) = {result:>15}  [{dt:.4f}s]  {match}", flush=True)

    print("\n" + "=" * 70)
    print("SCALING TEST")
    print("=" * 70, flush=True)

    for n in [20, 30, 40, 50, 60, 70, 76, 80, 90, 100]:
        t0 = perf_counter()
        result = a000568_batch_crt(n, verbose=True)
        dt = perf_counter() - t0
        s = str(result)
        print(f"  a({n:3d}) = {s[:40]}{'...' if len(s) > 40 else '':3s} "
              f"({len(s):4d} digits)  [{dt:8.3f}s]", flush=True)

    # Verify n=76 against OEIS
    print("\n" + "=" * 70, flush=True)
    print("VERIFICATION vs OEIS", flush=True)
    print("=" * 70, flush=True)
    oeis_76 = 457153791459731763873714717642097124535855865100705012123154296782586066392599174177717111708812242945048119323356894628050571858925201215481427701704774498054845450499257140261101951670038378795297655176726665694457285914072690084154974004045609532979011190165711892081116449736469812420129338461848172656662018517999846546136461519522876668617723623760713964531731492213163219661503132797160477050776498667637823035846487741247681442428085385480808545797005707547708306606906039466675985395538562415540138435374979216235780263415436443048770796035424757246769390586508396300108786258541736828941774383334172760307820994535160154251942298642998619446487984483719749964075470321740294755571922446326714092505683852691999256324957470064758984015872
    result_76 = a000568_batch_crt(76, verbose=False)
    print(f"  a(76) matches OEIS: {result_76 == oeis_76}", flush=True)
