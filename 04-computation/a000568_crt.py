#!/usr/bin/env python3
"""
a000568_crt.py — Ultra-fast A000568 via DP + Chinese Remainder Theorem.

The DP over odd part sizes is correct but slow for large n due to
Fraction arithmetic. Solution: compute a(n) mod p for many primes p
using machine-word arithmetic, then reconstruct via CRT.

Key insight: In the DP, the weight at each state is a rational number
whose denominator divides D = Π_{k odd} k^{floor(n/k)} · floor(n/k)!
Working mod p (prime > n), D is invertible, so we can use modular inverses.

Author: opus-2026-03-07-S46f
"""

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

def a000568_mod_p_dp(n, p):
    """Compute A000568(n) mod p using the DP method with modular arithmetic.

    Requires p > n (so factorials are invertible).
    """
    if n <= 1:
        return 1 % p

    phi_cache = {}
    for i in range(1, n + 1):
        phi_cache[i] = euler_totient(i)

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

    # Precompute modular inverses for k and m! that we'll need
    inv_cache = {}
    def mod_inv(x):
        x = x % p
        if x in inv_cache:
            return inv_cache[x]
        result = pow(x, p - 2, p)
        inv_cache[x] = result
        return result

    # Precompute 2^e mod p for exponents we'll encounter
    # Max exponent: roughly n^2/2

    active = active_divisors_at(odd_parts[0])
    div_index = {d: i for i, d in enumerate(active)}
    init_state = tuple(0 for _ in active)

    # dp[(total, state)] = weight mod p
    dp = {(0, init_state): 1}

    for ki, k in enumerate(odd_parts):
        if k > n:
            break

        divs_k = divisors_of[k]
        div_idx_k = [div_index[d] for d in divs_k if d in div_index]
        phi_k = [phi_cache[d] for d in divs_k if d in div_index]

        # Precompute k^m mod p and (m!)^{-1} mod p for m = 0..n//k
        max_m = n // k
        k_pow = [1] * (max_m + 1)  # k^m mod p
        for m in range(1, max_m + 1):
            k_pow[m] = k_pow[m-1] * k % p
        fact_inv = [1] * (max_m + 1)  # (m!)^{-1} mod p
        fact = 1
        for m in range(1, max_m + 1):
            fact = fact * m % p
            fact_inv[m] = mod_inv(fact)

        new_dp = {}
        for (total, state), weight in dp.items():
            max_m_here = (n - total) // k
            for m in range(0, max_m_here + 1):
                new_total = total + m * k

                # Compute Δt
                cross = 0
                for i in range(len(div_idx_k)):
                    cross += phi_k[i] * state[div_idx_k[i]]
                delta_t = m * cross + m * (m - 1) * k // 2 + m * (k - 1) // 2

                # Weight factor: 2^{delta_t} / (k^m · m!) mod p
                factor = pow(2, delta_t, p) * mod_inv(k_pow[m]) % p * fact_inv[m] % p

                # Update state
                new_state = list(state)
                for idx in div_idx_k:
                    new_state[idx] += m
                new_state = tuple(new_state)

                key = (new_total, new_state)
                if key in new_dp:
                    new_dp[key] = (new_dp[key] + weight * factor) % p
                else:
                    new_dp[key] = weight * factor % p

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
                        projected_dp[key] = (projected_dp[key] + val) % p
                    else:
                        projected_dp[key] = val
                dp = projected_dp
                active = new_active
                div_index = new_div_index

    # Sum over all states with total = n
    result = 0
    for (total, state), val in dp.items():
        if total == n:
            result = (result + val) % p
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


def a000568_crt(n, verbose=True):
    """Compute A000568(n) using CRT over multiple primes."""
    if n <= 1:
        return 1

    # Estimate number of bits in a(n)
    # Upper bound: 2^{C(n,2)} / n! ~ 2^{n(n-1)/2} / n!
    # log2(a(n)) ≈ n(n-1)/2 - n*log2(n) + n*log2(e) ≈ n(n-1)/2 - n*log2(n)
    import math
    est_bits = int(n * (n - 1) / 2 - n * math.log2(max(n, 2)) + n) + 100
    if est_bits < 100:
        est_bits = 100

    # Use ~30-bit primes so products fit in int64
    prime_bits = 30
    num_primes = est_bits // prime_bits + 5

    # Generate primes > n, starting near 2^30 for ~30-bit primes
    primes = primes_above(max(n, 2**30), num_primes)

    if verbose:
        print(f"  n={n}: using {num_primes} primes ({prime_bits}-bit), est {est_bits} bits")

    # Compute a(n) mod each prime
    residues = []
    for i, p in enumerate(primes):
        r = a000568_mod_p_dp(n, p)
        residues.append(r)

    # CRT reconstruction
    # Use incremental CRT for efficiency
    result = residues[0]
    modulus = primes[0]
    for i in range(1, len(residues)):
        r = residues[i]
        p = primes[i]
        # Combine: find x ≡ result (mod modulus) and x ≡ r (mod p)
        # x = result + modulus * t where t = (r - result) * modulus^{-1} mod p
        diff = (r - result % p) % p
        t = diff * pow(modulus, p - 2, p) % p
        result = result + modulus * t
        modulus = modulus * p

    # The result might need to be taken mod modulus (should already be < modulus)
    if result < 0:
        result += modulus

    return result


# Known values
known = {0: 1, 1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
         9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168}

if __name__ == "__main__":
    print("=" * 70)
    print("A000568 via DP + CRT")
    print("=" * 70)

    # Validate
    print("\nValidation:")
    for n in range(13):
        t0 = perf_counter()
        result = a000568_crt(n, verbose=False)
        dt = perf_counter() - t0
        expected = known.get(n, "?")
        match = "✓" if result == expected else f"✗ (got {result}, expected {expected})"
        print(f"  a({n:2d}) = {result:>15}  [{dt:.4f}s]  {match}")

    # Benchmark
    print("\n" + "=" * 70)
    print("SCALING TEST")
    print("=" * 70)

    for n in [20, 30, 40, 50, 60, 70, 76, 80, 90, 100]:
        t0 = perf_counter()
        result = a000568_crt(n, verbose=True)
        dt = perf_counter() - t0
        s = str(result)
        print(f"  a({n:3d}) = {s[:40]}{'...' if len(s) > 40 else '':3s} "
              f"({len(s):4d} digits)  [{dt:8.3f}s]")

    # Verify n=76 against OEIS
    print("\n" + "=" * 70)
    print("VERIFICATION vs OEIS")
    print("=" * 70)
    oeis_76 = 457153791459731763873714717642097124535855865100705012123154296782586066392599174177717111708812242945048119323356894628050571858925201215481427701704774498054845450499257140261101951670038378795297655176726665694457285914072690084154974004045609532979011190165711892081116449736469812420129338461848172656662018517999846546136461519522876668617723623760713964531731492213163219661503132797160477050776498667637823035846487741247681442428085385480808545797005707547708306606906039466675985395538562415540138435374979216235780263415436443048770796035424757246769390586508396300108786258541736828941774383334172760307820994535160154251942298642998619446487984483719749964075470321740294755571922446326714092505683852691999256324957470064758984015872
    result_76 = a000568_crt(76, verbose=False)
    print(f"  a(76) matches OEIS: {result_76 == oeis_76}")
