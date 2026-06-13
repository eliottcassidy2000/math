#!/usr/bin/env python3
"""
a000568_speedup.py — Fast computation of A000568 (nonisomorphic tournaments).

STATE OF THE ART: Davis formula (1954) sums over all partitions of n
into odd parts. Number of terms = p(n) ~ exp(pi*sqrt(2n/3)).

THIS SCRIPT explores multiple speedup strategies:

Strategy A: Direct Davis formula (baseline)
Strategy B: Polynomial convolution (process parts sequentially)
Strategy C: CRT + modular arithmetic
Strategy D: Dirichlet-accelerated gcd quadratic form

Author: opus-2026-03-07-S46f
"""

from math import gcd, factorial, isqrt
from functools import lru_cache
from collections import Counter
from time import perf_counter
import sys

# =====================================================================
# STRATEGY A: Direct Davis formula (baseline)
# =====================================================================

def euler_totient(n):
    """Euler's totient function φ(n)."""
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

def partitions_odd_parts(n):
    """Generate all partitions of n into odd parts (= all partitions of n by Euler)."""
    # We generate partitions of n into odd parts directly
    def _gen(remaining, max_part):
        if remaining == 0:
            yield []
            return
        # max_part must be odd
        if max_part % 2 == 0:
            max_part -= 1
        for k in range(min(remaining, max_part), 0, -2):  # odd parts only
            for rest in _gen(remaining - k, k):
                yield [k] + rest
    yield from _gen(n, n if n % 2 == 1 else n - 1)

def davis_t(partition):
    """Compute exponent t(λ) for a partition (given as sorted list of parts)."""
    # Count multiplicities
    mults = Counter(partition)
    parts = sorted(mults.keys())

    # t = (1/2) [Σ_{r,s} m_r m_s gcd(r,s) - Σ_r m_r]
    total_gcd_sum = 0
    for r in parts:
        for s in parts:
            total_gcd_sum += mults[r] * mults[s] * gcd(r, s)
    total_m = sum(mults.values())
    return (total_gcd_sum - total_m) // 2

def davis_z(partition):
    """Compute z_λ = Π_k k^{m_k} · m_k! for partition λ."""
    mults = Counter(partition)
    z = 1
    for k, m in mults.items():
        z *= k**m * factorial(m)
    return z

def a000568_davis(n):
    """Compute A000568(n) using direct Davis formula."""
    if n <= 1:
        return 1
    total = 0
    for partition in partitions_odd_parts(n):
        t = davis_t(partition)
        z = davis_z(partition)
        total += (1 << t) * factorial(n) // z
    return total // factorial(n)

# =====================================================================
# STRATEGY B: Polynomial convolution (sequential part processing)
# =====================================================================

def a000568_polyconv(n):
    """Compute A000568(n) via polynomial convolution over odd part sizes.

    Process parts k = 1, 3, 5, ... sequentially.
    Maintain poly[j] = weighted sum over partitions of j using parts ≤ k.

    For each new part size k, we need to account for:
    1. Self-interaction: m_k copies of k contribute m_k(m_k-1)/2 · k
       to the gcd sum (since gcd(k,k) = k)
    2. Cross-interaction: m_k copies of k interact with all previously
       processed parts via gcd(k, prev_parts)

    The cross-interaction prevents simple convolution. We handle it by
    maintaining auxiliary state: for each divisor d of future parts,
    track S_d = Σ_{j≤k, d|j} m_j.

    SIMPLIFIED VERSION: Just track the weighted sum directly.
    """
    if n <= 1:
        return 1

    # poly[j] = Σ over partitions of j (using allowed parts) of 2^{t(λ)} / z_λ
    # But t depends on cross-terms, so we can't do simple convolution...

    # ALTERNATIVE: Dynamic programming over the partition lattice.
    # State: (total_size, gcd_accumulator)
    # The gcd_accumulator needs to track enough info for future cross-terms.

    # For a practical speedup, we use a DP where:
    #   dp[j] = list of (weight, divisor_profile) pairs
    # where divisor_profile[d] = S_d for the partition so far.
    # This is exponential in the number of divisors, so let's limit it.

    # FALLBACK: Use the Dirichlet trick (Strategy D) within the Davis sum.
    # For now, implement a clean version that's correct and profile it.

    # Actually, the simplest speedup: avoid regenerating partitions,
    # and compute t(λ) using the Dirichlet trick.
    return a000568_dirichlet(n)

# =====================================================================
# STRATEGY C: CRT + modular arithmetic
# =====================================================================

def is_prime(n):
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def next_prime(n):
    n += 1
    while not is_prime(n):
        n += 1
    return n

def a000568_mod_p(n, p):
    """Compute A000568(n) mod p using Davis formula with modular arithmetic."""
    if n <= 1:
        return 1 % p

    n_fact = factorial(n) % p
    if n_fact == 0:
        # p ≤ n, need careful handling
        # Actually the formula a(n) = Σ 2^t / z_λ might have p in denominators
        # Skip for now — use primes > n
        return None

    n_fact_inv = pow(n_fact, p - 2, p)  # Fermat's little theorem

    total = 0
    for partition in partitions_odd_parts(n):
        t = davis_t(partition)
        z = davis_z(partition) % p
        if z == 0:
            continue
        z_inv = pow(z, p - 2, p)
        contrib = (pow(2, t, p) * (factorial(n) % p) % p * z_inv) % p
        total = (total + contrib) % p

    return (total * n_fact_inv) % p

def a000568_crt(n, num_primes=None):
    """Compute A000568(n) using CRT over multiple primes.

    Choose primes p > n to avoid factorial issues.
    """
    if n <= 1:
        return 1

    # Estimate number of digits in a(n): roughly n^2 * log(2) / (2 * log(10))
    est_bits = n * (n - 1) // 2  # upper bound: 2^{C(n,2)}
    if num_primes is None:
        num_primes = est_bits // 30 + 10  # 30-bit primes

    # Collect primes > n
    p = max(n + 1, 100)
    primes = []
    while len(primes) < num_primes:
        p = next_prime(p)
        primes.append(p)

    # Compute a(n) mod each prime
    residues = []
    for pr in primes:
        r = a000568_mod_p(n, pr)
        if r is not None:
            residues.append((r, pr))

    # CRT reconstruction
    result = 0
    M = 1
    for r, p in residues:
        M *= p

    for r, p in residues:
        Mi = M // p
        yi = pow(Mi, p - 2, p)  # Mi^{-1} mod p
        result = (result + r * Mi * yi) % M

    # The actual a(n) might be less than M
    return result

# =====================================================================
# STRATEGY D: Dirichlet-accelerated gcd computation
# =====================================================================

def a000568_dirichlet(n):
    """Davis formula with Dirichlet-accelerated gcd quadratic form.

    Uses: gcd(r,s) = Σ_{d|gcd(r,s)} φ(d)
    So: Σ_{r,s} m_r m_s gcd(r,s) = Σ_d φ(d) (Σ_{d|r} m_r)²

    For each partition, this replaces O(k²) gcd calls with O(n) divisor sums.
    """
    if n <= 1:
        return 1

    # Precompute Euler totients
    phi = [0] * (n + 1)
    for i in range(1, n + 1):
        phi[i] = euler_totient(i)

    total = 0
    for partition in partitions_odd_parts(n):
        mults = Counter(partition)
        parts = sorted(mults.keys())

        # Compute S_d = Σ_{d|r} m_r for each d
        S = {}
        for d in range(1, n + 1):
            sd = sum(mults.get(r, 0) for r in range(d, n + 1, d) if r % 2 == 1)
            if sd > 0:
                S[d] = sd

        # Compute Σ_{r,s} m_r m_s gcd(r,s) = Σ_d φ(d) S_d²
        gcd_sum = sum(phi[d] * S[d]**2 for d, sd in S.items() if d in S)

        total_m = sum(mults.values())
        t = (gcd_sum - total_m) // 2

        z = 1
        for k, m in mults.items():
            z *= k**m * factorial(m)

        total += (1 << t) * factorial(n) // z

    return total // factorial(n)

# =====================================================================
# STRATEGY E: DP over part sizes (avoid explicit partition enumeration)
# =====================================================================

def a000568_dp(n):
    """Compute A000568(n) by DP over odd part sizes.

    Instead of enumerating partitions, process parts k=1,3,5,...
    and accumulate the weighted sum.

    dp[j] = rational number = partial sum of 2^{t(λ)} / z_λ
    over partitions of j using only parts from the already-processed set.

    PROBLEM: The cross-terms in t(λ) depend on which specific multiplicities
    were chosen for previous parts, not just the total j. So we need to
    track more state.

    SOLUTION: Track the "divisor profile" — for each odd divisor d ≤ n,
    what is S_d (the sum of multiplicities of parts divisible by d).
    Then t(λ) = (1/2)[Σ_d φ(d) S_d² - Σ_r m_r].

    When adding m copies of part k, S_d increases by m for each d | k.
    """
    if n <= 1:
        return 1

    from fractions import Fraction

    # Precompute
    phi = [0] * (n + 1)
    for i in range(1, n + 1):
        phi[i] = euler_totient(i)

    odd_parts = list(range(1, n + 1, 2))

    # For each odd part k, find its odd divisors
    divisors_of = {}
    for k in odd_parts:
        divisors_of[k] = [d for d in range(1, k + 1, 2) if k % d == 0]

    # State: (total_used, tuple of S_d values for relevant divisors)
    # This is too large in general. Let's limit to small n for validation.

    # SIMPLER DP: process parts and track the gcd-relevant state.
    # For n ≤ 30 or so, we can use a dict-based DP.

    # State = (total_size, frozenset of (d, S_d) pairs)
    # Start: total=0, all S_d=0
    # For each part size k = 1, 3, 5, ...:
    #   For each multiplicity m = 0, 1, ..., floor(n/k):
    #     New total = old_total + m*k
    #     For each d | k: S_d += m
    #     Self-interaction of this part: m(m-1)/2 * k (within the gcd sum)
    #     Cross-interaction: m * Σ_{d|k} φ(d) * (old_S_d) ... wait, more complex

    # The exponent contribution from adding m copies of part k:
    # Δt = (1/2) [Σ_d φ(d) (S_d + m·[d|k])² - Σ_d φ(d) S_d² - m]
    #     = (1/2) [Σ_{d|k} φ(d) (2·m·S_d + m²) - m]
    #     = m · Σ_{d|k} φ(d) · S_d + (m²/2) · Σ_{d|k} φ(d) - m/2
    #     = m · Σ_{d|k} φ(d) · S_d + m(m-1)/2 · k + m·(k-1)/2

    # Note: Σ_{d|k} φ(d) = k (standard identity). So:
    # Δt = m · Σ_{d|k} φ(d) · S_d + m(m-1)/2 · k + m(k-1)/2

    # The weight contribution: 2^{Δt} / (k^m · m!)

    # So the DP is:
    # For each state (total, S_profile), weight w:
    #   For each multiplicity m of part k:
    #     Δt = m · Σ_{d|k} φ(d) · S_d + m(m-1)k/2 + m(k-1)/2
    #     new_weight = w * 2^{Δt} / (k^m · m!)
    #     Update S_d for d|k: S_d += m

    # The S_profile only needs to track S_d for ODD d that divide some
    # future part. As we process parts in order, fewer divisors matter.

    # For practical n, let's just implement this with a dictionary.

    # Relevant odd divisors (those that appear as divisors of parts 1..n)
    all_odd_divs = sorted(set(d for k in odd_parts for d in divisors_of[k]))

    # State: tuple(S_d for d in all_odd_divs)
    # Initial state: all zeros
    init_state = tuple(0 for _ in all_odd_divs)
    div_index = {d: i for i, d in enumerate(all_odd_divs)}

    # dp[state] = (total_size, accumulated_weight_as_Fraction)
    # Actually: dp[(total, state)] = accumulated weight
    # This is too big. Let's compress.

    # dp is a dict: (total, state_tuple) -> weight (Fraction)
    dp = {(0, init_state): Fraction(1)}

    for k in odd_parts:
        if k > n:
            break
        new_dp = {}
        divs_k = divisors_of[k]
        div_indices_k = [div_index[d] for d in divs_k]

        for (total, state), weight in dp.items():
            max_m = (n - total) // k
            for m in range(0, max_m + 1):
                new_total = total + m * k

                # Compute Δt
                cross = sum(phi[divs_k[i]] * state[div_indices_k[i]] for i in range(len(divs_k)))
                delta_t = m * cross + m * (m - 1) * k // 2 + m * (k - 1) // 2

                # Weight factor
                factor = Fraction(1 << delta_t, k**m * factorial(m))

                # New state: S_d += m for d | k
                new_state = list(state)
                for idx in div_indices_k:
                    new_state[idx] += m
                new_state = tuple(new_state)

                key = (new_total, new_state)
                if key in new_dp:
                    new_dp[key] = new_dp[key] + weight * factor
                else:
                    new_dp[key] = weight * factor

        dp = new_dp

    # Sum over all states with total = n
    result = sum(w for (total, state), w in dp.items() if total == n)
    return int(result)

# =====================================================================
# VALIDATION
# =====================================================================

known = {0: 1, 1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
         9: 191536, 10: 9733056}

def validate(method, name, max_n=10):
    print(f"\n{'='*60}")
    print(f"Validating: {name}")
    print(f"{'='*60}")
    for n in range(max_n + 1):
        t0 = perf_counter()
        result = method(n)
        dt = perf_counter() - t0
        expected = known.get(n)
        match = "✓" if result == expected else f"✗ (expected {expected})"
        print(f"  a({n:2d}) = {result:>15d}  [{dt:.4f}s]  {match}")

def benchmark(methods, n_range):
    print(f"\n{'='*60}")
    print(f"Benchmark: n = {min(n_range)}..{max(n_range)}")
    print(f"{'='*60}")
    for name, method in methods:
        times = []
        for n in n_range:
            t0 = perf_counter()
            result = method(n)
            dt = perf_counter() - t0
            times.append(dt)
            if n == max(n_range):
                print(f"  {name:25s}: a({n}) = {result}, time = {dt:.4f}s")
        total = sum(times)
        print(f"  {name:25s}: total time = {total:.4f}s")

if __name__ == "__main__":
    # Validate all methods
    validate(a000568_davis, "Davis (baseline)", max_n=8)
    validate(a000568_dirichlet, "Dirichlet-accelerated", max_n=8)
    validate(a000568_dp, "DP over part sizes", max_n=8)

    # Benchmark
    methods = [
        ("Davis (baseline)", a000568_davis),
        ("Dirichlet", a000568_dirichlet),
        ("DP", a000568_dp),
    ]
    benchmark(methods, range(6, 11))

    # Try larger n with DP
    print(f"\n{'='*60}")
    print(f"Pushing to larger n with DP method")
    print(f"{'='*60}")
    for n in [12, 14, 16, 18, 20]:
        t0 = perf_counter()
        result = a000568_dp(n)
        dt = perf_counter() - t0
        expected = known.get(n, "?")
        print(f"  a({n:2d}) = {result:>20d}  [{dt:.4f}s]")
