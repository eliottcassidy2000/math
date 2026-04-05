"""
a000568_fiber_recursion.py — opus-2026-04-04-S23

A COMPLETELY DIFFERENT APPROACH to A000568: use the fiber bundle recursion.

From S12 and S14:
  a(n) = (1/n) Σ_{classes C at n-1} orbits(C)
  orbits(C) = (1/|Aut(C)|) Σ_{α∈Aut(C)} 2^{cycles(α on n-1 vertices)}

BUT this overcounts (merging ratio ≈ n, not exact).

HOWEVER: there's an EXACT formula using Burnside at the fiber level:

  a(n) = (1/n!) Σ_{σ∈S_n, all cycles odd} 2^{e(σ)}

The standard Burnside. Can we compute this WITHOUT enumerating partitions?

NEW IDEA: "Binary Decision Burnside"
Instead of grouping by partition type, process VERTICES one at a time.

For each vertex v added to the tournament:
  - v creates (v-1) new arc decisions
  - The Burnside contribution of σ that maps v to position p depends
    on the cycle structure at that point

This is a SEQUENTIAL Burnside computation: build up the permutation
group action vertex-by-vertex, tracking partial contributions.

ALTERNATIVE: Direct DP on cycle types
Instead of enumerating partitions, use a DP that builds up the cycle type
by choosing how each "remaining size" is allocated.

State: (remaining_to_allocate, list_of_divisor_sums)
Transition: add one more part of odd size k

The key speedup over raw enumeration: we can BATCH partitions that
differ only in the order of their parts (since we enumerate in
decreasing order, this is already done).

PRACTICAL: Let me implement the "generating function" approach.
The Burnside sum is:
  n! * a(n) = Σ_{λ: all odd} (n!/z_λ) * 2^{e(λ)}
           = coefficient of x^n in: Π_{k odd} Σ_{m≥0} (x^{mk} * 2^{e_self(k,m)} / (k^m * m!))
             × cross-interaction terms

The self-contribution for part k with multiplicity m is:
  e_self(k,m) = m*(k-1)/2 + C(m,2)*k

The cross-contribution between parts (k₁,m₁) and (k₂,m₂) is:
  e_cross = m₁*m₂*gcd(k₁,k₂)

This is NOT a simple product because the cross terms couple different parts.
But we CAN absorb the cross terms into a DP by tracking divisor sums.

STATE: for each divisor d of parts seen so far, track S_d = Σ_{k placed: d|k} m_k

Then: e_cross(new part k with mult m) = m * Σ_{d|k} φ(d) * S_d

Wait, that's not quite right. The cross-edge formula is:
  Σ_{i<j} m_i * m_j * gcd(k_i, k_j)

If we define S = Σ_i m_i * k_i (total weight) and T = Σ_i m_i * Σ_{d|k_i} d,
then... this doesn't simplify easily.

Actually, the standard trick: gcd(a,b) = Σ_{d|gcd(a,b)} φ(d) = Σ_{d|a, d|b} φ(d).

So: Σ_{i<j} m_i m_j gcd(k_i,k_j) = (1/2)[( Σ_d φ(d) (Σ_{k_i: d|k_i} m_i)² ) - Σ_d φ(d) Σ_{k_i: d|k_i} m_i²]
  = (1/2) Σ_d φ(d) [S_d² - Σ_{i: d|k_i} m_i²]

where S_d = Σ_{i: d|k_i} m_i.

This is the key! The cross-edge contribution is a QUADRATIC FORM in the
divisor sums S_d. And we can track S_d incrementally in a DP.

DP STATE: (remaining, S_1, S_3, S_5, ..., S_{max_odd_divisor})
Each odd divisor d of any part k ≤ n needs its own counter.

Odd divisors of odd parts ≤ n: all odd numbers ≤ n.
So the state space has dimension = #{odd numbers ≤ n} = ceil(n/2).

For n=100: ceil(100/2) = 50 dimensions. State space is HUGE.

But most S_d values are small (bounded by n/d), so the actual state
count is much smaller. Let me estimate and test.

Actually, the STANDARD approach is:
  Process parts in decreasing order (k = n, n-2, ..., 3, 1).
  For each part k: choose multiplicity m = 0, 1, 2, ..., floor(remaining/k).
  Update: remaining -= m*k.
  Update: for each d|k: S_d += m.
  Accumulate: weight *= 2^{e_self(k,m)} / (k^m * m!)
  Cross contribution: += m * Σ_{d|k} φ(d) * S_d_before

This gives an O(p(n) * #divisors_per_part) algorithm, same as
raw enumeration but without explicit partition generation.

For MODULAR computation (CRT): this is perfect because all arithmetic
is mod p using native 64-bit integers.

Let me implement and benchmark.
"""
from math import gcd
from sympy import totient, nextprime
import time

def compute_a000568_dp_mod_p(n, p):
    """
    Compute a(n) mod p using DP over odd-part partitions
    with divisor-sum tracking for cross-edges.
    """
    if n <= 2: return 1 % p

    # Odd parts in decreasing order
    odd_parts = list(range(n if n % 2 == 1 else n-1, 0, -2))

    # Odd divisors needed
    odd_divs = list(range(1, n+1, 2))
    div_idx = {d: i for i, d in enumerate(odd_divs)}
    ndivs = len(odd_divs)

    # Precompute: for each odd part k, its odd divisors
    part_divs = {}
    for k in odd_parts:
        part_divs[k] = [d for d in odd_divs if k % d == 0]

    # Euler totient for odd numbers (precompute)
    phi = {}
    for d in odd_divs:
        phi[d] = int(totient(d))

    # n! mod p
    nfact = 1
    for k in range(2, n+1): nfact = nfact * k % p

    # Powers of 2 mod p
    max_e = n * (n-1) // 2
    pow2 = [1] * (max_e + 1)
    for i in range(1, max_e + 1): pow2[i] = pow2[i-1] * 2 % p

    # Recursive DP: process parts from largest to smallest
    total = [0]

    def dp(part_index, remaining, S_divs, edges_so_far, z_inv_so_far):
        """
        part_index: which part we're considering (index into odd_parts)
        remaining: sum left to allocate
        S_divs: list of S_d values (divisor sums of placed parts)
        edges_so_far: accumulated edge count
        z_inv_so_far: accumulated 1/z contribution mod p
        """
        if remaining == 0:
            contrib = nfact * z_inv_so_far % p * pow2[min(edges_so_far, max_e)] % p
            total[0] = (total[0] + contrib) % p
            return

        if part_index >= len(odd_parts):
            return

        k = odd_parts[part_index]
        if k > remaining:
            dp(part_index + 1, remaining, S_divs, edges_so_far, z_inv_so_far)
            return

        # m = 0: skip this part
        dp(part_index + 1, remaining, S_divs, edges_so_far, z_inv_so_far)

        # m = 1, 2, ... up to remaining // k
        max_m = remaining // k
        k_pow_inv = pow(k, p-2, p)  # k^{-1} mod p
        cum_z_factor = 1  # (k^m * m!)^{-1} accumulator

        for m in range(1, max_m + 1):
            # Update z: divide by k * m
            cum_z_factor = cum_z_factor * k_pow_inv % p
            cum_z_factor = cum_z_factor * pow(m, p-2, p) % p

            # Self-edge contribution: m*(k-1)/2 + C(m,2)*k
            e_self = m * (k - 1) // 2 + m * (m - 1) // 2 * k

            # Cross-edge contribution: m * Σ_{d|k} φ(d) * S_d
            e_cross = 0
            for d in part_divs[k]:
                e_cross += m * phi[d] * S_divs[div_idx[d]]

            new_edges = edges_so_far + e_self + e_cross

            # Update S_divs
            new_S = S_divs.copy()
            for d in part_divs[k]:
                new_S[div_idx[d]] += m

            dp(part_index + 1, remaining - m * k, new_S,
               new_edges, z_inv_so_far * cum_z_factor % p)

    initial_S = [0] * ndivs
    dp(0, n, initial_S, 0, 1)

    nfact_inv = pow(nfact, p-2, p)
    return total[0] * nfact_inv % p

# Test
if __name__ == "__main__":
    known = {3:2, 4:4, 5:12, 6:56, 7:456, 8:6880, 9:191536, 10:9733056}

    print("DP-CRT APPROACH TO A000568")
    print("=" * 60)

    p = nextprime(2**30)
    for n in [3, 4, 5, 6, 7, 8, 9, 10]:
        t0 = time.time()
        r = compute_a000568_dp_mod_p(n, p)
        elapsed = time.time() - t0
        expected = known[n] % p
        match = "✓" if r == expected else f"✗ (got {r}, expected {expected})"
        print(f"  a({n}) mod p = {r} {match} ({elapsed:.4f}s)")

    # Benchmark at n=30, 50
    for n in [20, 30, 40]:
        t0 = time.time()
        r = compute_a000568_dp_mod_p(n, p)
        elapsed = time.time() - t0
        print(f"  a({n}) mod p = {r} ({elapsed:.3f}s)")
