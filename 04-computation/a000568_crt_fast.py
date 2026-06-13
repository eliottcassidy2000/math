"""
a000568_crt_fast.py — opus-2026-04-04-S22

FAST A000568 via CRT: compute a(n) mod p for many primes p using
native 64-bit arithmetic, then reconstruct via Chinese Remainder Theorem.

THE BINARY DECISION PERSPECTIVE:
A tournament = C(n,2) binary choices. Burnside compresses this into
a sum over odd-part partitions of n. CRT compresses the bignum
arithmetic into independent modular computations.

Each prime p gives a(n) mod p using the SAME partition enumeration
but with O(1) modular ops instead of O(d) GMP ops per partition.

For n=200: a(n) has ~5400 digits = ~18000 bits.
Need ~18000/62 ≈ 290 primes near 2^62.
Each mod-p computation: ~seconds.
Total: ~290 × seconds = ~15 minutes (vs 11 minutes for GMP).

BUT: with 8-way parallelism: ~2 minutes. And each prime is INDEPENDENT.
"""
import time
from math import gcd, factorial, log2, ceil
from functools import reduce

def compute_a000568_mod_p(n, p):
    """
    Compute a(n) mod p using the Burnside formula restricted to odd-part partitions.

    a(n) = (1/n!) Σ_{λ: all odd parts} (n!/z_λ) * 2^{e(λ)} mod p

    We compute: n! * a(n) = Σ (n!/z_λ) * 2^{e(λ)} mod p
    Then: a(n) = [n! * a(n)] * inverse(n!, p) mod p

    The edge count e(λ) for cycle type λ = (l₁^{m₁}, l₃^{m₃}, ...):
    e(λ) = Σ_i m_i * (l_i - 1) / 2  +  Σ_i C(m_i, 2) * l_i  +  Σ_{i<j} m_i * m_j * gcd(l_i, l_j)
    """
    if n <= 2:
        return 1 % p

    # Precompute
    # Odd parts up to n
    odd_parts = list(range(1, n+1, 2))  # 1, 3, 5, ...
    num_odd = len(odd_parts)

    # GCD table
    gcd_tab = {}
    for i, a in enumerate(odd_parts):
        for j, b in enumerate(odd_parts):
            if j >= i:
                gcd_tab[(i,j)] = gcd(a, b)
                gcd_tab[(j,i)] = gcd(a, b)

    # Precompute n! mod p and inverse
    nfact_mod = 1
    for k in range(2, n+1):
        nfact_mod = (nfact_mod * k) % p

    # Inverse of n! mod p
    nfact_inv = pow(nfact_mod, p-2, p)  # Fermat's little theorem

    # Precompute powers of 2 mod p
    max_edges = n * (n-1) // 2
    pow2 = [1] * (max_edges + 1)
    for i in range(1, max_edges + 1):
        pow2[i] = (pow2[i-1] * 2) % p

    # Precompute k^m * m! mod p for small m and each odd k
    # z_λ = Π_i (l_i^{m_i} * m_i!)
    # n!/z_λ = n! / Π_i (l_i^{m_i} * m_i!)

    # We'll compute (n!/z_λ) mod p incrementally during enumeration

    # Inverse table for odd parts: inv_k[i] = inverse of parts[i] mod p
    inv_k = {}
    for i, k in enumerate(odd_parts):
        if k % p != 0:
            inv_k[i] = pow(k, p-2, p)
        else:
            inv_k[i] = 0  # Handle p | k case separately

    # Enumerate odd-part partitions via recursive backtracking
    total = [0]  # accumulator mod p

    def enum_parts(remaining, part_idx, edges, z_inv_mod):
        """
        remaining: sum left to allocate
        part_idx: current part index (descending order)
        edges: current edge count
        z_inv_mod: current value of (1/z_partial) mod p
        """
        if remaining == 0:
            # Contribution: (n!/z) * 2^edges mod p
            contrib = (nfact_mod * z_inv_mod % p) * pow2[edges] % p
            total[0] = (total[0] + contrib) % p
            return

        if part_idx < 0:
            return

        k = odd_parts[part_idx]
        if k > remaining:
            enumerate(remaining, part_idx - 1, edges, z_inv_mod)
            return

        # Try multiplicity m = 0
        enumerate(remaining, part_idx - 1, edges, z_inv_mod)

        # Try m = 1, 2, ... up to remaining // k
        max_m = remaining // k

        # Incremental computation of:
        # z_contribution: k^m * m! (divide from z)
        # edge_contribution from this part

        cum_k_pow = 1  # k^m mod p
        cum_fact = 1   # m! mod p
        cum_self_edges = 0  # Σ for this part: m*(k-1)/2 + C(m,2)*k

        for m in range(1, max_m + 1):
            cum_k_pow = cum_k_pow * k % p
            cum_fact = cum_fact * m % p

            # Self-edge contribution
            # m*(k-1)/2 from single-cycle arcs
            # + C(m,2)*k from same-size pairs
            self_single = m * (k - 1) // 2  # integer, exact
            self_pairs = m * (m - 1) // 2 * k  # integer, exact

            # Cross-edge contribution with previously placed parts
            # This is tracked implicitly... we need to track placed parts.
            # PROBLEM: the recursive enumeration doesn't track cross-edges incrementally.

            # ALTERNATIVE: compute edges from scratch at each leaf.
            # This is O(d^2) per partition where d = number of distinct parts.
            # For the DP version, we'd track divisor sums.

            # For now: just track the partition and compute edges at the leaf.
            pass

        # Actually, let me use a simpler approach: enumerate partitions,
        # compute edges at each leaf. Less efficient but correct.

    # SIMPLER: enumerate partitions explicitly, compute at each leaf
    def enum_simple(remaining, max_part_idx, partition):
        """partition = list of (part_idx, multiplicity) pairs"""
        if remaining == 0:
            # Compute e(λ) and z(λ)
            edges = 0
            z_mod = 1  # z_λ mod p

            parts_list = [(odd_parts[pi], m) for pi, m in partition]

            for i, (ki, mi) in enumerate(parts_list):
                # Self-cycle: mi * (ki-1) / 2
                edges += mi * (ki - 1) // 2
                # Same-size pairs: C(mi, 2) * ki
                edges += mi * (mi - 1) // 2 * ki
                # z contribution: ki^mi * mi!
                z_mod = z_mod * pow(ki, mi, p) % p
                for f in range(1, mi + 1):
                    z_mod = z_mod * f % p

                # Cross pairs
                for j in range(i+1, len(parts_list)):
                    kj, mj = parts_list[j]
                    edges += mi * mj * gcd(ki, kj)

            # Contribution: (n!/z) * 2^edges mod p
            if z_mod % p == 0:
                return  # p divides z, skip (very rare)
            z_inv = pow(z_mod, p-2, p)
            contrib = nfact_mod * z_inv % p * pow2[min(edges, max_edges)] % p
            total[0] = (total[0] + contrib) % p
            return

        for pi in range(max_part_idx, -1, -1):
            k = odd_parts[pi]
            if k > remaining:
                continue
            for m in range(1, remaining // k + 1):
                if m * k > remaining:
                    break
                enum_simple(remaining - m * k, pi - 1, partition + [(pi, m)])

    enum_simple(n, num_odd - 1, [])

    # Result: total = n! * a(n) mod p. Multiply by (n!)^{-1}
    result = total[0] * nfact_inv % p
    return result

def crt_combine(residues, primes):
    """Chinese Remainder Theorem: reconstruct integer from residues."""
    # Using successive combination
    result = residues[0]
    modulus = primes[0]
    for i in range(1, len(residues)):
        r = residues[i]
        p = primes[i]
        # Find x such that x ≡ result (mod modulus) and x ≡ r (mod p)
        # x = result + modulus * t where t ≡ (r - result) * modulus^{-1} (mod p)
        diff = (r - result % p) % p
        mod_inv = pow(modulus % p, p - 2, p)  # modulus^{-1} mod p
        t = diff * mod_inv % p
        result = result + modulus * t
        modulus = modulus * p
    return result

def compute_a000568_fast(n):
    """Compute a(n) exactly using CRT."""
    t0 = time.time()

    # Estimate number of bits
    # a(n) ~ 2^{C(n,2)} / n! → log2(a) ~ C(n,2) - log2(n!)
    est_log2 = n*(n-1)//2 - sum(log2(k) for k in range(1, n+1))
    est_bits = int(est_log2) + 10  # safety margin

    # Choose primes
    # Use primes near 2^62 for maximum range
    # For small n, fewer primes needed
    prime_bits = 30  # use ~30-bit primes for safety
    num_primes = max(2, est_bits // prime_bits + 2)

    # Generate primes
    from sympy import nextprime
    primes = []
    p = 2**prime_bits
    for _ in range(num_primes):
        p = nextprime(p)
        primes.append(p)

    print(f"  n={n}: est {est_bits} bits, using {num_primes} primes of ~{prime_bits} bits")

    # Compute residues
    residues = []
    for i, p in enumerate(primes):
        r = compute_a000568_mod_p(n, p)
        residues.append(r)

    t_compute = time.time() - t0

    # CRT reconstruction
    result = crt_combine(residues, primes)

    t_total = time.time() - t0
    print(f"  Time: compute={t_compute:.3f}s, total={t_total:.3f}s")

    return result

# ==============================================================
# BENCHMARK
# ==============================================================
if __name__ == "__main__":
    print("A000568 FAST CRT COMPUTATION")
    print("=" * 60)

    # Known values for verification
    known = {
        3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
        9: 191536, 10: 9733056
    }

    for n in [3, 4, 5, 6, 7, 8, 9, 10]:
        result = compute_a000568_fast(n)
        expected = known[n]
        match = "✓" if result == expected else f"✗ (expected {expected})"
        print(f"  a({n}) = {result} {match}")

    # Benchmark at larger n
    print(f"\nBenchmark:")
    for n in [15, 20]:
        t0 = time.time()
        result = compute_a000568_fast(n)
        elapsed = time.time() - t0
        print(f"  a({n}) = {result} ({elapsed:.3f}s)")
