#!/usr/bin/env python3
"""
a000568_hybrid_crt.py — Hybrid multi-process CRT for A000568.

For each prime p, runs the C enum binary (a000568_c_enum) to compute a(n) mod p.
All primes run in parallel across CPU cores. CRT reconstruction in Python.

Designed for the crossover point where CRT beats GMP enum
(very large n where GMP numbers become enormous).

Author: opus-2026-03-08-S48
"""

import os
import sys
import math
import subprocess
from time import perf_counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from sympy import nextprime  # for generating primes

_dir = os.path.dirname(os.path.abspath(__file__))
_binary = os.path.join(_dir, "a000568_c_enum")


def _compute_one(args):
    """Compute a(n) mod p using C enumeration."""
    n, p = args
    result = subprocess.check_output([_binary, str(n), str(p)],
                                     timeout=7200)
    return int(result.strip())


def estimate_bits(n):
    """Estimate number of bits in a(n).
    a(n) ~ 2^{C(n,2)} / n!, so log2(a(n)) ~ C(n,2) - log2(n!)."""
    bits = n * (n - 1) / 2 - sum(math.log2(k) for k in range(2, n + 1))
    return int(bits) + 100  # safety margin


def generate_primes(start, count):
    """Generate count primes starting above start."""
    primes = []
    p = start
    for _ in range(count):
        p = nextprime(p)
        primes.append(p)
    return primes


def a000568_hybrid_crt(n, num_workers=None, verbose=True):
    """Compute A000568(n) using parallel C enum + CRT."""
    if n <= 1:
        return 1

    est_bits = estimate_bits(n)
    prime_bits = 30  # each prime ~30 bits
    num_primes = est_bits // prime_bits + 5

    # Use primes above max(n, 2^30) to avoid issues with small primes
    prime_start = max(n, 2**30)
    primes = generate_primes(prime_start, num_primes)

    if num_workers is None:
        num_workers = os.cpu_count() or 4

    if verbose:
        print(f"  n={n}: {num_primes} primes, {num_workers} workers, est {est_bits} bits",
              flush=True)

    t0 = perf_counter()
    args_list = [(n, p) for p in primes]

    residues = [None] * num_primes
    completed = 0

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        future_to_idx = {executor.submit(_compute_one, args_list[i]): i
                         for i in range(num_primes)}
        for future in as_completed(future_to_idx):
            idx = future_to_idx[future]
            residues[idx] = future.result()
            completed += 1
            if verbose and completed % 10 == 0:
                elapsed = perf_counter() - t0
                print(f"    {completed}/{num_primes} primes done ({elapsed:.1f}s)",
                      flush=True)

    t_par = perf_counter() - t0
    if verbose:
        print(f"  Parallel enum: {t_par:.1f}s ({t_par/num_primes:.3f}s/prime avg)",
              flush=True)

    # CRT reconstruction
    t_crt = perf_counter()
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

    t_crt_elapsed = perf_counter() - t_crt
    if verbose:
        print(f"  CRT reconstruction: {t_crt_elapsed:.1f}s", flush=True)

    return result


# Known values for validation
known = {0: 1, 1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
         9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168}


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 a000568_hybrid_crt.py <n> [num_workers]")
        sys.exit(1)

    n = int(sys.argv[1])
    nw = int(sys.argv[2]) if len(sys.argv) > 2 else None

    t0 = perf_counter()
    result = a000568_hybrid_crt(n, num_workers=nw)
    dt = perf_counter() - t0

    s = str(result)
    print(f"\na({n}) = {s}")
    print(f"({len(s)} digits, {dt:.1f}s total)")

    if n in known:
        print(f"Validation: {'OK' if result == known[n] else 'FAIL'}")
