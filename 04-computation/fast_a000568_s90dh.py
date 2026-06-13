#!/usr/bin/env python3
"""
fast_a000568_s90dh.py — Fast formula for A000568 (non-isomorphic tournaments)
opus-2026-03-16-S90dh

T(n) = (1/n!) * sum over partitions lambda of n:
       |C_lambda| * 2^{c(lambda)}

where c(lambda) = sum_i floor(l_i/2) + sum_{i<j} gcd(l_i, l_j)
and |C_lambda| = n! / (prod m_i! * prod l_i)
"""

from math import gcd, factorial, comb, log2
from functools import lru_cache
from collections import Counter
import time

def partitions(n, max_part=None):
    """Generate all partitions of n."""
    if max_part is None:
        max_part = n
    if n == 0:
        yield ()
        return
    for first in range(min(n, max_part), 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest

def c_lambda(parts):
    """Number of orbits of a permutation with cycle type `parts` on C(n,2) pairs.
    c(lambda) = sum_i floor(l_i / 2) + sum_{i<j} gcd(l_i, l_j)
    """
    s = sum(p // 2 for p in parts)
    for i in range(len(parts)):
        for j in range(i+1, len(parts)):
            s += gcd(parts[i], parts[j])
    return s

def conj_class_size(n, parts):
    """Size of the conjugacy class in S_n with cycle type `parts`.
    = n! / (prod_{cycle length l} (l^{m_l} * m_l!))
    where m_l = multiplicity of l in the partition.
    """
    counts = Counter(parts)
    denom = 1
    for l, m in counts.items():
        denom *= (l ** m) * factorial(m)
    return factorial(n) // denom

def T(n):
    """Number of non-isomorphic tournaments on n vertices (A000568)."""
    if n <= 1:
        return 1
    total = 0
    for parts in partitions(n):
        c = c_lambda(parts)
        size = conj_class_size(n, parts)
        total += size * (2 ** c)
    return total // factorial(n)

# Self-complementary count: include the transpose involution
def T_self_dual(n):
    """Number of self-complementary tournament classes on n vertices."""
    # = (1/(2*n!)) * [sum_sigma 2^c(sigma) + sum_sigma 2^c'(sigma)]
    # where c'(sigma) counts orbits of sigma composed with pair-complement.
    # For the complement: c'(lambda) = number of orbits of sigma on pairs
    # where the action includes flipping all pairs.
    # This is more complex; for now just use the basic count.
    pass

# Compute and verify
print("A000568: Non-isomorphic tournaments on n vertices")
print(f"{'n':>4} | {'T(n)':>15} | {'time (ms)':>10} | {'log2(T)':>8} | {'verify':>8}")
print("-" * 55)

known = {0:1, 1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456, 8:6880}

for n in range(0, 21):
    t0 = time.time()
    tn = T(n)
    elapsed = (time.time() - t0) * 1000
    verify = ""
    if n in known:
        verify = "OK" if tn == known[n] else f"FAIL (expected {known[n]})"
    lg = log2(tn) if tn > 0 else 0
    print(f"{n:4d} | {tn:15d} | {elapsed:10.1f} | {lg:8.1f} | {verify}")

# Extended computation
print(f"\nExtended values:")
for n in range(21, 31):
    t0 = time.time()
    tn = T(n)
    elapsed = (time.time() - t0) * 1000
    print(f"  T({n}) = {tn}  ({elapsed:.0f} ms)")

# The formula breakdown for small n
print(f"\nFormula breakdown:")
for n in [3, 4, 5, 6]:
    print(f"\n  n = {n}:")
    total = 0
    for parts in partitions(n):
        c = c_lambda(parts)
        size = conj_class_size(n, parts)
        contrib = size * (2 ** c)
        total += contrib
        if contrib > 0:
            print(f"    lambda={parts}: |C|={size}, c={c}, 2^c={2**c}, contrib={contrib}")
    print(f"    Total = {total}, /n! = {total}/{factorial(n)} = {total//factorial(n)}")

# Speed comparison: how fast can we go?
print(f"\nSpeed test:")
for n in [10, 15, 20, 25, 30, 35, 40]:
    t0 = time.time()
    tn = T(n)
    elapsed = (time.time() - t0) * 1000
    nparts = sum(1 for _ in partitions(n))
    print(f"  T({n:2d}) = {tn:>30d}  ({elapsed:8.0f} ms, {nparts} partitions)")

print(f"""

THE FORMULA:

  T(n) = (1/n!) * sum_{{lambda partition of n}} |C_lambda| * 2^{{c(lambda)}}

  where:
    c(lambda) = sum_i floor(l_i/2) + sum_{{i<j}} gcd(l_i, l_j)
    |C_lambda| = n! / prod_l (l^m_l * m_l!)

  This runs in O(p(n) * n^2) time where p(n) = number of partitions.
  p(20) = 627, p(30) = 5604, p(40) = 37338.
  For n=40: computable in seconds.
  For n=100: p(100) ~ 190 million. Needs optimization.

  The FAST PATH: precompute gcd table, vectorize the partition sum.

THE EIGENVALUE CONNECTION:

  The transition eigenvalues of the flip chain are:
    n=3: 1, -1/3             (2 distinct values, denominator 3)
    n=4: 1, 1/3, 0, -1/3    (4 distinct, denominator 3)
    n=5: 1, 3/5, 2/5, 1/5, 0, -1/5, -3/5  (7 distinct, denominator 5)

  The NUMBER of distinct eigenvalues: 2, 4, 7 for n=3,4,5.
  The eigenvalues are k/q for some fixed q (= 3 for n=3,4; = 5 for n=5).
  The denominator q seems to be the largest odd prime <= n.

  The MULTIPLICITY structure of eigenvalues encodes the representation
  theory of S_n on tournament space. Each eigenvalue corresponds to
  an irreducible representation, and its multiplicity is the dimension.
""")

print("opus-2026-03-16-S90dh: FAST FORMULA FOR A000568")
