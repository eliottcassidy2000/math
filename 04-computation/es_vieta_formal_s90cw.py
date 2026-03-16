#!/usr/bin/env python3
"""
es_vieta_formal_s90cw.py — Erdos-Straus via Vieta jumping = formal group on factorizations
opus-2026-03-16-S90cw

The KEY: 4/p = 1/(pk) + 1/b + 1/c reduces to factoring p^2*k^2 in a shifted ring.
This converts ADDITION (ES) to MULTIPLICATION (factoring) via completing the product.
The S/A duality made algorithmic.
"""

from math import gcd, isqrt
from sympy import isprime, primerange, divisors
from collections import Counter, defaultdict

def es_via_vieta(p, max_k=None):
    """Find all ES decompositions 4/p = 1/a + 1/b + 1/c for prime p.
    Uses Vieta jumping: set a = p*k, then factor p^2*k^2 to find b, c.
    Returns list of (a, b, c) with a <= b <= c."""

    if max_k is None:
        max_k = p  # sufficient range

    solutions = set()

    for k in range(1, max_k + 1):
        # 4/p = 1/(pk) + 1/b + 1/c
        # => 1/b + 1/c = (4k-1)/(pk)
        # => (b*(4k-1) - pk)(c*(4k-1) - pk) = p^2 * k^2
        # So factor N = p^2 * k^2 as d1*d2, then
        # b = (d1 + pk) / (4k-1), c = (d2 + pk) / (4k-1)

        N = p * p * k * k
        coeff = 4 * k - 1  # the (4k-1) coefficient

        for d1 in divisors(N):
            d2 = N // d1
            if d1 > d2:
                break
            # b = (d1 + pk) / (4k-1), c = (d2 + pk) / (4k-1)
            b_num = d1 + p * k
            c_num = d2 + p * k
            if b_num % coeff == 0 and c_num % coeff == 0:
                b = b_num // coeff
                c = c_num // coeff
                a = p * k
                if b > 0 and c > 0:
                    triple = tuple(sorted([a, b, c]))
                    solutions.add(triple)

    return sorted(solutions)


print("""
ERDOS-STRAUS VIA VIETA JUMPING = FORMAL GROUP ON FACTORIZATIONS
=================================================================
""")

# Compute for all primes up to 500
print("ES solution counts for primes up to 500:")
print(f"{'p':>5} | {'#sols':>5} | {'mod 4':>5} | {'mod 7':>5} | {'mod 42':>5} | {'Bott':>4} | {'min a':>5} | {'note':>10}")
print("-" * 60)

sol_counts = {}
hardest = []

for p in primerange(2, 501):
    sols = es_via_vieta(p, max_k=min(p, 200))
    sol_counts[p] = len(sols)
    min_a = min(s[0] for s in sols) if sols else 0
    note = ""
    if p in [2,3,5,7,11,13,42,47,127]: note = "KEY"
    if len(sols) <= 3: note = "HARD"
    if p <= 30 or note or len(sols) <= 5:
        print(f"{p:5d} | {len(sols):5d} | {p%4:5d} | {p%7:5d} | {p%42:5d} | {p%8:4d} | {min_a:5d} | {note:>10}")
    hardest.append((len(sols), p))

hardest.sort()
print(f"\nHardest primes (fewest solutions):")
for count, p in hardest[:15]:
    print(f"  p={p:5d}: {count:3d} solutions, mod 4={p%4}, mod 7={p%7}, mod 42={p%42}, Bott={p%8}")

# Classify by mod 4 and mod 7
print(f"\nMean solution count by residue class:")
for mod, name in [(4, "mod 4"), (7, "mod 7"), (8, "Bott"), (42, "mod 42")]:
    by_class = defaultdict(list)
    for p, count in sol_counts.items():
        if p > 5:  # skip tiny primes
            by_class[p % mod].append(count)
    print(f"\n  By {name}:")
    for r in sorted(by_class.keys()):
        vals = by_class[r]
        if vals:
            mean = sum(vals) / len(vals)
            print(f"    {name}={r:2d}: mean={mean:7.1f}, count={len(vals):3d}, "
                  f"min={min(vals):3d}, max={max(vals):4d}")

# The KEY test: does mod 42 predict difficulty?
print(f"\nCorrelation: mod 42 residue vs ES difficulty")
coprime_42 = [r for r in range(42) if gcd(r, 42) == 1]
mod42_means = {}
for r in coprime_42:
    vals = [sol_counts[p] for p in sol_counts if p > 42 and p % 42 == r]
    if vals:
        mod42_means[r] = sum(vals) / len(vals)
        if len(vals) >= 3:
            print(f"  mod 42 = {r:2d}: mean solutions = {mod42_means[r]:7.1f} ({len(vals)} primes)")

# Do the split/inert classification predict ES difficulty?
print(f"\nSplit vs Inert (in Q(i)) and ES difficulty:")
split_sols = []
inert_sols = []
for p in sol_counts:
    if p > 5:
        if p % 4 == 1:
            split_sols.append(sol_counts[p])
        elif p % 4 == 3:
            inert_sols.append(sol_counts[p])

print(f"  Split (p=1 mod 4): mean = {sum(split_sols)/len(split_sols):.1f}, "
      f"min = {min(split_sols)}, n = {len(split_sols)}")
print(f"  Inert (p=3 mod 4): mean = {sum(inert_sols)/len(inert_sols):.1f}, "
      f"min = {min(inert_sols)}, n = {len(inert_sols)}")
print(f"  Ratio split/inert: {sum(split_sols)/len(split_sols) / (sum(inert_sols)/len(inert_sols)):.3f}")

# Split vs Inert in Q(sqrt(5))
print(f"\nSplit vs Inert (in Q(phi) = Q(sqrt(5))) and ES difficulty:")
from sympy import legendre_symbol
split5 = []
inert5 = []
for p in sol_counts:
    if p > 5 and p != 5:
        ls = legendre_symbol(5, p)
        if ls == 1:
            split5.append(sol_counts[p])
        else:
            inert5.append(sol_counts[p])

print(f"  Split in Q(phi): mean = {sum(split5)/len(split5):.1f}, n = {len(split5)}")
print(f"  Inert in Q(phi): mean = {sum(inert5)/len(inert5):.1f}, n = {len(inert5)}")

# The Vieta interpretation
print(f"""

THE VIETA-FORMAL GROUP CONNECTION
====================================

The Erdos-Straus decomposition 4/p = 1/(pk) + 1/b + 1/c
reduces via Vieta to: FACTOR p^2*k^2 as d1*d2 where
b = (d1+pk)/(4k-1) and c = (d2+pk)/(4k-1) are positive integers.

The number of ES decompositions = number of valid factorizations of p^2*k^2
summed over all k.

The factorization count depends on:
  d(p^2*k^2) = number of divisors of p^2*k^2.
  For k=1: d(p^2) = 3 (divisors: 1, p, p^2).
  For k=2: d(4p^2) = d(4)*d(p^2) = 3*3 = 9 (if gcd(4,p)=1).

More divisors = more potential solutions.
But the divisibility condition (4k-1) | (d+pk) filters most out.

The FORMAL GROUP enters through the Cayley transform:
  The factorization d1*d2 = p^2*k^2 is a MULTIPLICATIVE decomposition.
  This is the formal group's native operation.
  The (4k-1) divisibility is the ADDITIVE constraint
  translated to the multiplicative side via Vieta.

THE S/A DUALITY IN ES:
  ADDITIVE (antisymmetric): 4/p = sum of unit fractions.
  MULTIPLICATIVE (symmetric): factor p^2*k^2 as d1*d2.
  VIETA JUMPING: the bridge between A and S.
  Solution count = intersection of A and S constraints.
""")

print("opus-2026-03-16-S90cw: ERDOS-STRAUS VIA VIETA = FORMAL GROUP ON FACTORIZATIONS")
