#!/usr/bin/env python3
"""prime_telescope_21.py — The product of Q(1/p) over primes 2-19 = 21

DISCOVERY: PROD_{p prime, p<=19} (p+1)/(p-1) = 21 = 3*7 = THE FORBIDDEN NUMBER.

This investigates WHY this happens and what it means.
"""
from fractions import Fraction
from math import sqrt, log

def is_prime(n):
    if n < 2: return False
    for d in range(2, int(sqrt(n))+1):
        if n % d == 0: return False
    return True

primes = [p for p in range(2, 200) if is_prime(p)]

print("PRIME CAYLEY TELESCOPE: PROD Q(1/p_k) = PROD (p_k+1)/(p_k-1)")
print("="*60)
print()

# Compute the running product
prod = Fraction(1, 1)
for p in primes[:20]:
    prod *= Fraction(p+1, p-1)
    is_int = prod.denominator == 1
    marker = ""
    if prod == 21:
        marker = " <-- 21 = 3*7 = FORBIDDEN NUMBER"
    elif is_int:
        marker = f" <-- INTEGER"
    print(f"  PROD up to p={p:3d}: {str(prod):20s} = {float(prod):12.6f}{marker}")

print()

# WHERE does 21 come from? Factor the numerator and denominator
print("FACTORIZATION ANALYSIS")
print("-"*40)
print()
print("  PROD(p+1)/(p-1) for p=2,3,5,7,11,13,17,19")
print("  Numerator:   PROD(p+1) = 3*4*6*8*12*14*18*20")
num = 3*4*6*8*12*14*18*20
print(f"             = {num}")
print("  Denominator: PROD(p-1) = 1*2*4*6*10*12*16*18")
den = 1*2*4*6*10*12*16*18
print(f"             = {den}")
print(f"  Ratio: {num}/{den} = {num//den} (with remainder {num%den})")
print()

# Why does it telescope?
print("TELESCOPING STRUCTURE")
print("-"*40)
print()
print("  Write p+1 and p-1 for each prime:")
for p in primes[:8]:
    print(f"    p={p:2d}: (p+1)/(p-1) = {p+1}/{p-1}")

print()
print("  Key observation: many factors appear in BOTH numerator and denominator.")
print("  The massive cancellation leaves only the 'boundary' factors.")
print()

# Let's see which factors survive
from collections import Counter

def prime_factorize(n):
    """Return prime factorization as a Counter."""
    factors = Counter()
    d = 2
    while d*d <= n:
        while n % d == 0:
            factors[d] += 1
            n //= d
        d += 1
    if n > 1:
        factors[n] += 1
    return factors

num_factors = Counter()
den_factors = Counter()
for p in primes[:8]:
    num_factors += prime_factorize(p+1)
    den_factors += prime_factorize(p-1)

print("  Prime factorization of numerator PROD(p+1):")
for p, e in sorted(num_factors.items()):
    print(f"    {p}^{e}", end="")
print()
print("  Prime factorization of denominator PROD(p-1):")
for p, e in sorted(den_factors.items()):
    print(f"    {p}^{e}", end="")
print()
print()

# Cancel common factors
result_num = Counter(num_factors)
result_den = Counter(den_factors)
for p in set(list(num_factors.keys()) + list(den_factors.keys())):
    common = min(result_num.get(p, 0), result_den.get(p, 0))
    result_num[p] -= common
    result_den[p] -= common

print("  After cancellation:")
print("  Numerator:", end="")
prod_n = 1
for p, e in sorted(result_num.items()):
    if e > 0:
        print(f" {p}^{e}", end="")
        prod_n *= p**e
print(f" = {prod_n}")
print("  Denominator:", end="")
prod_d = 1
for p, e in sorted(result_den.items()):
    if e > 0:
        print(f" {p}^{e}", end="")
        prod_d *= p**e
print(f" = {prod_d}")
print(f"  Result: {prod_n}/{prod_d} = {prod_n//prod_d}")
print()

# The general product PROD_{p<=N} (p+1)/(p-1)
print("GENERAL PRIME CAYLEY PRODUCT")
print("-"*40)
print()
prod = Fraction(1, 1)
integer_hits = []
for i, p in enumerate(primes[:30]):
    prod *= Fraction(p+1, p-1)
    if prod.denominator == 1:
        integer_hits.append((p, int(prod)))
        print(f"  PROD up to p={p:3d}: {int(prod):10d}  *** INTEGER ***")
    else:
        print(f"  PROD up to p={p:3d}: {float(prod):14.6f}  (denom = {prod.denominator})")

print()
print("  Integer values of the prime Cayley product:")
for p, val in integer_hits:
    print(f"    Up to p={p}: product = {val}")
    if val == 21:
        print(f"      = 3 * 7 = the two FORBIDDEN-adjacent primes")
    # Factor val
    v = val
    factors = []
    d = 2
    while d*d <= v:
        while v % d == 0:
            factors.append(d)
            v //= d
        d += 1
    if v > 1:
        factors.append(v)
    if len(factors) > 1:
        print(f"      = {'*'.join(str(f) for f in factors)}")

print()

# Connection to Mertens' theorem
print("CONNECTION TO MERTENS' THEOREM")
print("-"*40)
print()
print("  Mertens' theorem: PROD_{p<=x} (1-1/p) ~ e^{-gamma}/ln(x)")
print("  where gamma = Euler-Mascheroni constant = 0.5772...")
print()
print("  Our product: PROD_{p<=x} (p+1)/(p-1) = PROD (1+1/p)/(1-1/p)")
print("  = PROD(1+1/p) / PROD(1-1/p)")
print("  = [PROD(1+1/p)] * [ln(x) * e^gamma + O(1)]")
print()
print("  By Mertens: PROD(1+1/p) ~ C * ln(x) as well (with different C).")
print("  So our product ~ C' * ln(x)^2 asymptotically.")
print()

# Compute C' = product / ln(p_max)^2 for large p_max
prod = Fraction(1, 1)
for p in primes[:30]:
    prod *= Fraction(p+1, p-1)
    ratio = float(prod) / log(p)**2
print(f"  At p=113: product/{log(113):.4f}^2 = {float(prod)/log(113)**2:.6f}")
print(f"  The constant C' ~ {float(prod)/log(113)**2:.6f}")
print()

# THE DEEPEST OBSERVATION
print("="*60)
print("THE DEEPEST OBSERVATION")
print("="*60)
print()
print("  PROD_{p<=19} Q(1/p) = 21 = 3 * 7.")
print("  3 = the curvature quantum (smallest odd cycle).")
print("  7 = the hyperbolic threshold (7 triangles per vertex).")
print("  21 = the SECOND forbidden H value.")
print()
print("  The primes involved: {2, 3, 5, 7, 11, 13, 17, 19}.")
print("  Remove {2}: {3, 5, 7, 11, 13, 17, 19} = the odd primes up to 19.")
print("  The Paley primes of our theory: {3, 7, 11} appear at positions 2, 4, 5.")
print()
print("  H(T_3)  = 3    = the first Paley prime")
print("  H(T_7)  = 189  = 3^3 * 7")
print("  H(T_11) = 95095 = 5 * 7 * 11 * 13 * 19")
print()
print("  The 8 primes in H(T_11) are {5, 7, 11, 13, 19}.")
print("  The 8 primes in the telescope are {2, 3, 5, 7, 11, 13, 17, 19}.")
print("  INTERSECTION: {5, 7, 11, 13, 19}.")
print("  DIFFERENCE: telescope has {2, 3, 17}, H(T_11) doesn't.")
print("  H(T_11) has no {2, 3, 17}. The telescope has no {empty}.")
print()
print("  17 is the Fermat prime MISSING from H(T_11).")
print("  The telescope INCLUDES 17, the Paley factorization EXCLUDES it.")
print("  The forbidden number 21 INCLUDES the factor 17 contributed.")
print()

# What about the Wallis-like product?
print("WALLIS-LIKE PRODUCT INTERPRETATION")
print("-"*40)
print()
print("  The Wallis product for pi/2:")
print("  pi/2 = PROD (2k/(2k-1)) * (2k/(2k+1))")
print("       = (2/1)(2/3) * (4/3)(4/5) * (6/5)(6/7) * ...")
print()
print("  Our product: PROD (p+1)/(p-1) over primes.")
print("  Compare to Wallis which runs over ALL integers.")
print("  Our product runs only over PRIMES — it's a 'prime Wallis product'.")
print()
print("  Wallis gives pi/2 (transcendental).")
print("  Prime Wallis gives 21 (at cutoff 19), then grows without bound.")
print("  The fact that the prime product hits an INTEGER at exactly p=19")
print("  and that integer is 21 is a number-theoretic coincidence of depth.")
