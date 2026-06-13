#!/usr/bin/env python3
"""
prime_polynomial_s257.py -- opus-2026-03-23-S257

POLYNOMIALS INVOLVING PRIMES AND THE RECURSIVE STRUCTURE

The user asks: what polynomial involving primes or their reciprocals
connects to the recursive tiling structure and the missing E(G_n) formula?

KEY CONNECTIONS:
1. The Burnside formula for T_n involves cycle types — only ODD cycle lengths contribute
2. The OCF is a FINITE EULER PRODUCT over odd cycles (primes of tournament theory)
3. The arc neutrality numerators 1, 3, 11, 79 have no known pattern — could primes help?
4. The twin SL formula involves C(fix, 2) × 2^{orbits - cycles} — prime factorization of fix matters

APPROACH: Consider the Euler product / Dirichlet series
  Z_n(s) = prod_{p odd prime, p ≤ n} (1 - p^{-s})^{-1}
restricted to primes ≤ n. This is the "tournament zeta function."

Also: the RECIPROCAL PRIME POLYNOMIAL
  R_n(x) = prod_{p prime, p ≤ n} (1 - x/p)
This encodes which primes "participate" at level n.

And: the CYCLE INDEX polynomial of S_n acting on pairs:
  Z(S_n, pairs)(x1, x2, ...) = (1/n!) sum_sigma prod x_i^{c_i(sigma)}
where c_i counts i-cycles of sigma acting on the C(n,2) pairs.

T_n = Z(S_n, pairs)(2, 2, 2, ...) — evaluate at all x_i = 2.
The PRIMES enter because only ODD cycle lengths of sigma contribute:
tournaments fixed by sigma require all EVEN cycles of sigma to have
consistent orientation, which forces the corresponding arc to be zero.

Author: opus-2026-03-23-S257
"""

import sys
from math import comb, factorial, gcd
from fractions import Fraction
from itertools import permutations
from collections import defaultdict, Counter
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def cycle_type(perm):
    n = len(perm); visited = [False]*n; cycles = []
    for i in range(n):
        if visited[i]: continue
        L = 0; j = i
        while not visited[j]:
            visited[j] = True; j = perm[j]; L += 1
        cycles.append(L)
    return tuple(sorted(cycles))

def pair_orbits(perm, n):
    """Count orbits of perm acting on unordered pairs."""
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    seen = set(); count = 0
    for p in pairs:
        if p in seen: continue
        orbit = set()
        a, b = p
        while (min(a,b), max(a,b)) not in orbit:
            orbit.add((min(a,b), max(a,b)))
            a, b = perm[a], perm[b]
        seen.update(orbit)
        count += 1
    return count

def pair_orbit_details(perm, n):
    """Return list of orbit lengths of perm on pairs, and which are 'free' vs 'forced'."""
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    seen = set(); orbits = []
    for p in pairs:
        if p in seen: continue
        orbit = []
        a, b = p
        while (min(a,b), max(a,b)) not in set(orbit):
            orbit.append((min(a,b), max(a,b)))
            a, b = perm[a], perm[b]
        seen.update(orbit)

        # Check if orbit reverses: does perm eventually swap a and b?
        a0, b0 = p
        reversed_at_some_point = False
        a, b = a0, b0
        for _ in range(len(orbit)):
            a, b = perm[a], perm[b]
            if a > b: reversed_at_some_point = True; break

        orbits.append({'length': len(orbit), 'reversed': reversed_at_some_point})
    return orbits

banner("THE BURNSIDE FORMULA AND PRIMES")

# T_n = (1/n!) * sum_sigma 2^{free_pair_orbits(sigma)}
# where free_pair_orbits = orbits of sigma on pairs that are "free"
# (not forced by the tournament constraint).

# A pair orbit of length L is:
# - FREE (2 choices) if sigma never reverses the pair in its orbit
# - FORCED (1 choice) if sigma reverses (L is even in the pair orbit)

# For sigma to fix ANY tournament: all pair orbits with reversal must
# have even length. If a pair orbit has odd length with reversal,
# it forces A[a][b] = A[b][a] = 1, impossible for tournaments.

# So: Fix(sigma) = 0 if sigma has any pair orbit that is odd AND reversed.
# Otherwise: Fix(sigma) = 2^{free_orbits} × 1^{forced_orbits}

# The CYCLE TYPE of sigma determines this:
# - Sigma with any even cycle has a pair orbit that reverses (the pair
#   of vertices in the same even cycle). This forces Fix = 0.
# - Sigma with ALL ODD cycles: all pair orbits are either free or
#   consistently forced. Fix = 2^{free_orbits}.

# THIS IS WHY ONLY ODD-CYCLE-TYPE PERMUTATIONS CONTRIBUTE TO T_n.

print("  The Burnside formula for T_n:\n")
print("  T_n = (1/n!) * sum_{sigma with all odd cycles} 2^{pair_orbits(sigma)}\n")

for n in [3, 4, 5, 6]:
    total = 0
    by_type = defaultdict(lambda: [0, 0])  # [count, sum_fix]

    for perm in permutations(range(n)):
        ct = cycle_type(list(perm))
        # Check: all cycles odd?
        all_odd = all(c % 2 == 1 for c in ct)
        if not all_odd: continue

        po = pair_orbits(list(perm), n)
        fix = 2**po

        by_type[ct][0] += 1
        by_type[ct][1] += fix
        total += fix

    T_n = total // factorial(n)
    print("  n=%d: T_n = %d" % (n, T_n))
    for ct in sorted(by_type.keys()):
        cnt, s = by_type[ct]
        # What primes appear in cycle lengths?
        primes_in_ct = set()
        for c in ct:
            for p in range(2, c+1):
                if c % p == 0 and all(c % q != 0 for q in range(2, p)):
                    primes_in_ct.add(p)
        print("    %s: %d perms, Fix_sum=%d, primes=%s" %
              (ct, cnt, s, sorted(primes_in_ct) if primes_in_ct else '{1}'))
    print()

banner("THE PRIME POLYNOMIAL FOR T_n")

# From the Burnside formula:
# T_n = (1/n!) * sum over partitions lambda of n into odd parts:
#   (n! / prod c_k^{a_k} a_k!) × 2^{pair_orbits(lambda)}
# where the sum is over all such lambda.

# The pair_orbits of a permutation with cycle type (1^a1, 3^a3, 5^a5, ...)
# can be computed from the cycle type alone.

# For a permutation with cycle lengths c1, c2, ...:
# pair_orbits = sum_{i<j} gcd(ci, cj) + sum_i floor(ci/2)
# (orbits from pairs within the same cycle + pairs across cycles)

# Wait: more carefully:
# Pairs within a k-cycle: C(k,2) pairs, forming orbits of size k
#   → floor(C(k,2)/k) = floor(k(k-1)/(2k)) = floor((k-1)/2) orbits for k odd
#   But some orbits reverse. For k odd: none reverse.
#   So: (k-1)/2 free orbits per k-cycle (k odd).
# Pairs across a k-cycle and an l-cycle: k*l pairs, orbits of size lcm(k,l)
#   → k*l / lcm(k,l) = gcd(k,l) orbits, all free (no reversal for odd k,l)

# TOTAL pair_orbits for cycle type (c1, c2, ...):
# = sum_i (c_i - 1)/2 + sum_{i<j} gcd(c_i, c_j)

def pair_orbits_from_type(ct):
    """Compute number of free pair orbits from cycle type (all odd)."""
    within = sum((c - 1) // 2 for c in ct)
    across = sum(gcd(ct[i], ct[j]) for i in range(len(ct)) for j in range(i+1, len(ct)))
    return within + across

# Verify
print("  Pair orbits from cycle type formula:\n")
for n in [3, 4, 5, 6]:
    total_fix = 0
    for perm in permutations(range(n)):
        ct = cycle_type(list(perm))
        if not all(c % 2 == 1 for c in ct): continue
        po_direct = pair_orbits(list(perm), n)
        po_formula = pair_orbits_from_type(ct)
        if po_direct != po_formula:
            print("  MISMATCH at n=%d, ct=%s: direct=%d, formula=%d" %
                  (n, ct, po_direct, po_formula))
        total_fix += 2**po_direct

    T_check = total_fix // factorial(n)
    T_known = {3:4, 4:16, 5:88, 6:704}[n]
    print("  n=%d: T from formula = %d (expected %d, match=%s)" %
          (n, T_check, T_known, T_check == T_known))

banner("THE PRIME CONTENT OF CYCLE TYPES")

# For each n, list the cycle types that contribute to T_n,
# along with the PRIME FACTORIZATION of each cycle length.

print("  Cycle types contributing to T_n and their prime content:\n")
for n in range(3, 10):
    # Partitions of n into odd parts
    partitions = []
    def gen_odd_partitions(remaining, max_part, current):
        if remaining == 0:
            partitions.append(tuple(sorted(current)))
            return
        for p in range(min(remaining, max_part), 0, -2):  # odd parts only
            gen_odd_partitions(remaining - p, p, current + [p])
    gen_odd_partitions(n, n if n % 2 == 1 else n-1, [])

    po_total = 0
    for partition in sorted(set(partitions)):
        po = pair_orbits_from_type(partition)
        # Number of permutations with this cycle type
        ct_counter = Counter(partition)
        n_perms = factorial(n) // reduce(lambda a,b: a*b,
                    [c**ct_counter[c] * factorial(ct_counter[c]) for c in ct_counter], 1)
        contribution = n_perms * 2**po
        po_total += contribution

        # Prime factorization of cycle lengths
        cycle_primes = set()
        for c in partition:
            if c == 1: continue
            for p in range(3, c+1, 2):
                if c % p == 0 and all(c % q != 0 for q in range(2, p)):
                    cycle_primes.add(p)

        if n <= 7:
            print("  n=%d: %s → po=%d, #perms=%d, Fix=%d, primes=%s" %
                  (n, partition, po, n_perms, 2**po, sorted(cycle_primes) if cycle_primes else '{1}'))

    T = po_total // factorial(n)
    T_known = {3:4, 4:16, 5:88, 6:704, 7:8912, 8:188288, 9:6847200}
    print("  → T_%d = %d (expected %s)\n" % (n, T, T_known.get(n, '?')))

banner("THE PRIME POLYNOMIAL P(n,x)")

# Define: P(n, x) = sum over odd partitions lambda of n:
#   (n! / prod c^a a!) × x^{pair_orbits(lambda)}
#
# Then T_n = P(n, 2) / n!

# The polynomial P(n, x) encodes ALL information about the Burnside sum.
# Its coefficients at each power of x count how many (weighted) cycle types
# contribute that many free orbits.

print("  The polynomial P(n, x) such that T_n = P(n, 2) / n!:\n")

for n in range(3, 8):
    coeffs = defaultdict(int)  # power of x -> coefficient
    partitions = []
    def gen_op(remaining, max_part, current):
        if remaining == 0:
            partitions.append(tuple(sorted(current)))
            return
        for p in range(min(remaining, max_part), 0, -2):
            gen_op(remaining - p, p, current + [p])
    partitions = []
    gen_op(n, n if n % 2 == 1 else n-1, [])

    for partition in sorted(set(partitions)):
        po = pair_orbits_from_type(partition)
        ct_counter = Counter(partition)
        n_perms = factorial(n) // reduce(lambda a,b: a*b,
                    [c**ct_counter[c] * factorial(ct_counter[c]) for c in ct_counter], 1)
        coeffs[po] += n_perms

    # Print polynomial
    terms = []
    for power in sorted(coeffs.keys()):
        terms.append("%d x^%d" % (coeffs[power], power))
    poly_str = " + ".join(terms)
    T_check = sum(c * 2**p for p, c in coeffs.items()) // factorial(n)

    print("  P(%d, x) = %s" % (n, poly_str))
    print("  P(%d, 2) / %d! = %d = T_%d  ✓\n" % (n, n, T_check, n))

banner("THE RECIPROCAL PRIME CONNECTION")

print("""
  NOW: How do primes connect to the RECURSIVE STRUCTURE?

  The recursion T(n) decomposes along the staircase:
  T(n) = bottom(n-1) + top(n-1) - overlap(n-2) + apex_correction

  In the Burnside polynomial P(n, x):
  - Each term x^k comes from a cycle type with k free pair orbits
  - The cycle type's odd primes determine WHICH arcs are free

  THE KEY QUESTION: Does P(n, x) have a RECURSIVE formula?

  P(n, x) depends on the partitions of n into odd parts.
  The generating function for odd partitions:
  prod_{k=0}^{inf} 1/(1 - z^{2k+1}) = 1 + z + z^2 + 2z^3 + 2z^4 + 3z^5 + ...

  And the pair_orbits function maps partition → integer.
  So P(n, x) = sum_lambda (multinomial coeff) × x^{pair_orbits(lambda)}.

  IF the pair_orbits function has a recursive structure:
  pair_orbits(lambda + [1]) = pair_orbits(lambda) + sum gcd(1, c_i) = pair_orbits(lambda) + len(lambda)
  (adding a 1-cycle adds gcd(1, c) = 1 for each existing cycle)

  pair_orbits(lambda + [3]) = pair_orbits(lambda) + 1 + sum gcd(3, c_i)
  (the within-orbit for a 3-cycle is (3-1)/2 = 1, plus cross-orbits)

  So: ADDING A PRIME CYCLE OF LENGTH p TO A PARTITION:
  Δ(pair_orbits) = (p-1)/2 + sum_{c in lambda} gcd(p, c)

  The gcd(p, c) = p if p | c, else gcd(p, c) depends on common factors.
  For p prime: gcd(p, c) = p if p | c, else 1.

  So: Δ = (p-1)/2 + #{cycles divisible by p} × (p-1) + #{other cycles}
      = (p-1)/2 + (p-1) × d_p + (len(lambda) - d_p)
      = (p-1)/2 + (p-1) × d_p + len(lambda) - d_p
      = (p-1)/2 + (p-2) × d_p + len(lambda)

  where d_p = #{existing cycles with length divisible by p}.

  For p = 1 (fixed point): Δ = 0 + len(lambda) = len(lambda).
  For p = 3: Δ = 1 + 2 × d_3 + len(lambda).
  For p = 5: Δ = 2 + 3 × d_5 + len(lambda).
  For p = 7: Δ = 3 + 5 × d_7 + len(lambda).

  THE PATTERN: Adding a p-cycle increases pair_orbits by
  (p-1)/2 + (p-2) × d_p + k  (where k = # existing cycles).

  This gives a RECURSIVE formula for P(n, x) in terms of P(n-p, x)
  for each odd prime p ≤ n.
""")

# Verify the recursion
print("  Verifying pair_orbits recursion:\n")
# partition (1,1,1) for n=3: po = 0 + 3*1 = 3? Let me check.
# Actually: pair_orbits((1,1,1)) = within: 3×0 = 0, across: 3×gcd(1,1) = 3. Total = 3.
# Adding a 1-cycle to (1,1): Δ = len((1,1)) = 2. New po = 1 + 2 = 3. Old po((1,1)) = 1 + gcd(1,1) = 1. Wait.
# pair_orbits((1,1)) = within: 0+0 = 0, across: gcd(1,1) = 1. Total = 1.
# pair_orbits((1,1,1)) = within: 0, across: 3. Total = 3.
# Δ = 3 - 1 = 2 = len((1,1)) = 2. ✓

# partition (3) for n=3: po = (3-1)/2 = 1. Total = 1.
# Correct: a single 3-cycle on 3 vertices has 1 pair orbit.

# partition (1,3) for n=4: po = 0 + (3-1)/2 + gcd(1,3) = 0 + 1 + 1 = 2.
# Check: pair_orbits((1,3)) should be 2.
print("  po((1,1,1)) = %d (expect 3)" % pair_orbits_from_type((1,1,1)))
print("  po((3,)) = %d (expect 1)" % pair_orbits_from_type((3,)))
print("  po((1,3)) = %d (expect 2)" % pair_orbits_from_type((1,3)))
print("  po((1,1,3)) = %d (expect ?)" % pair_orbits_from_type((1,1,3)))
print("  po((3,3)) = %d (expect ?)" % pair_orbits_from_type((3,3)))
print("  po((1,1,1,1,1)) = %d (expect C(5,2)=10)" % pair_orbits_from_type((1,1,1,1,1)))
print("  po((5,)) = %d" % pair_orbits_from_type((5,)))
print("  po((1,1,3)) = %d" % pair_orbits_from_type((1,1,3)))
print("  po((1,1,1,3)) = %d" % pair_orbits_from_type((1,1,1,3)))
print("  po((3,3)) = %d" % pair_orbits_from_type((3,3)))

banner("THE GENERATING FUNCTION WITH PRIMES")

print("""
  The generating function for T_n involves:

  Z(x, z) = sum_n (P(n, x) / n!) × z^n

  = sum over all sequences (a_1, a_3, a_5, ...) of non-neg integers:
    prod_{k odd} (z^k / k)^{a_k} / a_k! × x^{pair_orbits}

  This is an EXPONENTIAL generating function where each odd k contributes
  z^k/k per cycle, and the pair_orbits add multiplicatively.

  The pair_orbits from a single k-cycle: (k-1)/2.
  The pair_orbits from two cycles of lengths k and l: gcd(k, l).

  So the GF factorizes PARTLY:
  Within-cycle: prod_{k odd} sum_{a≥0} (z^k/(k × a!))^a × x^{a(k-1)/2}
  Cross-cycle: harder (gcd terms couple different cycle types).

  FOR PRIMES: if k is prime, gcd(k, l) = k if l ≡ 0 (mod k), else 1.
  So prime cycles interact with other cycles in a BINARY way:
  either full coupling (gcd = k) or minimal coupling (gcd = 1).

  This makes the GF SIMPLER for prime cycle types than composite ones!

  THE PRIME APPROXIMATION:
  If we only keep PRIME cycle types (1, 3, 5, 7, ...):
  Z_prime(x, z) ≈ exp(z × x^0) × prod_{p prime, p≥3} exp(z^p/p × x^{(p-1)/2})
  × (cross-terms depending on gcd)

  The cross-term for two prime cycles of different primes p, q:
  gcd(p, q) = 1 (since p ≠ q are both prime).
  So cross contribution = +1 per pair of cycles from different primes.

  The cross-term for two cycles of the SAME prime p:
  gcd(p, p) = p. So cross contribution = +p per pair.

  THIS IS WHERE THE PRIMES ENTER THE RECURSION.
""")

# Compute T_n using the cycle-type formula to verify at larger n
print("  T_n from cycle-type formula:\n")
T_known = {3:4, 4:16, 5:88, 6:704, 7:8912, 8:188288, 9:6847200}
for n in range(3, 10):
    partitions_set = set()
    def gen_op2(remaining, max_part, current):
        if remaining == 0:
            partitions_set.add(tuple(sorted(current)))
            return
        for p in range(min(remaining, max_part), 0, -2):
            gen_op2(remaining - p, p, current + [p])
    gen_op2(n, n if n % 2 == 1 else n-1, [])

    total = 0
    for partition in partitions_set:
        po = pair_orbits_from_type(partition)
        ct_counter = Counter(partition)
        n_perms = factorial(n) // reduce(lambda a,b: a*b,
                    [c**ct_counter[c] * factorial(ct_counter[c]) for c in ct_counter], 1)
        total += n_perms * 2**po

    T = total // factorial(n)
    expected = T_known.get(n, '?')
    print("  T_%d = %d (expected %s, match=%s)" % (n, T, expected, T == expected))

print("\nDone. opus-2026-03-23-S257")
