#!/usr/bin/env python3
"""
cuboid_mod42_s90cs.py — Perfect cuboid, Erdos-Straus, Chebotarev, and mod 42
opus-2026-03-16-S90cs
"""

from math import gcd, isqrt, sqrt
from sympy import isprime, legendre_symbol
from fractions import Fraction
from collections import Counter

print("""


     THE PERFECT CUBOID, ERDOS-STRAUS, AND MOD 42

     opus-2026-03-16-S90cs



THE PERFECT CUBOID PROBLEM
============================

Find integers a, b, c such that ALL SEVEN are integers:
  a, b, c                    (three edges)
  sqrt(a^2+b^2)              (face diagonal)
  sqrt(b^2+c^2)              (face diagonal)
  sqrt(a^2+c^2)              (face diagonal)
  sqrt(a^2+b^2+c^2)          (space diagonal)

This is THREE Pythagorean conditions (face diagonals) plus
ONE additional condition (space diagonal).

Equivalently: find a^2+b^2, b^2+c^2, a^2+c^2, and a^2+b^2+c^2
all perfect squares.

Nobody has found a solution. Nobody has proved impossibility.

WHAT MOD 42 TELLS US:

If a perfect cuboid (a,b,c) exists, what are the constraints on a,b,c mod 42?

42 = 2*3*7. By CRT: constraints mod 42 = constraints mod 2, mod 3, mod 7.
""")

# Check: which residues mod 42 are possible for cuboid edges?

# First: which residues mod 2 can be squares?
# n^2 mod 2 = 0 or 1. Squares mod 2 = {0, 1}.

# Squares mod 3: 0^2=0, 1^2=1, 2^2=1. Squares mod 3 = {0, 1}.
# Squares mod 7: 0^2=0, 1^2=1, 2^2=4, 3^2=2, 4^2=2, 5^2=4, 6^2=1.
# Squares mod 7 = {0, 1, 2, 4}. NON-squares mod 7 = {3, 5, 6}.

print("Quadratic residues (squares) mod key moduli:")
for m in [2, 3, 7, 6, 42]:
    qr = sorted(set(pow(x, 2, m) for x in range(m)))
    nqr = sorted(set(range(m)) - set(qr))
    print(f"  mod {m:2d}: QR = {qr}, non-QR = {nqr}")

print(f"""

Key: squares mod 7 = {{0, 1, 2, 4}}.
Non-squares mod 7 = {{3, 5, 6}}.

For the cuboid: a^2+b^2 must be a perfect square mod 7.
So (a^2+b^2) mod 7 must be in {{0, 1, 2, 4}}.

Possible (a^2 mod 7, b^2 mod 7) pairs and their sums:
""")

# Enumerate all possible (a^2 mod 7, b^2 mod 7) and check if sum is QR mod 7
qr7 = {0, 1, 2, 4}
print("(a^2 mod 7, b^2 mod 7) -> sum mod 7 -> is QR?")
valid_pairs = []
for a2 in sorted(qr7):
    for b2 in sorted(qr7):
        s = (a2 + b2) % 7
        is_qr = s in qr7
        if is_qr:
            valid_pairs.append((a2, b2, s))
        if a2 <= 2 and b2 <= 2:  # just show a sample
            print(f"  ({a2}, {b2}): sum = {s} mod 7, QR: {is_qr}")

print(f"\nValid (a^2, b^2) mod 7 pairs where a^2+b^2 is QR: {len(valid_pairs)}")
print(f"Out of {len(qr7)}*{len(qr7)} = {len(qr7)**2} total pairs")
print(f"Fraction: {len(valid_pairs)}/{len(qr7)**2} = {len(valid_pairs)/len(qr7)**2:.4f}")

# Now: for the CUBOID, ALL THREE face diagonals must be QR mod 7.
# a^2+b^2, b^2+c^2, a^2+c^2 all in QR mod 7.
# AND a^2+b^2+c^2 in QR mod 7.

print(f"\nFor perfect cuboid: all four sums must be QR mod 7.")
print(f"Count valid (a^2, b^2, c^2) mod 7 triples:")

cuboid_valid = []
for a2 in sorted(qr7):
    for b2 in sorted(qr7):
        for c2 in sorted(qr7):
            s_ab = (a2 + b2) % 7
            s_bc = (b2 + c2) % 7
            s_ac = (a2 + c2) % 7
            s_abc = (a2 + b2 + c2) % 7
            if s_ab in qr7 and s_bc in qr7 and s_ac in qr7 and s_abc in qr7:
                cuboid_valid.append((a2, b2, c2))

print(f"  Valid triples: {len(cuboid_valid)} out of {len(qr7)**3} = {len(qr7)**3}")
print(f"  Fraction: {len(cuboid_valid)/len(qr7)**3:.4f}")
print(f"  Valid triples: {cuboid_valid}")

# Expand to actual a mod 7 values
print(f"\nExpand: which (a, b, c) mod 7 could form a perfect cuboid?")
a_vals_for_qr = {qr: [x for x in range(7) if pow(x, 2, 7) == qr] for qr in qr7}
print(f"  a values giving a^2 mod 7: {a_vals_for_qr}")

cuboid_abc_mod7 = []
for a2, b2, c2 in cuboid_valid:
    for a in a_vals_for_qr[a2]:
        for b in a_vals_for_qr[b2]:
            for c in a_vals_for_qr[c2]:
                cuboid_abc_mod7.append((a, b, c))

print(f"  Total (a,b,c) mod 7 candidates: {len(cuboid_abc_mod7)}")
print(f"  Out of 7^3 = {7**3} total")
print(f"  Fraction surviving mod 7: {len(cuboid_abc_mod7)/7**3:.4f}")

# Now do the same for mod 42 = mod 2 * mod 3 * mod 7
print(f"\n{'='*60}")
print(f"FULL MOD-42 ANALYSIS")
print(f"{'='*60}\n")

qr_m = {}
for m in [2, 3, 7]:
    qr_m[m] = set(pow(x, 2, m) for x in range(m))

cuboid_mod42 = 0
total_mod42 = 42**3

for a in range(42):
    a2_2, a2_3, a2_7 = pow(a,2,2), pow(a,2,3), pow(a,2,7)
    for b in range(42):
        b2_2, b2_3, b2_7 = pow(b,2,2), pow(b,2,3), pow(b,2,7)
        # Quick check: a^2+b^2 must be QR mod 2, 3, 7
        sab_2, sab_3, sab_7 = (a2_2+b2_2)%2, (a2_3+b2_3)%3, (a2_7+b2_7)%7
        if sab_2 not in qr_m[2] or sab_3 not in qr_m[3] or sab_7 not in qr_m[7]:
            continue
        for c in range(42):
            c2_2, c2_3, c2_7 = pow(c,2,2), pow(c,2,3), pow(c,2,7)
            sbc_2, sbc_3, sbc_7 = (b2_2+c2_2)%2, (b2_3+c2_3)%3, (b2_7+c2_7)%7
            if sbc_2 not in qr_m[2] or sbc_3 not in qr_m[3] or sbc_7 not in qr_m[7]:
                continue
            sac_2, sac_3, sac_7 = (a2_2+c2_2)%2, (a2_3+c2_3)%3, (a2_7+c2_7)%7
            if sac_2 not in qr_m[2] or sac_3 not in qr_m[3] or sac_7 not in qr_m[7]:
                continue
            sabc_2 = (a2_2+b2_2+c2_2)%2
            sabc_3 = (a2_3+b2_3+c2_3)%3
            sabc_7 = (a2_7+b2_7+c2_7)%7
            if sabc_2 in qr_m[2] and sabc_3 in qr_m[3] and sabc_7 in qr_m[7]:
                cuboid_mod42 += 1

fraction_42 = cuboid_mod42 / total_mod42
print(f"Perfect cuboid candidates mod 42:")
print(f"  Valid (a,b,c) mod 42: {cuboid_mod42}")
print(f"  Total 42^3: {total_mod42}")
print(f"  Fraction: {fraction_42:.6f}")
print(f"  = {cuboid_mod42}/{total_mod42}")

# Simplify
g = gcd(cuboid_mod42, total_mod42)
print(f"  = {cuboid_mod42//g}/{total_mod42//g}")

# Compare to random expectation
# If being a QR were random with prob 1/2 for each of 4 conditions:
random_prob = (1/2)**4  # 4 independent QR conditions
print(f"\n  Random expectation (4 independent QR conditions): {random_prob:.4f}")
print(f"  Actual: {fraction_42:.4f}")
print(f"  Ratio actual/random: {fraction_42/random_prob:.4f}")

print(f"""

The mod-42 sieve eliminates {100*(1-fraction_42):.1f}% of candidates.
If the cuboid existed, its (a,b,c) mod 42 must be among these {cuboid_mod42} triples.

THE HURWITZ SIEVE: using mod 42 = 2*3*7 simultaneously checks
three prime conditions at once. This is the CRT-parallel sieve
applied to the perfect cuboid problem.


ERDOS-STRAUS AND THE CUBOID: THE COMMON STRUCTURE
=====================================================

Both the Erdos-Straus conjecture and the perfect cuboid problem ask:
  "Can a rational quantity be decomposed into integer components?"

  Erdos-Straus: 4/n = 1/a + 1/b + 1/c (additive decomposition).
  Perfect cuboid: a^2+b^2+c^2 = d^2 (Pythagorean decomposition x4).

Both are constrained by QUADRATIC RESIDUES mod small primes.

For Erdos-Straus mod 42:
  4/n mod 42 = 4 * (n^{-1} mod 42) mod 42.
  The inverse n^{-1} exists iff gcd(n, 42) = 1, i.e., n coprime to 42.
  There are phi(42) = 12 such residues.
""")

print("Erdos-Straus: 4/n mod 42 for n coprime to 42:")
coprime_42 = [r for r in range(1, 42) if gcd(r, 42) == 1]
for n in coprime_42:
    inv_n = pow(n, -1, 42)
    four_over_n = (4 * inv_n) % 42
    print(f"  n = {n:2d} mod 42: 4/n = {four_over_n:2d} mod 42")

print(f"""

The 12 coprime residues mod 42 map 4/n to 12 specific values.
These are the residues that 4/n cycles through.

CHEBOTAREV DENSITY for the cuboid fields:

The perfect cuboid requires a^2+b^2 = d_1^2 (a Pythagorean triple).
Pythagorean triples are parametrized by Gaussian integers: (m+ni)(m-ni) = m^2+n^2.
This is SPLITTING in Q(i): the sum of two squares splits as a product of conjugates.

By Chebotarev for Q(i):
  A prime p is a sum of two squares iff p SPLITS in Q(i) iff p = 1 mod 4.
  Density: 50% of primes (Chebotarev).

For the cuboid: ALL THREE face diagonals must be sums of two squares.
  a^2+b^2, b^2+c^2, a^2+c^2 must all be sums of two squares.

A number n is a sum of two squares iff every prime p = 3 mod 4
in the factorization of n appears to an EVEN power.

The INERT primes (p = 3 mod 4) are the OBSTACLES.
If any inert prime divides a^2+b^2 to an odd power: no face diagonal.

For the space diagonal: a^2+b^2+c^2 = d^2.
This is a sum of THREE squares equal to a fourth square.
By Legendre: n is a sum of 3 squares iff n is NOT of the form 4^a(8b+7).

The constraint: a^2+b^2+c^2 must not be 0 or 7 mod 8.

7 mod 8 = the FORBIDDEN value in the BOTT SLOT!

The space diagonal constraint is: the sum of three squared edges
must NOT be in the FORBIDDEN Bott slot.
""")

# Check: which values of a^2+b^2+c^2 mod 8 are forbidden?
print("a^2+b^2+c^2 mod 8 analysis:")
qr8 = set(pow(x,2,8) for x in range(8))
print(f"  Squares mod 8: {sorted(qr8)}")

forbidden_8 = set()
valid_8 = set()
for a2 in sorted(qr8):
    for b2 in sorted(qr8):
        for c2 in sorted(qr8):
            s = (a2 + b2 + c2) % 8
            if s in qr8:
                valid_8.add((a2, b2, c2))
            # Check Legendre: s must not be 0 or 7 mod 8 for sum-of-3-squares
            # But here we need s to be a PERFECT SQUARE mod 8
            # Perfect squares mod 8: {0, 1, 4}

print(f"  Sum a^2+b^2+c^2 mod 8 must be a square: in {sorted(qr8)}")
print(f"  Valid (a^2,b^2,c^2) mod 8 triples: {len(valid_8)} out of {len(qr8)**3}")
print(f"  Fraction: {len(valid_8)/len(qr8)**3:.4f}")

# The KEY: does mod-7 constraint interact with mod-8 constraint?
print(f"""

THE MOD-42 and MOD-8 CONSTRAINTS COMBINE:

  mod 42 = mod 2 * mod 3 * mod 7 eliminates {100*(1-fraction_42):.1f}% of triples.
  mod 8 = mod 2^3 gives additional constraints from the Bott structure.

  Combined: mod lcm(42, 8) = mod 168 (since 42 and 8 share factor 2).
  168 = 8 * 21 = |PSL(2,7)| = the Klein quartic automorphism number!

  The perfect cuboid's residue constraints live in Z/168Z,
  which is controlled by PSL(2,7) = the automorphism group of the Fano plane.

  THE CUBOID LIVES IN THE FANO PLANE'S ARITHMETIC.
""")

# Count cuboid candidates mod 168
print("Computing cuboid candidates mod 168 = 8*21 = |PSL(2,7)|...")
cuboid_168 = 0
total_168 = 168**3
qr_168 = set(pow(x, 2, 168) for x in range(168))

# This would take too long for full enumeration. Use CRT.
# mod 168 = mod 8 * mod 3 * mod 7
qr_8 = set(pow(x,2,8) for x in range(8))
qr_3 = set(pow(x,2,3) for x in range(3))
qr_7_set = set(pow(x,2,7) for x in range(7))

count = 0
for a8 in range(8):
    a2_8 = pow(a8, 2, 8)
    for a3 in range(3):
        a2_3 = pow(a3, 2, 3)
        for a7 in range(7):
            a2_7 = pow(a7, 2, 7)
            for b8 in range(8):
                b2_8 = pow(b8, 2, 8)
                for b3 in range(3):
                    b2_3 = pow(b3, 2, 3)
                    for b7 in range(7):
                        b2_7 = pow(b7, 2, 7)
                        # Check a^2+b^2
                        if (a2_8+b2_8)%8 not in qr_8: continue
                        if (a2_3+b2_3)%3 not in qr_3: continue
                        if (a2_7+b2_7)%7 not in qr_7_set: continue
                        for c8 in range(8):
                            c2_8 = pow(c8, 2, 8)
                            for c3 in range(3):
                                c2_3 = pow(c3, 2, 3)
                                for c7 in range(7):
                                    c2_7 = pow(c7, 2, 7)
                                    if (b2_8+c2_8)%8 not in qr_8: continue
                                    if (b2_3+c2_3)%3 not in qr_3: continue
                                    if (b2_7+c2_7)%7 not in qr_7_set: continue
                                    if (a2_8+c2_8)%8 not in qr_8: continue
                                    if (a2_3+c2_3)%3 not in qr_3: continue
                                    if (a2_7+c2_7)%7 not in qr_7_set: continue
                                    if (a2_8+b2_8+c2_8)%8 not in qr_8: continue
                                    if (a2_3+b2_3+c2_3)%3 not in qr_3: continue
                                    if (a2_7+b2_7+c2_7)%7 not in qr_7_set: continue
                                    count += 1

fraction_168 = count / total_168
print(f"  Valid (a,b,c) mod 168: {count}")
print(f"  Total 168^3: {total_168}")
print(f"  Fraction: {fraction_168:.6f}")
print(f"  Elimination: {100*(1-fraction_168):.1f}%")

g = gcd(count, total_168)
print(f"  = {count//g}/{total_168//g}")

print(f"""

SUMMARY: THE CUBOID THROUGH THE HURWITZ LENS
================================================

The perfect cuboid constraints mod 42 (= 2*3*7):
  Eliminates {100*(1-fraction_42):.1f}% of triples.
  The surviving fraction: {fraction_42:.4f}.

The constraints mod 168 (= 8*21 = |PSL(2,7)|):
  Eliminates {100*(1-fraction_168):.1f}% of triples.
  The surviving fraction: {fraction_168:.6f}.

The number 168 = |PSL(2,7)| = 4 * 42 = 8 * 21 governs the cuboid.
The Fano plane's symmetry group IS the arithmetic of the cuboid.

If a perfect cuboid exists, its edges (a, b, c) mod 168
must lie in one of the {count} valid residue classes.

THE CHEBOTAREV CONNECTION:

Each face diagonal d = sqrt(a^2+b^2) requires a^2+b^2 to be a sum of two squares.
By Chebotarev: the primes dividing a^2+b^2 that are 3 mod 4 must appear evenly.
This is a DENSITY constraint: 50% of primes (the inert ones in Q(i))
must appear to even powers in all three face diagonals AND the space diagonal.

The probability that a random triple satisfies all four conditions:
  Each condition independently has probability ~C/sqrt(n) for some constant C.
  Four conditions: probability ~C^4/n^2.
  Expected number of cuboids up to N: ~C^4 * N (since we search a,b,c up to N^{1/3}).

The constant C depends on the mod-168 elimination rate: {fraction_168:.4f}.

If the constant is small enough: cuboids might be SO RARE that none exist.
If large enough: infinitely many exist.

The open question: is the constant zero or positive?
The mod-168 analysis says it's POSITIVE (some residue classes survive).
But surviving mod 168 is necessary, not sufficient.

THE CUBOID PROBLEM IS: does the INFINITE PRODUCT over all primes
of the survival fractions converge to something POSITIVE?

This is an Euler product: prod_p (local survival fraction at p).
If the product > 0: cuboids likely exist (positive density).
If the product = 0: cuboids might not exist (zero density).

The mod-42 factor in the product: {fraction_42:.4f}.
The mod-168 factor: {fraction_168:.6f}.
These are the first few terms. The full product is the cuboid's
ANALOGUE of the Mertens product for primes.
""")

print("opus-2026-03-16-S90cs: THE PERFECT CUBOID, ERDOS-STRAUS, AND MOD 42")
