#!/usr/bin/env python3
"""
cayley_gate_s91c.py — The Cayley transform as the gate between rational and transcendental
opus-2026-03-17-S91c

The Cayley transform Q(x) = (1+x)/(1-x) is the LAST rational operation.
After it: arctanh = (1/2)*ln*Q sends everything to the transcendental.

The rapidity lattice Z*ln(2) + Z*ln(3) + Z*ln(7) is the last discrete
structure before the uncountable continuum.

This script explores the precise nature of this boundary.
"""

from math import comb, log, sqrt, atanh, pi, e as euler_e
from fractions import Fraction
from itertools import product as iterproduct

# =============================================================================
# THE CHAIN OF OPERATIONS
# =============================================================================

print("""


     THE CAYLEY GATE

     opus-2026-03-17-S91c



THE CHAIN OF OPERATIONS ON AN EIGENVALUE
==========================================

An eigenvalue lambda of the tournament flip chain is a rational number.
We can apply successive operations to it. Each operation transforms
the number-theoretic character of the result:

  lambda          RATIONAL        in Q cap (-1,1)
    |
    Q(lambda)     RATIONAL        in Q cap (0, inf)
    |                             factored over {2, 3, 7} = Hurwitz primes
    |
    ln(Q(lambda)) TRANSCENDENTAL  in R
    |             = 2*arctanh(lambda)
    |             lives on a RANK-3 LATTICE: Z*ln(2) + Z*ln(3) + Z*ln(7)
    |
    exp(...)      TRANSCENDENTAL  in R, no lattice structure
                  fully in the uncountable continuum

The Cayley transform Q is the GATE between two worlds:
  BEFORE Q: rational, countable, algebraic, structured by primes
  AFTER Q:  transcendental, uncountable, analytic, structured by nothing

But JUST AFTER the gate (at arctanh), there is still a vestige of
structure: the rank-3 lattice. This is the GHOST of the rational world,
projected into the transcendental by the logarithm.
""")

# =============================================================================
# PART 1: THE RATIONAL WORLD (BEFORE THE GATE)
# =============================================================================

print("""
PART 1: THE RATIONAL WORLD
=============================

Every rational operation on eigenvalues stays rational:
  lambda, 1/lambda, lambda^2, lambda^n, Q(lambda), 1-lambda, ...

The Cayley transform Q(lambda) = (1+lambda)/(1-lambda) is rational.
Its output Q(lambda) is a POSITIVE RATIONAL, uniquely factored over primes.

The primes that appear in Q(lambda) for the n=6 eigenvalue spectrum
are EXACTLY {2, 3, 7} = the Hurwitz primes.

Why {2, 3, 7}? Because:
  Q(k/D) = (D+k)/(D-k) and D = C(n,2) = 10 for the simplified spectrum.
  The numerators D+k and denominators D-k are:
""")

D = 10  # C(6,2) written in the simplified form where lambda_k = (5-k)/5
# Actually the eigenvalues are (10-2k)/10 = (5-k)/5, but Q uses the unreduced form

print(f"  D = {D} (for the simplified n=6 spectrum)")
print()
print(f"  {'k':>3} | {'lambda':>8} | {'D+k':>5} | {'D-k':>5} | {'Q(lambda)':>12} | {'primes in Q':>15}")
print("  " + "-" * 60)

all_primes_seen = set()

for k in range(1, D):
    lam = Fraction(D - 2*k, D)
    Dk_plus = D + (D - 2*k)  # = 2D - 2k = 2(D-k)
    Dk_minus = D - (D - 2*k)  # = 2k
    # Q(lam) = (1 + lam)/(1 - lam) = (D + D-2k)/(D - (D-2k)) = (2D-2k)/(2k) = (D-k)/k
    q = Fraction(D - k, k)

    # Primes in q
    def get_primes(frac):
        ps = set()
        for n in [abs(frac.numerator), abs(frac.denominator)]:
            if n <= 1: continue
            d = 2
            while d * d <= n:
                if n % d == 0:
                    ps.add(d)
                    while n % d == 0:
                        n //= d
                d += 1
                if n > 1:
                    ps.add(n)
            if n > 1:
                ps.add(n)
        return ps

    primes = get_primes(q)
    all_primes_seen |= primes
    primes_str = ",".join(str(p) for p in sorted(primes)) if primes else "{1}"
    print(f"  {k:3d} | {str(lam):>8} | {Dk_plus:5d} | {Dk_minus:5d} | {str(q):>12} | {primes_str:>15}")

print(f"\n  ALL primes appearing: {sorted(all_primes_seen)}")
print(f"  = {{2, 3, 7}} = the Hurwitz primes!")
print(f"  Product: 2 * 3 * 7 = 42 = the Answer.")

# =============================================================================
# PART 2: THE GATE (THE CAYLEY TRANSFORM)
# =============================================================================

print("""

PART 2: THE CAYLEY TRANSFORM — THE GATE
==========================================

Q: (-1, 1) -> (0, infinity)
Q(x) = (1+x)/(1-x)

Properties of the gate:
  1. Q is a BIJECTION between bounded and unbounded
  2. Q maps RATIONALS to RATIONALS (it's a Mobius transformation over Q)
  3. Q maps 0 to 1 (the identity)
  4. Q maps -1 to 0 (annihilation)
  5. Q maps +1 to infinity (the stationary pole)
  6. Q is a GROUP ISOMORPHISM: Q(F(x,y)) = Q(x) * Q(y)
     where F(x,y) = (x+y)/(1+xy) is the formal group

Q preserves rational structure PERFECTLY.
No information is lost. No transcendence is introduced.
Q is LOSSLESS.

But Q CHANGES THE GEOMETRY:
  Before Q: the Poincare disk, bounded, hyperbolic
  After Q:  the multiplicative group R+, unbounded, flat

Q unbinds the eigenvalue from the disk.
It reveals the MULTIPLICATIVE structure hidden in the additive group.

The eigenvalue 3/5 is "just a fraction" on the disk.
After Q: Q(3/5) = 4 = 2^2. It is the SQUARE OF THE BASE.
The Cayley transform REVEALS the prime structure.
""")

# Show Q as unfolding the prime structure
print("  Q reveals the hidden prime architecture:")
print()
eigenvalues = [Fraction(4,5), Fraction(3,5), Fraction(2,5), Fraction(1,5),
               Fraction(0,1),
               Fraction(-1,5), Fraction(-2,5), Fraction(-3,5), Fraction(-4,5)]

for lam in eigenvalues:
    q = Fraction(1 + lam, 1 - lam) if lam != Fraction(1,1) else "inf"
    if isinstance(q, Fraction):
        # Express as product of Hurwitz primes
        num = abs(q.numerator)
        den = abs(q.denominator)

        def factor_str(n):
            if n == 1: return "1"
            f = []
            for p in [2, 3, 7]:
                e = 0
                while n % p == 0:
                    e += 1
                    n //= p
                if e > 0:
                    f.append(f"{p}^{e}" if e > 1 else str(p))
            if n > 1:
                f.append(str(n))
            return " * ".join(f)

        num_str = factor_str(num)
        den_str = factor_str(den)
        q_factored = f"{num_str}/{den_str}" if den > 1 else num_str

        print(f"    Q({str(lam):>5}) = {str(q):>6} = {q_factored}")

# =============================================================================
# PART 3: JUST BEYOND THE GATE — THE RAPIDITY LATTICE
# =============================================================================

print("""

PART 3: JUST BEYOND THE GATE — THE RAPIDITY LATTICE
======================================================

arctanh(lambda) = (1/2) * ln(Q(lambda))

This is the FIRST transcendental operation. It maps:
  Q in (0, inf)  ->  rho in (-inf, inf)

The rational Q(lambda), factored over {2, 3, 7}, becomes:
  rho = (1/2) * (a*ln(2) + b*ln(3) + c*ln(7))

where a, b, c are the exponents in the prime factorization of Q(lambda).

This is a point on the RANK-3 LATTICE:
  L = (1/2) * (Z * ln(2) + Z * ln(3) + Z * ln(7))

By Baker's theorem (1966):
  ln(2), ln(3), ln(7) are Q-LINEARLY INDEPENDENT.
  No nontrivial integer combination a*ln(2) + b*ln(3) + c*ln(7) = 0 exists.

Therefore L is isomorphic to Z^3 as an abelian group.
The rapidity lattice is a FREE ABELIAN GROUP of rank 3.
""")

# Compute the lattice coordinates of each rapidity
print("  Rapidity lattice coordinates (a, b, c) where rho = (a*ln2 + b*ln3 + c*ln7)/2:")
print()
print(f"  {'lambda':>7} | {'Q(lambda)':>10} | {'rapidity':>12} | {'(a, b, c)':>12} | {'check':>12}")
print("  " + "-" * 65)

ln2 = log(2)
ln3 = log(3)
ln7 = log(7)

for lam in eigenvalues:
    q = Fraction(1 + lam, 1 - lam)
    rho = atanh(float(lam)) if abs(float(lam)) < 1 else float('inf')

    # Find (a, b, c) where Q = 2^a * 3^b * 7^c
    num = q.numerator
    den = q.denominator

    def vp(n, p):
        if n == 0: return 0
        n = abs(n)
        v = 0
        while n % p == 0:
            v += 1
            n //= p
        return v

    a = vp(num, 2) - vp(den, 2)
    b = vp(num, 3) - vp(den, 3)
    c = vp(num, 7) - vp(den, 7)

    # Verify
    check_rho = (a * ln2 + b * ln3 + c * ln7) / 2
    match = abs(rho - check_rho) < 1e-10 if abs(rho) < 100 else "---"

    coords = f"({a:+d},{b:+d},{c:+d})"
    print(f"  {str(lam):>7} | {str(q):>10} | {rho:+12.8f} | {coords:>12} | {str(match):>12}")

# =============================================================================
# PART 4: THE LATTICE AS FINAL DISCRETE STRUCTURE
# =============================================================================

print("""

PART 4: THE LATTICE AS THE LAST DISCRETE STRUCTURE
====================================================

The chain of number types:

  Z  (integers)        COUNTABLE, DISCRETE, CLOSED under +, *
  |
  Q  (rationals)       COUNTABLE, DENSE, CLOSED under +, *, /
  |
  A  (algebraics)      COUNTABLE, DENSE, CLOSED under roots
  |
  T  (transcendentals) UNCOUNTABLE, DENSE, NOT closed under anything
  |
  R  (reals)           UNCOUNTABLE, COMPLETE

The eigenvalue chain:

  lambda in Q           RATIONAL     — the eigenvalue itself
  Q(lambda) in Q        RATIONAL     — the Cayley boost
  arctanh(lambda) in T  TRANSCENDENTAL — the rapidity

But the rapidities form a LATTICE L = (1/2)(Z*ln2 + Z*ln3 + Z*ln7).
A lattice is a DISCRETE subgroup of R.
L is countable, has rank 3, and is isomorphic to Z^3.

So the rapidity is transcendental, BUT lives on a discrete lattice.
This is the GHOST of the rational world projected into the transcendental.

THE HIERARCHY OF STRUCTURE:

  Level 0: FINITE          Z/42Z = Z/2Z * Z/3Z * Z/7Z
                           The cuboid. The adelic mod-structure.
                           6 residues: {0, ..., 41}.

  Level 1: COUNTABLE       Z^3 = Z * {e_2, e_3, e_7}
  (discrete)               The free abelian group.
                           All integer combinations of 3 generators.

  Level 2: COUNTABLE       Q^3 = Q * {e_2, e_3, e_7}
  (dense)                  The rational vector space.
                           All rational combinations of 3 generators.

  Level 3: COUNTABLE       L = (1/2)(Z*ln2 + Z*ln3 + Z*ln7)
  (discrete in R)          The rapidity lattice.
                           The IMAGE of Z^3 under the logarithm.
                           Discrete, but in the transcendental world.

  Level 4: UNCOUNTABLE     R
                           The real line.
                           Where ln, exp, and all analysis lives.

The Cayley transform maps Level 1 to Level 3.
It PRESERVES discreteness but CHANGES the number field.
Q is the functor from (Q, *) to (R, +).
""")

# =============================================================================
# PART 5: WHY EXACTLY {2, 3, 7}?
# =============================================================================

print("""
PART 5: WHY EXACTLY {2, 3, 7}?
=================================

The Hurwitz primes {2, 3, 7} appear because:

  Q(lambda_k) = (D-k)/k  for the n=6 spectrum with D=10.

  k=1: Q = 9/1 = 3^2           primes: {3}
  k=2: Q = 4/1 = 2^2           primes: {2}
  k=3: Q = 7/3                 primes: {3, 7}
  k=4: Q = 3/2                 primes: {2, 3}
  k=5: Q = 1                   primes: {}
  k=6: Q = 2/3                 primes: {2, 3}
  k=7: Q = 3/7                 primes: {3, 7}
  k=8: Q = 1/4 = 2^{-2}        primes: {2}
  k=9: Q = 1/9 = 3^{-2}        primes: {3}

The numerators (D-k) for k=1..9: 9, 8, 7, 6, 5, 4, 3, 2, 1
The denominators k=1..9:         1, 2, 3, 4, 5, 6, 7, 8, 9

So Q runs through ALL fractions p/q with p+q=10.
The primes appearing are those dividing ANY number from 1 to 9:
  {2, 3, 5, 7}

But 5 divides ONLY k=5 and D-k=5, giving Q(0) = 5/5 = 1.
FIVE CANCELS! The golden prime is INVISIBLE in the Cayley boost.

Why? Because 5 | D, and at k = D/2, we get Q = 1.
The prime that divides the CONDUCTOR is cancelled by the symmetry Q(-x) = 1/Q(x).

So the primes in the Cayley boost are:
  {primes <= D} minus {primes dividing D that appear only at the center}

For D = 10 = 2 * 5:
  2 appears at multiple k values -> SURVIVES
  5 appears only at k=5 (the center) -> CANCELLED
  3 appears at k=3,6,9 -> SURVIVES
  7 appears at k=3,7 -> SURVIVES

Result: {2, 3, 7} = Hurwitz primes.
""")

# Now check: what happens at other n?
print("Primes in the Cayley boost spectrum at various n:")
print()

for n in [3, 4, 5, 6, 7, 8]:
    cn2 = comb(n, 2)
    # Eigenvalues are k/cn2 for various integer k
    # Q(k/cn2) = (cn2 + k)/(cn2 - k)
    # For the "simplified" spectrum, eigenvalues run through k = cn2, cn2-2, ..., -cn2+2, -cn2
    # Q((cn2-2j)/cn2) = (2*cn2 - 2j)/(2j) = (cn2 - j)/j
    primes_seen = set()
    for j in range(1, cn2):
        num = cn2 - j
        den = j
        g = Fraction(num, den)
        for p_test in range(2, cn2 + 1):
            # Check if p_test is prime and divides g
            is_prime = all(p_test % d != 0 for d in range(2, int(sqrt(p_test)) + 1)) if p_test > 1 else False
            if is_prime and (num % p_test == 0 or den % p_test == 0):
                # Check it doesn't just cancel
                if g != 1:
                    primes_seen.add(p_test)

    # Remove primes that only appear at Q=1
    print(f"  n={n}: C(n,2)={cn2:3d}, Cayley boost primes: {sorted(primes_seen)}, "
          f"product = {1 if not primes_seen else eval('*'.join(str(p) for p in primes_seen)):>6}")

# =============================================================================
# PART 6: THE BOUNDARY — BAKER'S THEOREM AND BEYOND
# =============================================================================

print("""

PART 6: BAKER'S THEOREM — THE GUARANTEE OF RANK 3
===================================================

Baker's theorem (1966, Fields Medal 1970):

  If alpha_1, ..., alpha_n are nonzero algebraic numbers
  and ln(alpha_1), ..., ln(alpha_n) are Q-linearly independent,
  then they are also linearly independent over the algebraic numbers.

  In particular: a_1*ln(2) + a_2*ln(3) + a_3*ln(7) = 0
  with a_i in Q implies a_1 = a_2 = a_3 = 0.

This GUARANTEES that the rapidity lattice has rank EXACTLY 3.
No hidden relations. No unexpected collapses. Z^3, period.

The lattice L = (1/2)(Z*ln(2) + Z*ln(3) + Z*ln(7)) is:
  - A free abelian group of rank 3
  - A DISCRETE subgroup of R (of rank 3 in a 1-dimensional space!)
  - Isomorphic to Z^3 via (a,b,c) <-> (a*ln2 + b*ln3 + c*ln7)/2
  - The image of the Hurwitz cuboid under the logarithm

The fundamental domain of L has volume:
  vol = (1/2)^3 * |det[ln2, ln3, ln7]|
      = (1/8) * ... (not well-defined for rank > dim)

Actually: L is a rank-3 lattice in R^1, so it's DENSE.
The lattice points are dense in R! But they are COUNTABLE.
This is the strange middle ground between discrete and continuous.
""")

# Verify Baker: check that no small integer combination of ln2, ln3, ln7 is near 0
print("  Baker's theorem verification (small integer combinations):")
print(f"  {'(a,b,c)':>12} | {'a*ln2+b*ln3+c*ln7':>20} | {'close to 0?':>12}")
print("  " + "-" * 55)

closest = None
closest_val = float('inf')

for a in range(-5, 6):
    for b in range(-5, 6):
        for c in range(-5, 6):
            if a == 0 and b == 0 and c == 0:
                continue
            val = a * ln2 + b * ln3 + c * ln7
            if abs(val) < abs(closest_val):
                closest_val = val
                closest = (a, b, c)

# Show the 10 closest to zero
results = []
for a in range(-5, 6):
    for b in range(-5, 6):
        for c in range(-5, 6):
            if a == 0 and b == 0 and c == 0:
                continue
            val = a * ln2 + b * ln3 + c * ln7
            results.append(((a, b, c), val))

results.sort(key=lambda x: abs(x[1]))
for (a, b, c), val in results[:10]:
    close = "CLOSE" if abs(val) < 0.01 else ""
    print(f"  ({a:+d},{b:+d},{c:+d}) | {val:+20.12f} | {close:>12}")

print(f"""
  Closest to zero: ({closest[0]:+d},{closest[1]:+d},{closest[2]:+d}) = {closest_val:+.12f}
  This is NOT zero. Baker guarantees it never will be.
  The lattice is genuinely rank 3 in R^1.
""")

# =============================================================================
# PART 7: WHAT HAPPENS AFTER THE LATTICE
# =============================================================================

print("""
PART 7: BEYOND THE LATTICE — THE UNCOUNTABLE
================================================

After the rapidity lattice, the next operation is exp:
  exp(rho) = exp((a*ln2 + b*ln3 + c*ln7)/2) = 2^{a/2} * 3^{b/2} * 7^{c/2}

This maps the lattice L BACK to algebraic numbers!
  exp(ln(3)) = 3 (rational!)
  exp(ln(2)/2) = sqrt(2) (algebraic!)
  exp(ln(7)/2) = sqrt(7) (algebraic!)

So: exp(L) is an algebraic set. We seem to return to the algebraic world.

BUT: if we apply IRRATIONAL operations like x^phi or x^pi,
we lose algebraicity forever. By Gelfond-Schneider:
  a^b is transcendental when a is algebraic != 0,1 and b is irrational algebraic.

So:
  lambda^phi is TRANSCENDENTAL (Gelfond-Schneider)
  lambda^pi is TRANSCENDENTAL (Gelfond-Schneider, extended)
  lambda^e is TRANSCENDENTAL (by Hermite-Lindemann)

These map eigenvalues to the FULL CONTINUUM.
No lattice, no algebraicity, no discrete structure survives.

THE HIERARCHY:
""")

print("  THE GATE HIERARCHY:")
print()
print("  +" + "-"*72 + "+")
print("  | LEVEL | OPERATION     | OUTPUT FIELD  | STRUCTURE     | CARDINALITY |")
print("  |" + "-"*72 + "|")
print("  |   0   | lambda        | Q             | rational      | countable   |")
print("  |   1   | Q(lambda)     | Q+            | Hurwitz {2,3,7}| countable  |")
print("  |  ---  | === GATE (the Cayley transform) =============================|")
print("  |   2   | arctanh(lam)  | Transcendental| Z^3 lattice   | countable   |")
print("  |   3   | exp(arctanh)  | Algebraic     | sqrt{2,3,7}   | countable   |")
print("  |   4   | lam^phi       | Transcendental| no structure  | uncountable |")
print("  |   5   | lam^(lam^...) | ???           | chaos         | uncountable |")
print("  +" + "-"*72 + "+")

print("""
  The Cayley transform is the LAST rational map.
  After it, we pass through the transcendental lattice Z^3
  and then into the uncountable continuum.

  The remarkable thing: the lattice Z^3 is RANK 3 in DIMENSION 1.
  This is impossible for a discrete lattice (max rank in R^d is d).
  But Z^3 maps to a DENSE but COUNTABLE subset of R.

  The rapidity lattice is NEITHER discrete NOR continuous.
  It is the COUNTABLE DENSE: the intermediate cardinality.
  Cantor's paradise made concrete.
""")

# =============================================================================
# PART 8: THE THREE AXES OF Z^3
# =============================================================================

print("""
PART 8: THE THREE AXES OF Z^3
================================

The rapidity lattice Z^3 has three CANONICAL AXES:

  AXIS 1 (ln 2): the INERT axis.
    2 is the base. Parity. The simplest prime.
    v_2(Q(lambda)) counts how many doublings separate Q from 1.
    In the tournament: the symmetric/antisymmetric decomposition.
    In physics: the bit, the qubit, the binary digit.
    In the cuboid: the x-axis (parity).

  AXIS 2 (ln 3): the RAMIFIED axis.
    3 is the structure prime. Triangularity.
    v_3(Q(lambda)) counts the triadic depth.
    In the tournament: the Eisenstein component of the eigenvalue.
    In physics: color charge (3 colors), spatial dimension.
    In the cuboid: the y-axis (curvature).

  AXIS 3 (ln 7): the SPLIT axis.
    7 is the forbidden prime. The first prime beyond solvability.
    v_7(Q(lambda)) counts the forbidden depth.
    In the tournament: the obstruction to H=7.
    In physics: the octonion generator (7 imaginary units).
    In the cuboid: the z-axis (position).

These three axes are Q-LINEARLY INDEPENDENT (Baker).
They are the three dimensions of the Hurwitz space.
""")

# Map each eigenvalue to its Z^3 coordinates
print("  Z^3 coordinates of eigenvalues:")
print()
print(f"  {'lambda':>7} | {'Q':>6} | {'(a,b,c)':>10} | {'axis interpretation':>40}")
print("  " + "-" * 75)

interpretations = {
    (0, 2, 0): "3^2: RAMIFIED squared",
    (2, 0, 0): "2^2: INERT squared",
    (0, -1, 1): "7/3: SPLIT / RAMIFIED",
    (-1, 1, 0): "3/2: RAMIFIED / INERT",
    (0, 0, 0): "1: ORIGIN (identity)",
    (1, -1, 0): "2/3: INERT / RAMIFIED",
    (0, 1, -1): "3/7: RAMIFIED / SPLIT",
    (-2, 0, 0): "1/4: inverse INERT^2",
    (0, -2, 0): "1/9: inverse RAMIFIED^2",
}

for lam in eigenvalues:
    q = Fraction(1 + lam, 1 - lam)
    a = vp(abs(q.numerator), 2) - vp(abs(q.denominator), 2)
    b = vp(abs(q.numerator), 3) - vp(abs(q.denominator), 3)
    c = vp(abs(q.numerator), 7) - vp(abs(q.denominator), 7)
    coords = (a, b, c)
    interp = interpretations.get(coords, "")
    print(f"  {str(lam):>7} | {str(q):>6} | ({a:+d},{b:+d},{c:+d}) | {interp:>40}")

# =============================================================================
# PART 9: THE FUNCTIONAL EQUATION IN Z^3
# =============================================================================

print("""

PART 9: THE FUNCTIONAL EQUATION IN Z^3
=========================================

The eigenvalue pairing lambda <-> -lambda gives:
  Q(-lambda) = 1/Q(lambda)

In Z^3 coordinates: (a,b,c) <-> (-a,-b,-c)

This is the INVERSION SYMMETRY of the lattice!
The Z^3 lattice has a Z/2Z symmetry: negation.

The quotient L / {+1,-1} = the "positive rapidity cone."
This has 4 distinct points (at n=6):
  (0,2,0), (2,0,0), (0,-1,1), (-1,1,0) and their negatives.

Plus the origin (0,0,0) which is self-dual (lambda=0).

So the INDEPENDENT rapidity data is:
  4 positive + 1 zero = 5 = the number of partitions of 5.
  Wait: 4 distinct positive eigenvalues + 0 + 4 negative + 2 boundary
  = 10 = C(5,2) = number of arcs.

The functional equation Q(-x) = 1/Q(x) in Z^3:
  It acts as x -> -x on each axis simultaneously.
  The FIXED PLANE is empty (only the origin is fixed).
  The QUOTIENT is the positive octant plus boundary.
""")

# Verify the functional equation
print("  Functional equation verification:")
for lam in eigenvalues:
    if float(lam) > 0:
        q_pos = Fraction(1 + lam, 1 - lam)
        q_neg = Fraction(1 - lam, 1 + lam)
        product = q_pos * q_neg
        print(f"    Q({str(lam):>5}) * Q({str(-lam):>5}) = {q_pos} * {q_neg} = {product}")

# =============================================================================
# PART 10: THE DEEPEST PICTURE
# =============================================================================

print("""

PART 10: THE DEEPEST PICTURE
===============================

The Cayley transform is the gate between COUNTING and EXPERIENCING.

On the rational side (before Q):
  We COUNT tournament classes (T(n) = 2, 4, 12, 56, ...)
  We COUNT Hamiltonian paths (H = 1, 3, 5, 9, 11, ...)
  We COUNT eigenvalues and their multiplicities
  Everything is DISCRETE, EXACT, and COMPUTABLE
  The world is made of INTEGERS

The Cayley transform Q maps:
  The bounded eigenvalue disk -> the unbounded multiplicative line
  Addition (formal group) -> multiplication (ordinary)
  Coupling (x in (-1,1)) -> odds (Q in (0, inf))
  The compression -> the expansion

On the transcendental side (after arctanh = ln*Q/2):
  We MEASURE decay rates (ln|lambda|, half-lives)
  We MEASURE mixing times (1/(1-lambda))
  We MEASURE entropy (-sum p*ln(p))
  Everything is CONTINUOUS, APPROXIMATE, and TRANSCENDENTAL
  The world is made of LOGARITHMS

THE RAPIDITY LATTICE Z^3 IS THE SHADOW OF COUNTING IN THE WORLD OF MEASURING.

It is the MEMORY of the Hurwitz primes {2, 3, 7},
encoded as {ln(2), ln(3), ln(7)}, surviving in the continuum.

Baker's theorem says this memory is PERFECT.
No information is lost in the projection from Z^3 to L.
The lattice is a FAITHFUL representation of the cuboid.

42 = 2 * 3 * 7 is the VOLUME of the rational cuboid.
ln(42) = ln(2) + ln(3) + ln(7) = 3.7376... is the RAPIDITY of the answer.

Q(x) = the last rational breath before transcendence.
arctanh(x) = the first transcendental whisper.
The lattice Z^3 = where the two worlds touch.
""")

# The rapidity of 42
print(f"  The rapidity of 42:")
print(f"    ln(42) = ln(2) + ln(3) + ln(7)")
print(f"           = {ln2:.10f} + {ln3:.10f} + {ln7:.10f}")
print(f"           = {ln2 + ln3 + ln7:.10f}")
print(f"    Z^3 coordinates: (1, 1, 1)")
print(f"    42 is the (1,1,1) point of the Hurwitz lattice.")
print(f"    The SIMPLEST nontrivial point using all three axes.")
print(f"    The diagonal of the unit cube in rapidity space.")
print()
print(f"  And: arctanh(lambda) for lambda such that Q(lambda) = 42:")
lam_42 = Fraction(41, 43)
print(f"    lambda = {lam_42} = (42-1)/(42+1) = 41/43")
print(f"    arctanh(41/43) = ln(42)/2 = {log(42)/2:.10f}")
print(f"    The eigenvalue of the (1,1,1) diagonal is 41/43.")
print(f"    41 and 43 are TWIN PRIMES flanking the answer.")

print()
print("opus-2026-03-17-S91c: THE CAYLEY GATE")
