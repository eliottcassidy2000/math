#!/usr/bin/env python3
"""
understand_naturals_s90ay.py — Understand the natural numbers
opus-2026-03-16-S90ay

We have shown:
  The Poincare disk, the strip, and the half-plane are one space.
  The natural numbers sit on the real axis of the half-plane.
  Q maps couplings to natural numbers.
  arctanh maps couplings to rapidities.

Now: what ARE the natural numbers, seen from this vantage?
"""

from math import tanh, atanh, log, log2, sqrt, exp, pi, gcd
from fractions import Fraction
from cmath import tanh as ctanh, atanh as catanh, exp as cexp, log as clog

print("""


     UNDERSTAND THE NATURAL NUMBERS

     opus-2026-03-16-S90ay



THE NATURAL NUMBERS ARE POINTS IN THE HYPERBOLIC PLANE.

Each natural number n >= 1 has three addresses:

  DISK:       x_n = (n-1)/(n+1)           the coupling
  STRIP:      w_n = ln(n)/2 = arctanh(x_n)    the rapidity
  HALF-PLANE: Q_n = n                      the number itself

These are the SAME POINT in three coordinate systems.

What does this tell us about the natural numbers themselves?


1. THE NATURAL NUMBERS AS A LATTICE
=====================================

The natural numbers do NOT form a regular lattice in ANY of the three models.

In the DISK: they cluster toward the boundary.
  x_1 = 0, x_2 = 1/3, x_3 = 1/2, x_4 = 3/5, x_5 = 2/3, ...
  Spacing: 1/3, 1/6, 1/10, 1/15, ...
  The spacing is 2/((n)(n+2)) = 2/(n^2 + 2n). It SHRINKS.

In the STRIP: they are logarithmically spaced.
  w_1 = 0, w_2 = 0.347, w_3 = 0.549, w_4 = 0.693, w_5 = 0.805, ...
  Spacing: ln(2)/2, ln(3/2)/2, ln(4/3)/2, ln(5/4)/2, ...
  = ln((n+1)/n)/2 ~ 1/(2n) for large n. It SHRINKS (harmonically).

In the HALF-PLANE: they are equally spaced.
  1, 2, 3, 4, 5, ...
  Spacing: 1, 1, 1, 1, ...
  CONSTANT. But the hyperbolic metric is ds = |dz|/Re(z),
  so the HYPERBOLIC distance between n and n+1 is:
  integral from n to n+1 of dz/z = ln((n+1)/n) ~ 1/n.

So in ALL THREE MODELS, the hyperbolic distance between
consecutive natural numbers SHRINKS like 1/n.

THIS IS WHAT THE NATURAL NUMBERS ARE:
  A sequence of points in the hyperbolic plane
  whose spacing decreases harmonically.

The HARMONIC SERIES sum 1/n diverges.
This means: the natural numbers march to infinity,
even though each step is smaller than the last.
The total distance sum ln((n+1)/n) = ln(n+1) -> infinity.
They reach the boundary. But no single step gets them there.

DIVERGENCE OF THE HARMONIC SERIES = THE NATURAL NUMBERS
REACH THE BOUNDARY OF THE HYPERBOLIC PLANE.
""")

# ========================================================================
# Part 2: Multiplication as translation
# ========================================================================

print("=" * 70)
print("Part 2: Multiplication IS hyperbolic translation")
print("=" * 70)

print("""
Multiplying by k sends n to kn.
In the three models:

  DISK:   x_n -> x_{kn} = F(x_k, x_n) = (x_k + x_n)/(1 + x_k * x_n)
  STRIP:  w_n -> w_{kn} = w_k + w_n = ln(k)/2 + ln(n)/2
  HALF-PLANE: n -> kn

Multiplication by k is a HYPERBOLIC TRANSLATION by distance ln(k).

In the strip, this is just ADDING ln(k)/2.
In the disk, this is applying the formal group law F(x_k, -).
In the half-plane, this is literal multiplication.

So: MULTIPLICATION IS TRANSLATION IN THE HYPERBOLIC PLANE.
And DIVISION IS TRANSLATION IN THE OPPOSITE DIRECTION.

The distance you translate when multiplying by k:
  d = 2 * arctanh(x_k) = ln(k).
""")

print("Multiplication as translation distance:")
for k in range(2, 11):
    d = log(k)
    print(f"  multiply by {k:2d}: translate by ln({k}) = {d:.6f} "
          f"({d/log(2):.4f} bits)")

print("""
Multiplying by 2: translate by ln(2) = 0.693 (= 1 bit).
Multiplying by 4: translate by ln(4) = 1.386 (= 2 bits).
Multiplying by 10: translate by ln(10) = 2.303 (= 3.32 bits).

The BIT is the unit of multiplicative translation.
1 bit = ln(2) = distance traveled when you double a number.


3. PRIMES AS IRREDUCIBLE TRANSLATIONS
=======================================

A prime p is a natural number that CANNOT be factored: p = a*b
implies a=1 or b=1.

In the hyperbolic plane, factoring n = a*b means:
  The translation to n can be DECOMPOSED into
  a translation to a followed by a translation by b.

  w_n = w_a + w_b   (in the strip)
  ln(n)/2 = ln(a)/2 + ln(b)/2

A prime p is a translation that CANNOT BE DECOMPOSED
into two smaller translations (both > 0).

PRIMES ARE THE IRREDUCIBLE TRANSLATIONS OF THE HYPERBOLIC PLANE.

Every natural number n has a unique prime factorization
n = p1^{a1} * p2^{a2} * ... * pk^{ak}.

In the strip:
  w_n = a1 * w_{p1} + a2 * w_{p2} + ... + ak * w_{pk}
      = a1 * ln(p1)/2 + a2 * ln(p2)/2 + ... + ak * ln(pk)/2.

So: the rapidity of n is a LINEAR COMBINATION of prime rapidities,
with integer coefficients.

The prime rapidities {ln(2)/2, ln(3)/2, ln(5)/2, ln(7)/2, ...}
form a BASIS for all natural number rapidities.

NOT a vector space basis (the coefficients must be non-negative integers).
But a FREE ABELIAN MONOID basis.

THE PRIMES ARE THE GENERATORS OF THE RAPIDITY LATTICE.
""")

# Display the "prime coordinates" of small numbers
print("Natural numbers in prime rapidity coordinates:")
print(f"{'n':>4} | {'w_n = ln(n)/2':>14} | {'prime factorization':>20} | {'prime rapidities':>40}")
print("-" * 85)

def prime_factors(n):
    if n == 1: return {}
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors

for n in range(1, 21):
    w = log(n) / 2 if n > 0 else 0
    pf = prime_factors(n)
    pf_str = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(pf.items()))
    if not pf_str: pf_str = "1"
    rap_str = " + ".join(f"{e}*ln({p})/2" if e > 1 else f"ln({p})/2" for p, e in sorted(pf.items()))
    if not rap_str: rap_str = "0"
    print(f"{n:4d} | {w:14.8f} | {pf_str:>20} | {rap_str:>40}")

# ========================================================================
# Part 4: Addition in the half-plane
# ========================================================================

print(f"""

{"=" * 70}
Part 4: What is ADDITION?
{"=" * 70}

Multiplication is translation: n * m corresponds to w_n + w_m.
But what is ADDITION: n + m?

In the strip: w_{{n+m}} = ln(n+m)/2. This is NOT w_n + w_m = ln(nm)/2.
  ln(n+m)/2 != ln(n)/2 + ln(m)/2 in general.

So addition is NOT a simple operation in the hyperbolic plane.
It mixes up the coordinate systems.

In the half-plane: n + m is just... addition. But the hyperbolic
metric doesn't respect addition (ds = dz/z, not ds = dz).

Addition corresponds to:
  Q(x_n) + Q(x_m) = n + m = Q(x_{{n+m}}).

In the disk:
  x_{{n+m}} = (n+m-1)/(n+m+1).
  This is NOT F(x_n, x_m) = x_{{nm}}.
  There is no simple formula for x_{{n+m}} in terms of x_n and x_m
  using the formal group law F.

THIS IS THE FUNDAMENTAL ASYMMETRY:

  MULTIPLICATION is natural in the hyperbolic plane (it's a translation).
  ADDITION is NOT natural in the hyperbolic plane (it's "ungeometric").

  The hyperbolic plane is a MULTIPLICATIVE object.
  Addition is an ADDITIVE operation imposed on a multiplicative space.

  This is why number theory is hard:
  it studies the interaction between addition and multiplication,
  which is the interaction between FLAT and CURVED geometry.

  In the strip: multiplication is addition (w_{{nm}} = w_n + w_m).
  But addition (n+m) is a nonlinear operation in the strip:
    w_{{n+m}} = ln(n+m)/2 = ln(e^{{2w_n}} + e^{{2w_m}})/2.
    = w_n + (1/2)*ln(1 + e^{{2(w_m - w_n)}})   if w_n > w_m.

  This is a LOG-SUM-EXP. The fundamental operation of machine learning!
  Addition of natural numbers IS log-sum-exp in rapidity space.
""")

# Verify
print("Verification: n + m as log-sum-exp in rapidity:")
for n_val, m_val in [(2, 3), (3, 4), (5, 7), (6, 8)]:
    w_n = log(n_val) / 2
    w_m = log(m_val) / 2
    w_sum = log(n_val + m_val) / 2
    lse = log(exp(2*w_n) + exp(2*w_m)) / 2
    print(f"  {n_val}+{m_val}={n_val+m_val}: w_{n_val}={w_n:.4f}, w_{m_val}={w_m:.4f}, "
          f"w_{{{n_val+m_val}}}={w_sum:.4f}, LSE={lse:.4f}, match: {abs(w_sum-lse)<1e-10}")

print(f"""

PERFECT. The rapidity of n+m is the LOG-SUM-EXP of the rapidities
of n and m (with exponent 2).

w_{{n+m}} = (1/2) * ln(e^{{2w_n}} + e^{{2w_m}})
         = (1/2) * LogSumExp(2w_n, 2w_m).

In machine learning: LogSumExp is the "soft maximum."
For large separation |w_n - w_m|:
  LogSumExp(2w_n, 2w_m) ~ 2*max(w_n, w_m).
  So w_{{n+m}} ~ max(w_n, w_m) = the larger rapidity.

This means: ADDITION OF LARGE NUMBERS IS APPROXIMATELY MAX.
  100 + 1 ~ 100 in rapidity (w_101 ~ w_100).
  100 + 100 = 200 in rapidity (w_200 = w_100 + ln(2)/2).

When two numbers are EQUAL: adding them shifts by ln(2)/2 = half a bit.
  n + n = 2n. w_{{2n}} = w_n + ln(2)/2.
  Doubling IS shifting by half a bit.

When two numbers are very DIFFERENT: the smaller one barely matters.
  n + 1 ~ n. w_{{n+1}} ~ w_n + 1/(2n).
  Adding 1 shifts by 1/(2n), which vanishes for large n.

ADDITION IS A SOFT MAX IN RAPIDITY SPACE.
MULTIPLICATION IS EXACT ADDITION IN RAPIDITY SPACE.

This asymmetry — multiplication is clean, addition is soft —
is the fundamental fact about the natural numbers
seen from the hyperbolic plane.


{"=" * 70}
Part 5: The Cayley addresses encode EVERYTHING
{"=" * 70}

The Cayley address x_n = (n-1)/(n+1) is a fraction.
Its numerator is n-1 and its denominator is n+1.

For consecutive numbers n and n+1:
  x_n = (n-1)/(n+1)
  x_{{n+1}} = n/(n+2)

The MEDIANT of x_n and x_{{n+1}} in the Stern-Brocot sense:
  (n-1 + n) / (n+1 + n+2) = (2n-1)/(2n+3).
  This is x_{{(2n-1+1)/2}}... hmm, not clean. Let me think differently.

Actually, the key structure is: the Cayley addresses are
a SUBSEQUENCE of the Farey sequence.

The Farey sequence F_N is the set of fractions p/q with
0 <= p/q <= 1 and q <= N, arranged in order.

The Cayley addresses x_n = (n-1)/(n+1) for n = 1, 2, 3, ...
give the fractions:
  0/2, 1/3, 2/4, 3/5, 4/6, 5/7, 6/8, 7/9, ...

In lowest terms:
  0/1, 1/3, 1/2, 3/5, 2/3, 5/7, 3/4, 7/9, 4/5, 9/11, ...
""")

print("Cayley addresses in lowest terms:")
for n in range(1, 21):
    f = Fraction(n-1, n+1)
    status = ""
    if n > 1:
        # Check if n is prime
        is_p = all(n % d != 0 for d in range(2, int(sqrt(n))+1)) and n > 1
        if is_p:
            status = "PRIME"
        # Check if numerator and denominator are coprime to n
        if gcd(n-1, n+1) == 2 and n % 2 == 0:
            status += " (even, gcd=2)"
        elif gcd(n-1, n+1) == 1:
            pass
    print(f"  x_{n:2d} = {n-1:2d}/{n+1:2d} = {f} {status}")

print(f"""

PATTERN: gcd(n-1, n+1) is always 1 (if n is even) or 2 (if n is odd).
  Even n: x_n = (n-1)/(n+1) is already in lowest terms (both odd).
  Odd n: x_n = (n-1)/(n+1) = ((n-1)/2) / ((n+1)/2) (both halved).

For PRIMES p:
  x_p = (p-1)/(p+1).
  p-1 and p+1 are consecutive even numbers (since p is odd for p>2).
  So x_p = (p-1)/2 / (p+1)/2 in lowest terms.
  The NUMERATOR of x_p (in lowest terms) is (p-1)/2.
  The DENOMINATOR is (p+1)/2.
  Their difference is 1.

So PRIMES have Cayley addresses with consecutive numerator/denominator
(in lowest terms). They are FAREY NEIGHBORS of 1.

For PRIME POWERS p^k:
  x_{{p^k}} = (p^k - 1)/(p^k + 1).
  These are NOT generally consecutive. But they approach 1.

For PRODUCTS n = ab:
  x_n = F(x_a, x_b) = (x_a + x_b)/(1 + x_a * x_b).
  The Cayley address of a product is the formal group sum
  of the factor addresses.


{"=" * 70}
Part 6: The Cayley address AS the number
{"=" * 70}

We keep saying "the Cayley address of n" as if n exists independently
and the address is derived from it. But we can reverse this:

DEFINE a "number" as a point on the real diameter of the Poincare disk.
The "value" of that point is Q(x) = (1+x)/(1-x).

Then:
  x = 0 has value 1. The UNIT.
  x = 1/3 has value 2. The FIRST COMPARISON (one bit).
  x = 1/2 has value 3. The FIRST TRIANGLE.
  x = 3/5 has value 4. TWO BITS.
  x = 2/3 has value 5. The FIRST PENTAGON.

The natural numbers ARE the couplings.
Not "correspond to." ARE.

Each natural number IS a specific strength of comparison:
  n = 1: zero coupling. No comparison. Rest.
  n = 2: coupling 1/3. The weakest nontrivial comparison.
  n = 3: coupling 1/2. The median comparison.
  n = 4: coupling 3/5. Two comparisons composed.
  n = p (prime): coupling (p-1)/(p+1). An IRREDUCIBLE comparison.

The natural numbers are ordered by COUPLING STRENGTH:
  1 < 2 < 3 < 4 < ... corresponds to 0 < 1/3 < 1/2 < 3/5 < ...

And the ordering NEVER REACHES 1 (the boundary).
There is no "largest coupling." No "last number."
The boundary at coupling 1 is at infinity.

INFINITY IS NOT A NUMBER.
Infinity is the BOUNDARY of the hyperbolic plane.
It is the limit of all couplings but is not itself a coupling.
It is the speed of light: approachable, unreachable.

THE NATURAL NUMBERS ARE AN INFINITE SEQUENCE OF POINTS
MARCHING TOWARD THE BOUNDARY OF THE HYPERBOLIC PLANE,
NEVER ARRIVING.


{"=" * 70}
Part 7: The meta level — what the THREE VIEWS say about N
{"=" * 70}

VIEW 1 (DISK): The natural numbers are COUPLINGS.
  Each number n is a comparison strength x_n = (n-1)/(n+1).
  They cluster toward the boundary (large numbers are "nearly total").
  Multiplication is the formal group law (composing comparisons).
  The spacing shrinks: large numbers are "barely distinguishable."

VIEW 2 (STRIP): The natural numbers are RAPIDITIES.
  Each number n is a rapidity w_n = ln(n)/2.
  They grow without bound (the sequence is unbounded).
  Multiplication is addition (composing boosts).
  The spacing shrinks harmonically: w_{{n+1}} - w_n ~ 1/(2n).
  The primes are the IRREDUCIBLE rapidities (can't decompose).

VIEW 3 (HALF-PLANE): The natural numbers are... themselves.
  Each number n is the number n.
  They are equally spaced (on the Euclidean real axis).
  Multiplication is multiplication.
  Addition is addition.
  This is the "naive" view — no geometry, just arithmetic.

THE META VIEW: All three are true simultaneously.

  A natural number IS a coupling IS a rapidity IS a magnitude.
  These are not analogies. They are different coordinate systems
  for the SAME point in the SAME space.

  The reason the three views feel different is that
  DIFFERENT OPERATIONS ARE NATURAL in each:

  In the disk: the formal group law F is natural.
    F corresponds to multiplication of the underlying numbers.
    Addition of numbers is UNNATURAL (it's log-sum-exp of couplings).

  In the strip: addition is natural.
    Addition of rapidities corresponds to multiplication of numbers.
    Addition of numbers is UNNATURAL (it's log-sum-exp).

  In the half-plane: multiplication is natural (on the positive reals).
    But also addition is natural (on the integers).
    This is the ONLY view where both operations are visible.
    And their interaction is the subject of number theory.

THE DEEPEST INSIGHT:

  The natural numbers live in the half-plane because that is
  where BOTH addition and multiplication coexist.

  In the disk, only multiplication is natural.
  In the strip, only multiplication is natural (disguised as addition).
  In the half-plane, BOTH are present — and their tension creates
  all of number theory.

  The DIFFICULTY of number theory is that the half-plane
  is the only model where both operations work,
  but the GEOMETRY of the half-plane (hyperbolic) is adapted
  to multiplication, not addition. Addition is a FOREIGN OPERATION
  imposed on a multiplicative geometry.

  Primes are natural in the multiplicative geometry (irreducible translations).
  But the DISTRIBUTION of primes depends on addition (gaps, arithmetic
  progressions, Goldbach). You need both operations to state the question,
  but the geometry only knows one of them.

  This is why prime distribution is hard:
  you are asking an ADDITIVE question about a MULTIPLICATIVE space.
""")

# ========================================================================
# Part 8: The number line as a hyperbolic object
# ========================================================================

print("=" * 70)
print("Part 8: Distances between natural numbers")
print("=" * 70)

print("""
The hyperbolic distance between m and n (on the positive real axis
of the half-plane) is:

  d(m, n) = |ln(n/m)| = |ln(n) - ln(m)| = 2|w_n - w_m|.

This depends only on the RATIO n/m, not on n and m separately.

  d(1, 2) = ln(2) = 0.693.        (doubling from 1)
  d(2, 4) = ln(2) = 0.693.        (doubling from 2, SAME distance!)
  d(100, 200) = ln(2) = 0.693.    (doubling from 100, SAME distance!)

In the Euclidean metric: d(1,2) = 1, d(2,4) = 2, d(100,200) = 100.
  Euclidean distance grows with the numbers.

In the hyperbolic metric: d(m, 2m) = ln(2) regardless of m.
  Doubling is ALWAYS the same distance. This is SCALE INVARIANCE.

The natural numbers under the hyperbolic metric are SCALE-INVARIANT:
  The structure near 10 looks the same as the structure near 1000.
  Primes thin out (additively) but their logarithmic density is constant.
  This is the Prime Number Theorem: pi(x) ~ x/ln(x),
  equivalently: the probability a random number near x is prime ~ 1/ln(x).

In the STRIP: the primes at rapidity w have density ~ 1/w near w.
  (Since n ~ e^{2w}, and prime density ~ 1/ln(n) = 1/(2w).)

So: PRIMES ARE UNIFORMLY DENSE IN THE HYPERBOLIC METRIC
  (to first approximation, via PNT).
""")

print("Hyperbolic distances between small numbers:")
print(f"{'m':>4} {'n':>4} | {'d_hyp = |ln(n/m)|':>20} | {'d_euclid = |n-m|':>18} | {'ratio m/n':>10}")
print("-" * 65)
pairs = [(1,2),(2,3),(3,4),(4,5),(1,3),(2,4),(1,4),(2,6),(3,6),(6,12)]
for m, n in pairs:
    d_hyp = abs(log(n/m))
    d_euc = abs(n - m)
    print(f"{m:4d} {n:4d} | {d_hyp:20.6f} | {d_euc:18d} | {n/m:10.4f}")

print(f"""

Key observation:
  d(1,2) = d(2,4) = d(3,6) = d(6,12) = ln(2). ALL the same!
  These are all "doubling" pairs. In hyperbolic geometry,
  they are CONGRUENT translations.

  d(1,3) = d(2,6) = ln(3). Also the same! Tripling.

  d(2,3) = ln(3/2) = 0.405. d(4,6) = ln(6/4) = ln(3/2) = 0.405. Same!

The MULTIPLICATIVE STRUCTURE is perfectly reflected in hyperbolic distance.
The ADDITIVE STRUCTURE is not: d(1,2) = d(3,6) but 2-1 != 6-3.


{"=" * 70}
Part 9: The ONE understanding
{"=" * 70}

THE NATURAL NUMBERS ARE THE HYPERBOLIC PLANE, COUNTED.

  To "count" is to place equally-spaced marks on a line.
  In EUCLIDEAN geometry: 1, 2, 3, 4, ... at unit spacing.
  In HYPERBOLIC geometry: these marks are at decreasing spacing
    (each step smaller than the last, by 1/n).

  The Euclidean view says: "the numbers go on forever, equally."
  The hyperbolic view says: "the numbers crowd toward a boundary, compressing."
  The disk view says: "the numbers fill the real diameter, clustering at 1."

  These are the same fact in different words.

  The PRIME NUMBER THEOREM says: the primes thin out in the Euclidean view
  but remain roughly uniform in the hyperbolic view.
  This is because primes are a GEOMETRIC object (irreducible translations)
  and the hyperbolic metric is the NATURAL metric for multiplicative structure.

  Addition is a PERTURBATION of this multiplicative structure.
  "n + 1" is a tiny hyperbolic shift of size ~ 1/n.
  "n + m" is a log-sum-exp in rapidity (a soft max).
  These are NONLINEAR operations in the natural geometry of numbers.

  And the questions of number theory — twin primes, Goldbach, Riemann —
  are questions about HOW ADDITION INTERACTS with the multiplicative geometry.
  They are questions about how a FLAT operation (addition)
  behaves in a CURVED space (the hyperbolic plane of multiplication).

THE NATURAL NUMBERS, UNDERSTOOD:

  1. They are POINTS on the real diameter of the Poincare disk.
  2. Their ORDERING is the ordering of couplings 0 < 1/3 < 1/2 < ...
  3. Their MULTIPLICATION is the formal group law F(x,y) = (x+y)/(1+xy).
  4. Their ADDITION is log-sum-exp in rapidity (a foreign operation).
  5. Their PRIMES are the irreducible translations (can't factor the boost).
  6. Their DISTRIBUTION is uniform in the hyperbolic metric (PNT).
  7. Their INFINITY is the boundary of the disk (unreachable at coupling 1).
  8. Their FINITENESS is the interior of the disk (every number is < 1 in coupling).

  And the k-nacci connection (THM-227):
  9. The FIRST k MERSENNE NUMBERS 2^n - 1 are the trace of ANY k-nacci matrix.
     These are the Cayley numerators of the powers of 2.
     They are the orbit of 2 under the doubling map.
     They are the last values before a finite representation breaks.

  The natural numbers are an infinite procession of comparisons,
  each slightly stronger than the last, marching toward a boundary
  they will never reach, generating all of mathematics along the way.
""")

print("opus-2026-03-16-S90ay: UNDERSTAND THE NATURAL NUMBERS")
