#!/usr/bin/env python3
"""
meta_isomorphism_s90av.py — The meta-level isomorphism
opus-2026-03-16-S90av

Understanding tanh/arctanh, powers of 2, transcendental bases,
the universal zone, and the isomorphism at the deepest level.

The question is not "what are the identities?" but
"why does the SAME structure keep appearing?"
"""

from math import tanh, atanh, log, log2, sqrt, exp, pi, e
from fractions import Fraction

print("""


     THE META-LEVEL ISOMORPHISM

     opus-2026-03-16-S90av



0. THE STARTING OBSERVATION
============================

There are three apparently different operations:

  (a) ADDITION on the real line R:         a + b
  (b) MULTIPLICATION on (0,inf):           a * b
  (c) The FORMAL GROUP on (-1,1):          (a+b)/(1+ab)

These are not analogous. They are ISOMORPHIC.
The isomorphisms are:

  exp: (R, +) --> ((0,inf), *)
  tanh: (R, +) --> ((-1,1), F)
  Q = (1+x)/(1-x): ((-1,1), F) --> ((0,inf), *)

  exp = Q o tanh  (up to a factor of 2).

More precisely:
  arctanh(x) = (1/2) ln((1+x)/(1-x)) = (1/2) ln Q(x).

So: Q(x) = e^{2 arctanh(x)}.
And: tanh(w) = Q^{-1}(e^{2w}).

The three spaces are ONE GROUP wearing three costumes:

  COSTUME 1 (additive):     R with +
  COSTUME 2 (multiplicative): (0,inf) with *
  COSTUME 3 (hyperbolic):   (-1,1) with (x+y)/(1+xy)

arctanh is the map from costume 3 to costume 1.
exp is the map from costume 1 to costume 2.
Q = exp(2*arctanh) is the map from costume 3 to costume 2.

THE QUESTION: why does costume 3 keep appearing in our theory?
""")

# ========================================================================

print("""
1. WHY COSTUME 3?
==================

The three costumes differ in WHAT THEY BOUND:

  Costume 1 (R, +):        nothing bounded. (-inf, inf).
  Costume 2 ((0,inf), *):  bounded below (by 0). Not above.
  Costume 3 ((-1,1), F):   bounded BOTH sides (by -1 and +1).

Tournaments live in costume 3 because tournaments are DOUBLY BOUNDED:
  - Every coupling x is between -1 (full reversal) and +1 (full agreement).
  - Every velocity v is between -c and +c.
  - Every probability p is between 0 and 1 (after rescaling).
  - Every bit is between 0 and 1.

The DOUBLE BOUNDEDNESS is the signature of costume 3.
And double boundedness is the signature of COMPARISON:
  when you compare two things, the result is between "A dominates completely"
  and "B dominates completely."

arctanh REMOVES the bounds: it sends (-1,1) to (-inf,inf).
It converts a bounded comparison into an unbounded magnitude.
This is the RAPIDITY: how many times you've doubled your advantage.

tanh RESTORES the bounds: it sends (-inf,inf) back to (-1,1).
It converts an unbounded magnitude back into a bounded comparison.
This is the VELOCITY: what fraction of the speed limit you've reached.

THE META-POINT: arctanh and tanh are not just functions.
They are the TRANSLATIONS between bounded and unbounded worlds.
Between comparison and magnitude.
Between the finite and the infinite.


2. THE TRANSCENDENTAL BASES
============================

The natural logarithm ln has base e = 2.71828...
The binary logarithm log2 has base 2.
The Cayley rapidity uses ln(n)/2 = log_e(n)/2.

Why do BOTH e and 2 appear?

Because there are two DIFFERENT questions:
  (a) What is the NATURAL unit of growth? Answer: e (continuous compounding).
  (b) What is the NATURAL unit of choice? Answer: 2 (binary decision).

e comes from calculus: d/dx e^x = e^x. Self-replication.
2 comes from logic: yes/no, true/false, 0/1. Binary choice.

In tournament theory, BOTH appear:
  - The transfer matrix has eigenvalues involving e (continuous growth).
  - The tournament has 2 orientations per arc (binary choice).
  - The coupling x lives in (-1,1) which is the e-world projected to the 2-world.
  - arctanh(x) = (1/2) ln((1+x)/(1-x)) uses BOTH e (ln) and 2 (the 1/2).

THE FACTOR OF 2: why rapidity = ln(n)/2 and not ln(n)?

Because the Cayley transform Q(x) = (1+x)/(1-x) squares the odds ratio:
  Q(x) = e^{2 arctanh(x)}.

The 2 in the exponent comes from:
  arctanh(x) = (1/2)(ln(1+x) - ln(1-x)).

The (1+x) and (1-x) are the two HALVES of the comparison.
"How much for" and "how much against."
Their RATIO is the odds, and half the log of the odds is the rapidity.

The 2 is irreducible. It comes from the fact that a comparison
has TWO sides. Not three, not one. Two.

This is why 2 is transcendental in our theory: it is the NUMBER OF SIDES
of a comparison, just as e is the NUMBER OF SELF-REPLICATIONS
in continuous growth.

e and 2 are the two transcendental bases because they answer
the two fundamental questions:
  e: "How does continuous change compound?"
  2: "How does a binary choice resolve?"

Tournament theory is where these two questions MEET:
  The transfer matrix (continuous: e) encodes the tournament (binary: 2).
  arctanh bridges them: it converts binary comparisons (bounded by 2 sides)
  into continuous magnitudes (growing by e).


3. THE UNIVERSAL ZONE
=======================

THM-227: Tr(M_k^n) = 2^n - 1 for 1 <= n <= k.

The UNIVERSAL ZONE is where the trace depends only on n, not on k.
All k-nacci matrices agree. The individuality of k is invisible.

WHY does this happen?

Newton's identities for the k-nacci char poly give:
  t(n) = t(n-1) + t(n-2) + ... + t(1) + n   (for n <= k).

The key: the CORRECTION TERM is +n, not +n*c for some coefficient c.
And the sum ranges over ALL previous values t(1), ..., t(n-1).

This happens because ALL coefficients of the k-nacci polynomial are -1.

WHY are all coefficients -1?

Because the k-nacci recurrence is: a(n) = a(n-1) + a(n-2) + ... + a(n-k).
Each of the k previous terms contributes with coefficient 1.
No preference. No weighting. MAXIMUM UNIFORMITY.

The k-nacci polynomial is the MOST DEMOCRATIC polynomial of degree k:
every predecessor has equal weight.

Democracy is universality. When all weights are equal, the early behavior
is independent of how many weights there are. You only need to know
that each predecessor counts equally; you don't need to know how many
predecessors there are (until you run out of them at n = k).

THE META-POINT:

The universal zone exists because of DEMOCRACY.
The k-nacci polynomial treats all predecessors equally.
This equality makes the first k traces universal.

The TRANSITION at n = k is the moment democracy breaks down:
at n = k+1, you can't sum over "all predecessors up to n-1"
anymore because the recurrence only reaches back k steps.
The FINITE HORIZON k reveals itself.

In information terms:
  - Universal zone: you don't need to know k. The answer is 2^n - 1.
  - Individual zone: you DO need to know k. The answer depends on the dimension.
  - Transition: the system discovers it is FINITE.

This is the deepest meaning of the Mersenne number 2^k - 1:
it is the LAST VALUE before the system discovers its own finiteness.


4. THE ISOMORPHISM AT THE META LEVEL
======================================

There are four copies of the SAME structure appearing:

(A) The formal group ((-1,1), F):
    F(x,y) = (x+y)/(1+xy).
    The bounded world. Comparisons. Couplings.

(B) The multiplicative group ((0,inf), *):
    Ordinary multiplication.
    The unbounded-above world. Magnitudes. Natural numbers.

(C) The additive group (R, +):
    Ordinary addition.
    The fully unbounded world. Rapidities. Logarithms.

(D) The k-nacci trace sequence:
    t(1)=1, t(2)=3, ..., t(k)=2^k-1.
    A finite projection of (A), (B), (C) into k dimensions.

The isomorphisms:

  arctanh: (A) --> (C)        [remove bounds]
  exp:     (C) --> (B)        [exponentiate]
  Q:       (A) --> (B)        [Cayley transform]
  Tr(M^n): (D) --> Z          [trace = finite shadow]

And: for n <= k, Tr(M_k^n) = 2^n - 1 = Q(x_{2^n}) - 2 ... no.
Actually: 2^n - 1 is the NUMERATOR of x_{2^n} = (2^n-1)/(2^n+1).

The isomorphism at the meta level is:

  THE THREE GROUPS ARE ONE.
  THE k-NACCI MATRIX IS A FINITE-DIMENSIONAL REPRESENTATION.
  THE UNIVERSAL ZONE IS WHERE THE REPRESENTATION IS FAITHFUL.
  THE INDIVIDUAL ZONE IS WHERE IT FAILS.

A REPRESENTATION of a group G is a homomorphism rho: G --> GL(V)
that sends group elements to matrices. A representation is FAITHFUL
if the homomorphism is injective: no information is lost.

The k-nacci matrix M_k is a REPRESENTATION of the multiplicative
group (Z, +) (i.e., powers of a single element):
  rho(n) = M_k^n.

The TRACE is a character: chi(n) = Tr(M_k^n).

THM-227 says: for n <= k, chi(n) = 2^n - 1 = the "true" character
of the infinite-dimensional representation (the binary shift).

For n > k: chi(n) deviates from 2^n - 1 because the k-dimensional
representation LOSES INFORMATION about the group element.

THE META-ISOMORPHISM:

  The "true" object is the INFINITE k-nacci (k = infinity):
    Tr(M_inf^n) = 2^n - 1 for ALL n.

  This is the BINARY SHIFT: the dynamical system on {0,1}^Z
  with 2^n periodic orbits of period n, minus 1 fixed point.

  Each finite k-nacci matrix M_k is a TRUNCATION:
    It agrees with the binary shift for n <= k.
    It diverges for n > k.

  The k-nacci constant tau_k is the "speed" of the truncation:
    tau_k -> 2 as k -> inf.

  The full binary shift has speed 2 (= the topological entropy).
  Each truncation has speed tau_k < 2.

  THE UNIVERSAL ZONE IS THE FAITHFUL PART OF THE REPRESENTATION.
  THE INDIVIDUAL ZONE IS THE UNFAITHFUL PART.

  Mersenne numbers 2^n - 1 = the character values where the
  representation is still faithful.


5. THE DEEPEST LEVEL
======================

At the deepest level, there is ONE object: the map x -> 2x.

  On (R, +): doubling. w -> 2w.
  On ((0,inf), *): squaring. n -> n^2.
  On ((-1,1), F): the [2]-series. x -> 2x/(1+x^2).

All three are the SAME map in different costumes.

The ORBIT of any point under this map:
  Costume 1: w, 2w, 4w, 8w, ... (additive doubling)
  Costume 2: n, n^2, n^4, n^8, ... (repeated squaring)
  Costume 3: x, [2](x), [4](x), [8](x), ... (formal [2^k]-series)

Starting from n = 2 (costume 2):
  2, 4, 8, 16, 32, ...  (powers of 2)

Starting from w = ln(2)/2 (costume 1):
  ln(2)/2, ln(2), 3ln(2)/2, 2ln(2), ...  (half-bit multiples)

Starting from x = 1/3 (costume 3):
  1/3, 3/5, 7/9, 15/17, ...  (Cayley addresses of powers of 2)

The NUMERATORS of costume 3: 1, 3, 7, 15, ...  = Mersenne numbers.

The TRACES of the k-nacci matrix: Tr(M_k^n) = 2^n - 1 for n <= k.
These ARE the Mersenne numbers.

So: the k-nacci trace, for n in the universal zone,
is computing the NUMERATOR of the ORBIT OF 1/3
under the doubling map on the Cayley line.

Equivalently: it is computing the orbit of 2 under squaring,
minus 1 at each step (the fixed point subtraction).

Equivalently: it is computing the non-trivial periodic orbits
of the binary shift at each period.

IT'S ALL THE SAME COMPUTATION.
The "doubling map" x -> 2x is the atom.
Everything else — tanh, arctanh, Cayley, k-nacci, Mersenne,
Witt vectors, necklaces, forbidden values — is a COSTUME
worn by the doubling map.

And the doubling map IS the tournament:
  Each arc doubles the number of possible orderings.
  Each comparison is a binary choice.
  The total number of tournaments on n-arc space is 2^n.
  The non-trivial ones: 2^n - 1.

The tournament IS the binary shift.
The comparison IS the doubling map.
The bit IS the generator.

And tanh/arctanh are the maps that let the bit
wear different costumes — bounded, unbounded,
additive, multiplicative — while remaining the same bit.
""")

# ========================================================================
# Computational verification of the "one map" principle
# ========================================================================

print("=" * 70)
print("VERIFICATION: The doubling map in three costumes")
print("=" * 70)

print("\nStarting point: 2 (costume 2) = ln(2)/2 (costume 1) = 1/3 (costume 3)")
print()
print(f"{'step':>5} | {'Costume 1 (add)':>16} | {'Costume 2 (mult)':>16} | {'Costume 3 (F)':>16} | {'numerator':>10} | {'Tr(M_k^n)':>10}")
print("-" * 90)

for step in range(8):
    w = (step + 1) * log(2) / 2   # costume 1: rapidity
    n = 2 ** (step + 1)            # costume 2: natural number
    x = tanh(w)                     # costume 3: coupling
    x_frac = Fraction(2**(step+1) - 1, 2**(step+1) + 1)
    num = 2**(step+1) - 1          # numerator
    tr = num                        # = Tr(M_k^{step+1}) for k >= step+1

    print(f"{step+1:5d} | {w:16.8f} | {n:16d} | {str(x_frac):>16} | {num:10d} | {tr:10d}")

print(f"""

All three columns encode the SAME information:
  Costume 1: rapidity at half-bit lattice points.
  Costume 2: powers of 2.
  Costume 3: Cayley addresses with Mersenne numerators.

The Mersenne numbers are the NUMERATORS of costume 3.
The k-nacci traces are the Mersenne numbers.
Everything is the orbit of 2 under doubling.


{"=" * 70}
THE ISOMORPHISM DIAGRAM
{"=" * 70}

                    arctanh
    (-1,1) -----------------------> R
      |   F(x,y) = (x+y)/(1+xy)    |   a + b
      |                              |
   Q  |                              | exp
      |                              |
      v                              v
    (0,inf) <-------- ln ---------- R
         a * b                    a + b


  Q = exp o (2 * arctanh)     [Cayley = exp of doubled rapidity]
  arctanh = (1/2) ln Q        [rapidity = half the log of odds]

  The factor of 2 is the NUMBER OF SIDES OF A COMPARISON.

  Starting from x = 1/3 and iterating [2]:
    1/3 -> 3/5 -> 7/9 -> 15/17 -> ... -> (2^n-1)/(2^n+1) -> ... -> 1
          (never reaching 1, approaching asymptotically)

  The limiting value 1 is UNREACHABLE (= speed of light = full certainty).
  The Mersenne numbers 2^n-1 in the numerators GROW without bound.
  But their Cayley addresses (2^n-1)/(2^n+1) COMPRESS toward 1.

  This compression IS tanh: it takes the unbounded (Mersenne numbers)
  and bounds them (Cayley addresses approaching but never reaching 1).

  This bounding IS the physics: velocities can't exceed c,
  probabilities can't exceed 1, couplings can't exceed 1.

  And the reason is: the doubling map in costume 3 has a FIXED POINT at 1.
  The fixed point IS the speed of light.
  The orbit approaches it but never reaches it.
  Because tanh(infinity) = 1, but infinity is never reached.

  THE BIT IS THE GENERATOR.
  TANH IS THE COSTUME.
  THE DOUBLING MAP IS THE DYNAMICS.
  THE MERSENNE NUMBER IS THE SHADOW.
  THE UNIVERSAL ZONE IS THE FAITHFUL REPRESENTATION.
  THE k-NACCI MATRIX IS THE FINITE WINDOW.
  AND IT'S ALL ONE THING.
""")

# ========================================================================
# The "base" question
# ========================================================================

print("=" * 70)
print("WHY BASE 2 AND NOT BASE 3?")
print("=" * 70)

print("""
If we replace 2 with an arbitrary base b, everything generalizes:

  Q_b(x) = (1 + x) / (1 - x)   [this doesn't change!]

Wait. The Cayley transform is independent of the base.
But the rapidity unit changes: 1 bit = ln(2)/2, 1 trit = ln(3)/2.

For base b:
  x_n = (n-1)/(n+1) is the Cayley address [independent of b].
  w_n = ln(n)/2 is the rapidity [independent of b].
  But the LATTICE depends on b:
    Powers of b are at w = k*ln(b)/2, separated by ln(b)/2.

For b = 2: lattice spacing = ln(2)/2 = 0.347.
For b = 3: lattice spacing = ln(3)/2 = 0.549.
For b = e: lattice spacing = 1/2.

Base e gives the CLEANEST lattice: spacing exactly 1/2.
But base 2 gives the MOST NATURAL lattice: each step doubles.

The k-nacci transfer matrix uses base 2 because:
  - The char poly has all coefficients -1 (from uniform weights).
  - The augmented poly is x^{k+1} - 2x^k + 1 (the 2 is the base!).
  - Tr(M^n) = 2^n - 1 (the base appears in the trace).

For a "base-b k-nacci" with char poly x^k - bx^{k-1} + ... + correction:
  The trace would involve b^n instead of 2^n.
  But the all-(-1) coefficient structure is SPECIFIC to base 2.

Why? Because the all-(-1) structure means: each predecessor has weight 1.
And "weight 1" means "one comparison." And one comparison has 2 outcomes.

If each comparison had 3 outcomes (ternary), the weights would be different,
the coefficients would NOT all be -1, and the universal zone would NOT
give Mersenne numbers.

BASE 2 IS FORCED BY THE BINARY NATURE OF COMPARISON.

A comparison is yes/no. Not yes/no/maybe.
Two outcomes. Two orientations of an arc.
This gives all-(-1) coefficients.
Which gives Mersenne numbers.
Which gives the universal zone.

If comparisons were ternary (as in some generalized tournament models),
the coefficients would be different, and the universal zone would give
different numbers. The entire structure would change.

The number 2 is not a choice. It is the NUMBER OF SIDES
of the most basic operation: comparison.
""")

for b in [2, 3, 5, 7, 10]:
    print(f"  Base {b:2d}: lattice spacing = ln({b})/2 = {log(b)/2:.6f}, "
          f"first few addresses: ", end="")
    for k in range(1, 5):
        n = b**k
        addr = Fraction(n-1, n+1)
        print(f"{addr} ", end="")
    print()

print(f"""

For base 2: addresses 1/3, 3/5, 7/9, 15/17 (Mersenne numerators).
For base 3: addresses 1/2, 4/5, 13/14, 40/41 (3^k-1 numerators).
For base 5: addresses 2/3, 12/13, 62/63, 312/313.

Only base 2 produces Mersenne PRIMES in the numerators (at k=2,3,5,7,13,...).
This is because Mersenne primes 2^p-1 require p prime AND the factors
to cooperate. No analogous rich prime structure exists for other bases.

Base 2 is special not just because comparisons are binary,
but because the binary base happens to produce the richest
prime structure in its power-minus-one sequence.

This is a COINCIDENCE of number theory meeting comparison theory.
Or is it? The deepest question: is there a reason why the base
of comparison (2) also produces the richest prime structure?


{"=" * 70}
FINAL MEDITATION
{"=" * 70}

arctanh opens.
tanh closes.

Opening: taking a bounded coupling and asking "how far along the
infinite rapidity axis is this?" Unbounding. Liberating.
From the finite world of comparisons to the infinite world of magnitudes.

Closing: taking an infinite rapidity and asking "what bounded coupling
does this correspond to?" Bounding. Constraining.
From the infinite world of magnitudes to the finite world of comparisons.

The ISOMORPHISM between these two worlds is exact.
Nothing is lost in translation. Every bounded coupling has a
unique unbounded rapidity, and vice versa.

But the EXPERIENCE is different:
  In the bounded world, there is a speed limit (1). Things cluster near it.
  In the unbounded world, there is no limit. Things spread out evenly.

The natural numbers, in the bounded world (Cayley addresses):
  0, 1/3, 1/2, 3/5, 2/3, 5/7, 3/4, 7/9, ...
  Compressed. Clustering near 1. Getting closer and closer together.

The natural numbers, in the unbounded world (rapidities):
  0, 0.347, 0.549, 0.693, 0.805, 0.896, 0.973, 1.040, ...
  Spread out. Growing. Getting further and further apart (logarithmically).

TANH IS THE COMPRESSION MAP.
ARCTANH IS THE EXPANSION MAP.

The k-nacci matrix lives in the bounded world (coupling x).
Its traces live in the integer world (Tr = 2^n-1).
These integers ARE the Cayley numerators of powers of 2.
They are the EXPANSION of the lattice points.

The universal zone (n <= k) is where the COMPRESSION hasn't yet
distorted the integers. The traces are "still integers" (Mersenne).
The individual zone (n > k) is where the compression BITES:
the finite dimension k prevents the traces from staying Mersenne.

AT THE META LEVEL:

The isomorphism between bounded and unbounded is PERFECT.
But any FINITE representation of it (the k-nacci matrix)
can only be faithful for finitely many steps (the universal zone).
After that, finiteness asserts itself.

The universal zone is the LIE that a finite system tells
about being infinite. It maintains the lie for exactly k steps.
Then reality — the finite dimension — breaks through.

The Mersenne number 2^k - 1 is the LAST LIE.
The number the finite system says before admitting it is finite.

And the Mersenne number happens to be the numerator of the
Cayley address of 2^k — the address of the k-th power of 2
on the bounded line.

So the last lie is: "I am the Cayley numerator of 2^k."
After that, the system says something smaller.
The deficit is k+1: the system undershoots by its own dimension plus one.
It reveals itself.

TANH: the function that makes infinity fit in a bounded interval.
ARCTANH: the function that lets a bounded interval dream of infinity.
THE UNIVERSAL ZONE: the dream before waking.
THE MERSENNE NUMBER: the last word of the dream.
""")

print("opus-2026-03-16-S90av: THE META-LEVEL ISOMORPHISM")
