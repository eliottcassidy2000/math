#!/usr/bin/env python3
"""
even_odd_s90am.py — Even and Odd
opus-2026-03-16-S90am
"""

from math import log, sqrt, pi

print("""


                           Even and Odd



H is always odd.

Why?

H = 1 + 2α₁ + 4α₂ + 8α₃ + ...
  = 1 + (something even)
  = odd.

The 1 comes from α₀ = 1: the empty set. Something exists.
Everything else — α₁, α₂, α₃ — is multiplied by a power of 2.
Powers of 2 are even (for k ≥ 1).
Even times anything is even.
So everything added after existence is even.

The ODDNESS of H comes entirely from the 0th axiom.
Existence is odd. Correction is even.
The universe is 1 + (adjustments).

                              *

The number line splits:

  ODD:   1,  3,  5,  7,  9, 11, 13, ...
  EVEN:  2,  4,  6,  8, 10, 12, 14, ...

In our theory:
  ODD  = the tournament numbers. H(T) ∈ {1, 3, 5, 9, 11, ...}.
  EVEN = the adjustments. 2α₁, 4α₂, 8α₃, ...

The forbidden values 7 and 21 are ODD.
They are forbidden not because they're even (they're not)
but because they require an EVEN ADJUSTMENT of exactly 6 or 20,
and no tournament can produce that specific even number.

  H = 7 requires: even part = 6 = 2·3.  But α₁=3 with α₂=0: impossible.
  H = 21 requires: even part = 20. Many decompositions, all impossible.

The forbidden values are odd numbers whose EVEN PART is unbuildable.

                              *

On the Cayley line:

  ODD n has address (n-1)/2 / (n+1)/2 where numerator + denominator = n.
  EVEN n has address (n-1)/(n+1) which does NOT simplify (both odd).

  Odd numbers DECOMPOSE on the Cayley line: a + b = n.
  Even numbers DO NOT decompose: (n-1) and (n+1) are coprime odd numbers.

  The Stern-Brocot tree: odd numbers are CHILDREN (decomposable).
  Even numbers are MEDIANTS (indecomposable at this level).
  The odd/even partition IS the parent/child partition of Stern-Brocot.

                              *

In the cosh/sinh decomposition:

  Q^m = cosh(2m·arctanh(x)) + sinh(2m·arctanh(x))

  cosh = the EVEN part (symmetric in x: cosh(-x) = cosh(x)).
  sinh = the ODD part (antisymmetric: sinh(-x) = -sinh(x)).

  The even part is UNDIRECTED (doesn't care about orientation).
  The odd part is DIRECTED (reverses under complement T^op).

  The tournament H(T) is the SUM of even and odd:
    H = cosh + sinh = total (both directions count).

  The SIGNED tournament count F(T, -1) is the DIFFERENCE:
    F(-1) = cosh - sinh = net (direction matters).

  cosh + sinh = e^x  (the exponential: GROWTH).
  cosh - sinh = e^{-x} (the inverse: DECAY).

  The even and odd parts combine into growth and decay.
  Growth = having more timelines than anti-timelines.
  Decay = having fewer.

  H > 0 always (Rédei: growth always wins for the total count).
  F(-1) can be positive or negative (net direction can go either way).

                              *

Even and odd in the eigenvalues:

  The transfer matrix M(x) has char poly λ³ - λ² - xλ - x.

  This factors as: λ²(λ-1) - x(λ+1) = 0.
  Or: (λ-1)(λ²) = x(λ+1).

  The LEFT side: (λ-1) · λ². One ODD factor (λ-1) times one EVEN power (λ²).
  The RIGHT side: x · (λ+1). The coupling times one ODD-shifted factor.

  At the eigenvalues:
    (λ-1)·λ² = x·(λ+1)
    (λ-1)/(λ+1) = x/λ²

  The LEFT side: (λ-1)/(λ+1) = 1 - 2/(λ+1).
  This is 1 minus an EVEN correction. ODD base, EVEN adjustment.
  Same structure as H = 1 + even.

                              *

Even and odd in the discriminant:

  Δ(x) = 4x(x² - 11x - 1).

  4 = even (the coefficient).
  x = the coupling (can be even or odd at integer values).
  x² - 11x - 1 = a quadratic. At odd x: x² is odd, 11x is odd, 1 is odd.
    odd - odd - odd = odd. So (x²-11x-1) is odd at odd x.
  At even x: x² is even, 11x is even, 1 is odd.
    even - even - odd = odd. So (x²-11x-1) is ALWAYS ODD at integers!

  Therefore: Δ(x) = 4 · x · (odd) = 4x · (odd).
  At odd x: Δ = 4 · odd · odd = 4 · odd. Always divisible by 4, never by 8.
  At even x: Δ = 4 · even · odd = 8 · (x/2) · odd.

  Δ(1) = 4 · 1 · (-11) = -44 = -4 · 11.    (4 × odd)
  Δ(2) = 4 · 2 · (-19) = -152 = -8 · 19.    (8 × odd)
  Δ(3) = 4 · 3 · (-25) = -300 = -4 · 75.    (4 × odd, since 3 odd)

  The discriminant at ODD coupling: divisible by 4 but not 8.
  The discriminant at EVEN coupling: divisible by 8.

  This even/odd split in the discriminant controls the
  PARITY OF THE EIGENVALUE STRUCTURE.

                              *

Even and odd in the Bott period:

  Tr(M^n) mod 8 has period 8.
  The period IS 8 = 2³ (a power of 2 = an even number).

  The trace sequence mod 2: [1, 1, 1, 1, ...] (always 1 = always odd).
  The trace sequence mod 4: period 4 = [3, 1, 3, 3, ...].
  The trace sequence mod 8: period 8 = [3, 1, 3, 7, 3, 5, 7, 7].

  The ODD modular behavior (mod 2): TRIVIAL (always 1).
  The EVEN modular behavior (mod 4, 8, ...): NONTRIVIAL (periodic).

  The interesting structure is in the EVEN part of the modular arithmetic.
  The odd part just says "exists" (= 1 mod 2 = Rédei).
  The even part says "how" (= the specific remainder mod 4, 8, ...).

                              *

The deepest even/odd:

  2 = the binary base. The number of outcomes. EVEN.
  1 = existence. The 0th axiom. ODD.

  Every tournament is:
    1 (odd: existence) + C(n,2) arcs (each contributing 1 bit).
    C(n,2) = n(n-1)/2.
    For odd n: C(n,2) = odd·even/2 = odd·(even/2). Could be odd or even.
    For even n: C(n,2) = even·odd/2 = (even/2)·odd. Always odd? No.
      C(4,2) = 6 (even). C(6,2) = 15 (odd). C(8,2) = 28 (even).

  The NUMBER OF ARCS alternates between even and odd:
    n: 3  4  5  6  7  8  9 10 11
    C: 3  6 10 15 21 28 36 45 55
    p: O  E  E  O  O  E  E  O  O
       ^     ^        ^
    The parity of C(n,2) follows a pattern with period 4:
      O, E, E, O, O, E, E, O, ...
    (determined by n mod 4).

  When C(n,2) is even: the tournament space has 2^(even) = (2^2)^k elements.
    The number of tournaments is a perfect power of 4.
  When C(n,2) is odd: the tournament space has 2^(odd) elements.
    The number is 2 times a perfect power of 4.

  This even/odd split in the dimension of the tournament space
  affects EVERYTHING: which H values are achievable, which β₁ patterns
  appear, which score sequences are possible.

                              *

The ONE SENTENCE:

  ODD is what exists.
  EVEN is what changes.

  H(T) is odd because existence (α₀ = 1) is odd,
  and change (2ᵏαₖ for k ≥ 1) is always even.

  The 0th axiom (something is) contributes 1 = odd.
  The 1st axiom (things differ) contributes 2α₁ = even.
  Every comparison adds an even amount to an odd base.

  The universe is odd. Its dynamics are even.
  The total is always odd.

  And the forbidden values 7 and 21 are odd numbers
  whose even adjustments (6 and 20) are SELF-CONTRADICTORY:
  they cannot be assembled from the even building blocks
  (2ᵏαₖ) available to tournament conflict graphs.

  7 = 1 + 6, but 6 cannot be built as 2α₁ + 4α₂ + ... with valid αₖ.
  21 = 1 + 20, but 20 cannot be built with valid αₖ.

  The universe refuses to adjust by exactly 6 or 20.
  These adjustments require circuits that contradict themselves.
  The even part fails. The odd base persists.


""")

print("opus-2026-03-16-S90am: EVEN AND ODD")
