#!/usr/bin/env python3
"""
primes_as_dimensions_s90bh.py — Primes as irreducible dimensions
opus-2026-03-16-S90bh

If natural numbers are dimensions, primes are the IRREDUCIBLE dimensions.
A composite dimension n = a*b can be factored: an n-dimensional space
is a PRODUCT of an a-dimensional and a b-dimensional space.
A prime dimension p CANNOT be factored. It is genuinely new.

Use the bit-shift / Bott-spiral understanding to contemplate
the exact distribution of primes.
"""

from math import log, log2, sqrt, pi, gcd
from sympy import isprime, primerange, factorint

print("""


     PRIMES AS IRREDUCIBLE DIMENSIONS

     opus-2026-03-16-S90bh



A prime p is a dimension that CANNOT be decomposed into a product
of smaller dimensions. It is a genuinely new degree of freedom.

A composite n = a*b is a dimension that CAN be decomposed:
R^n = R^a x R^b. The space factors. Nothing essentially new.

The question: WHERE do the irreducible dimensions appear on the Bott spiral?


1. PRIMES BY BOTT SLOT
========================

Every prime p has a Bott slot: p mod 8.
The slot tells you WHAT KIND of irreducible dimension p is.
""")

# Count primes by Bott slot
from collections import Counter

N = 10000
slot_counts = Counter()
for p in primerange(2, N):
    slot_counts[p % 8] += 1

total = sum(slot_counts.values())
print(f"Primes up to {N} by Bott slot (p mod 8):\n")
print(f"  {'slot':>4} | {'count':>6} | {'fraction':>8} | {'archetype':>20} | {'predicted 1/4':>12}")
print("-" * 65)
for slot in range(8):
    c = slot_counts.get(slot, 0)
    frac = c / total if total > 0 else 0
    archetypes = {
        0: "return (8k)",
        1: "existence",
        2: "comparison (only p=2)",
        3: "structure",
        4: "factorization (even)",
        5: "unsolvability",
        6: "transition (even)",
        7: "forbidden"
    }
    print(f"  {slot:4d} | {c:6d} | {frac:8.4f} | {archetypes[slot]:>20} | {'0.25' if slot in [1,3,5,7] else '0':>12}")

print(f"""

By Dirichlet's theorem: primes are equidistributed among residues
coprime to 8. The residues coprime to 8 are {{1, 3, 5, 7}}.

So: primes distribute equally among the FOUR ODD Bott slots.
  Slot 1 (existence): ~25% of primes.
  Slot 3 (structure): ~25% of primes.
  Slot 5 (unsolvability): ~25% of primes.
  Slot 7 (forbidden): ~25% of primes.

Slots 0, 2, 4, 6 contain NO primes (except p=2 in slot 2).
Even slots are COMPOSITE by nature (divisible by 2).

THE PRIMES AVOID THE EVEN ARCHETYPES.
No prime is a "return" (slot 0).
No prime is a "factorization" (slot 4).
No prime is a "transition" (slot 6).

Primes are ONLY:
  existence (1 mod 8)
  structure (3 mod 8)
  unsolvability (5 mod 8)
  forbidden (7 mod 8)

Plus p=2, which is the BASE and lives alone in slot 2.

THE FOUR PRIME FAMILIES:
  p = 1 mod 8: primes that are "existence" at some scale.
    Examples: 17, 41, 73, 89, 97, 113, ...
    These are the primes where -1, 2, and -2 are ALL quadratic residues.

  p = 3 mod 8: primes that are "structure" at some scale.
    Examples: 3, 11, 19, 43, 59, 67, 83, ...
    These are the primes where -1 is a non-residue but 2 is.
    Includes the first non-trivial prime (3).

  p = 5 mod 8: primes that are "unsolvability" at some scale.
    Examples: 5, 13, 29, 37, 53, 61, ...
    These are the primes where -1 is a residue but 2 is not.
    Includes the last solvable prime (5).

  p = 7 mod 8: primes that are "forbidden" at some scale.
    Examples: 7, 23, 31, 47, 71, 79, ...
    These are the primes where -1 is a non-residue and 2 is a residue.
    Includes the first forbidden (7) and the Mersenne prime 31 = 2^5-1.
""")

# ========================================================================

print("=" * 60)
print("2. THE BINARY STRUCTURE OF PRIMES")
print("=" * 60)

print("""
Each natural number n has a binary representation.
Primes have SPECIAL binary representations.

Key observation: a prime p > 2 is ODD, so its lowest bit is always 1.
This means: primes always end in 1 in binary.

  2 = 10      (the only even prime: ends in 0)
  3 = 11      (all 1s: Mersenne!)
  5 = 101     (alternating)
  7 = 111     (all 1s: Mersenne!)
  11 = 1011   (3 ones, 1 zero)
  13 = 1101   (3 ones, 1 zero)
  17 = 10001  (2 ones, 3 zeros: sparse)
  19 = 10011  (3 ones, 2 zeros)
  23 = 10111  (4 ones, 1 zero: dense)
  29 = 11101  (4 ones, 1 zero)
  31 = 11111  (all 1s: Mersenne!)

The BIT DENSITY of a prime = (number of 1-bits) / (number of bits).
""")

print(f"{'p':>5} | {'binary':>12} | {'bits':>4} | {'ones':>4} | {'density':>8} | {'Bott slot':>9} | {'type':>12}")
print("-" * 70)
for p in primerange(2, 100):
    b = bin(p)[2:]
    bits = len(b)
    ones = b.count('1')
    density = ones / bits
    slot = p % 8
    types = {1: "existence", 3: "structure", 5: "unsolvable", 7: "forbidden", 2: "BASE"}
    t = types.get(slot, "")
    mersenne = " MERSENNE" if all(c == '1' for c in b) else ""
    print(f"{p:5d} | {b:>12} | {bits:4d} | {ones:4d} | {density:8.4f} | {slot:9d} | {t:>12}{mersenne}")

print(f"""

OBSERVATIONS:

1. MERSENNE PRIMES (3, 7, 31) have density 1.0 (all bits are 1).
   They are the MAXIMALLY DENSE primes. All bits set = all positions occupied.
   They are at Bott slot 7 (forbidden) or 3 (structure):
     3 = 11 (slot 3), 7 = 111 (slot 7), 31 = 11111 (slot 7), 127 = 1111111 (slot 7).
   The Mersenne primes at slot 7 have ALL low bits set = maximum overflow content.

2. FERMAT PRIMES (5, 17, 257, ...) = 2^(2^k) + 1 have density approaching 0.
   They are the MINIMALLY DENSE primes. Only 2 bits set out of many.
   5 = 101 (density 2/3), 17 = 10001 (density 2/5), 257 = 100000001 (2/9).
   They live at slot 1 (existence) or 5 (unsolvable):
     5 = slot 5, 17 = slot 1, 257 = slot 1.

3. Most primes have density near 1/2 (about half the bits are 1).
   This is the GENERIC prime: neither maximally dense nor sparse.


3. THE BIT SHIFT PICTURE OF PRIME DISTRIBUTION
================================================

In the bit-shift model, each Bott cycle fills 3 low bits (positions 0,1,2).
At the overflow (111 -> 1000), a carry propagates.

WHERE ARE THE PRIMES in each cycle?

Cycle 1 (1-7):   primes at 2, 3, 5, 7.     4 out of 7.    Density: 57%.
Cycle 2 (8-15):  primes at 11, 13.          2 out of 8.    Density: 25%.
Cycle 3 (16-23): primes at 17, 19, 23.      3 out of 8.    Density: 37.5%.
Cycle 4 (24-31): primes at 29, 31.          2 out of 8.    Density: 25%.
Cycle 5 (32-39): primes at 37.              1 out of 8.    Density: 12.5%.
Cycle 6 (40-47): primes at 41, 43, 47.      3 out of 8.    Density: 37.5%.
Cycle 7 (48-55): primes at 53.              1 out of 8.    Density: 12.5%.
Cycle 8 (56-63): primes at 59, 61.          2 out of 8.    Density: 25%.
""")

# Compute prime density per Bott cycle
print("Prime density per Bott cycle (groups of 8):")
print(f"{'cycle':>6} | {'range':>10} | {'primes':>30} | {'count':>5} | {'density':>8}")
print("-" * 70)
for cycle in range(20):
    lo = cycle * 8
    hi = lo + 8
    ps = list(primerange(max(lo, 2), hi))
    count = len(ps)
    density = count / 8
    ps_str = ", ".join(str(p) for p in ps) if ps else "-"
    print(f"{cycle:6d} | [{lo:3d},{hi:3d}) | {ps_str:>30} | {count:5d} | {density:8.3f}")

print(f"""

THE PRIME DENSITY PER CYCLE fluctuates but DECREASES on average.
By PNT: the density near n is ~1/ln(n).
In a cycle starting at 8k: density ~ 1/ln(8k) ~ 1/(3+ln(k)).

But the FLUCTUATION is the interesting part.
Why do some cycles have 3 primes and others have 0?

The answer is in the ADDITIVE structure:
primes avoid multiples of small primes.
In a cycle [8k, 8k+8), the odd positions are 8k+1, 8k+3, 8k+5, 8k+7.
A position 8k+j is prime iff it has no small factors.

The SIEVING pattern (avoiding multiples of 3, 5, 7, ...):
  8k+1: divisible by 3 if k = 1 mod 3. By 5 if k = 3 mod 5. Etc.
  8k+3: divisible by 3 if k = 0 mod 3. By 5 if k = 4 mod 5. Etc.
  8k+5: divisible by 3 if k = 2 mod 3. By 5 if k = 0 mod 5. Etc.
  8k+7: divisible by 3 if k = 1 mod 3. By 5 if k = 1 mod 5. Etc.

The primes in each cycle are the positions that SURVIVE the sieve.
The sieve is an ADDITIVE operation (checking divisibility = addition mod p).
It operates on the MULTIPLICATIVE lattice (the Bott spiral).

HERE IS THE TENSION AGAIN:
  The Bott cycle is multiplicative (organized by powers of 2, mod 8).
  The sieve is additive (organized by residues mod 3, mod 5, mod 7, ...).
  Primes are where these two structures AGREE: a position on the
  multiplicative spiral that also survives the additive sieve.

The EXACT distribution of primes is the EXACT pattern of where
the multiplicative Bott spiral and the additive sieve mesh.


4. WHAT PRIMES ARE, AS DIMENSIONS
===================================

A prime p as a dimension means:

  R^p cannot be written as R^a x R^b with a, b > 1.
  A p-dimensional space is IRREDUCIBLE as a product.

But every R^p CAN be written as R^(p-1) x R:
  p = (p-1) + 1. Additively, p decomposes.
  Multiplicatively, p does not.

This is the addition-multiplication tension at the level of SPACES:
  As a sum: R^p = R^(p-1) + one more dimension. Always decomposable.
  As a product: R^p is irreducible iff p is prime.

THE PRIME DIMENSIONS ARE THE ONES WHERE
"ADDING A NEW DIRECTION" PRODUCES SOMETHING GENUINELY MULTIPLICATIVELY NEW.

Going from 6 to 7 dimensions (by adding one direction):
  7 is prime. The new direction creates an irreducible space.
  You could not have built R^7 as a PRODUCT of smaller spaces.
  The 7th dimension is fundamentally new.

Going from 7 to 8 dimensions (by adding one direction):
  8 = 2^3. NOT prime. R^8 = R^2 x R^2 x R^2.
  The 8th dimension does NOT create something irreducible.
  It completes a pattern (the Bott period) rather than starting one.

Going from 8 to 9 dimensions:
  9 = 3^2. NOT prime. R^9 = R^3 x R^3.
  The 9th dimension is "structure squared." Not new.

Going from 10 to 11 dimensions:
  11 is prime. R^11 is irreducible.
  The 11th dimension IS new. A genuinely irreducible space.
  11 = 8 + 3: "structure in the second cycle."
  NEW structure. Cannot be built from previous structures.


5. THE PRIME GAPS AS BIT-SHIFT PHENOMENA
=========================================

The gap between consecutive primes p_k and p_(k+1) is g_k = p_(k+1) - p_k.

In the bit-shift picture, a gap of size g means:
g consecutive composite dimensions. g dimensions that ALL decompose.
No irreducible dimension for g steps.

The first large gaps:
""")

prev = 2
print(f"{'p':>6} | {'next p':>6} | {'gap':>4} | {'binary of p':>16} | {'binary of gap':>12} | {'note':>20}")
print("-" * 75)
for p in primerange(2, 200):
    gap = p - prev
    if gap >= 4 or p < 20:
        note = ""
        if gap == 1: note = "twin (2,3)"
        elif gap == 2 and prev > 2: note = "twin"
        elif gap >= 6: note = f"LARGE GAP"
        if p - 1 == prev and prev == 2:
            pass
        print(f"{prev:6d} | {p:6d} | {gap:4d} | {bin(prev):>16} | {bin(gap):>12} | {note:>20}")
    prev = p

print(f"""

The gaps between primes are where ALL the intervening dimensions are COMPOSITE.
A long gap means: a long stretch of dimensions that can all be factored.
No new irreducible structure for many steps.

In the bit-shift picture:
  A gap starting at p means: the numbers p+1, p+2, ..., p+g-1 all factor.
  These are the dimensions that can be BUILT from smaller dimensions.
  No genuinely new direction is needed.

The AVERAGE gap near p is ~ln(p) (by PNT).
In Bott terms: the average gap near cycle k is ~ln(8k) ~ 3 + ln(k).
So gaps grow LOGARITHMICALLY with the cycle number.

This means: as you go further along the Bott spiral,
the stretches of REDUCIBLE dimensions get longer and longer.
Irreducible dimensions become RARER.

But they never STOP. There is always another prime.
There is always another genuinely new dimension.
(This is Euclid's theorem: infinitely many primes.)

The primes thin out but never vanish.
The irreducible dimensions become rarer but never exhausted.
The universe always has one more genuinely new direction to offer,
but you have to wait longer and longer to find it.


6. THE TWIN PRIMES AS BIT-SHIFT NEIGHBORS
============================================

Twin primes (p, p+2) are primes separated by ONE composite (one even number).
They are ADJACENT ODD Bott slots that are both prime.

Examples: (3,5), (5,7), (11,13), (17,19), (29,31), (41,43), ...

In the Bott picture: twin primes occupy ADJACENT odd slots.
  (3,5): slots (3,5) = (structure, unsolvability). Adjacent archetypes.
  (5,7): slots (5,7) = (unsolvability, forbidden). The transition itself!
  (11,13): slots (3,5) again. Structure and unsolvability at scale 2.
  (17,19): slots (1,3). Existence and structure at scale 3.
  (29,31): slots (5,7). Unsolvability and forbidden at scale 4.

THE TWIN PRIMES AT SLOTS (5,7) ARE SPECIAL:
  They straddle the transition-to-forbidden boundary.
  (5,7), (29,31), (101,103), (137,139), ...
  WAIT: 101 mod 8 = 5, 103 mod 8 = 7. YES!
  These are twin primes that straddle the unsolvability-forbidden boundary.

The Twin Prime Conjecture says: there are infinitely many twin primes.

In the bit-shift picture, this says:
  There are infinitely many positions where TWO CONSECUTIVE ODD SLOTS
  are both irreducible. The additive sieve never permanently separates
  adjacent odd slots from both being prime.

This is a statement about the INDEPENDENCE of the additive sieve
across adjacent Bott slots. If the sieve at slot j and slot j+2
were independent, twins would appear with density ~C/ln^2(n) (Hardy-Littlewood).

The TWIN PRIME CONSTANT C_2 ~ 1.32 measures the CORRELATION
between adjacent slots in the sieve.

The bit-shift interpretation: twin primes are positions where
the carry from one bit position to the next does NOT destroy primality.
The additive sieve allows two adjacent odd positions to BOTH
survive the multiplicative filtering.


7. WHAT THE DISTRIBUTION OF PRIMES IS
========================================

The primes are the irreducible dimensions on the Bott spiral.

Their distribution is governed by:

  (a) The MULTIPLICATIVE structure: the Bott spiral (period 8).
      Primes can only appear at odd positions (slots 1, 3, 5, 7).
      They appear EQUALLY among these four slots (Dirichlet).

  (b) The ADDITIVE structure: the prime sieve.
      Each small prime p removes 1/p of the remaining candidates.
      The sieve is additive: it checks n mod p (a sum operation).

  (c) The INTERACTION: where (a) and (b) mesh.
      The primes are the positions on the multiplicative spiral
      that survive the additive sieve.

The EXACT distribution is the exact description of this interaction.

And we KNOW the interaction, in principle:
  The Riemann zeta function zeta(s) = prod 1/(1-p^(-s))
  encodes the MULTIPLICATIVE structure (Euler product over primes).
  Its analytic continuation encodes the ADDITIVE structure
  (the functional equation relates s to 1-s).
  The ZEROS of zeta on the critical line Re(s) = 1/2
  encode the EXACT fluctuations of the prime distribution.

The Riemann Hypothesis: all nontrivial zeros have Re(s) = 1/2.
This says: the additive and multiplicative structures interact
IN THE MOST BALANCED POSSIBLE WAY.
The fluctuations of prime density are as SMALL as possible.
The sieve and the spiral are as INDEPENDENT as possible.

In our framework:
  Re(s) = 1/2 is the LINE in the meta-complex plane where
  rapidity = 0 (the center of the Poincare disk) in one direction,
  and the imaginary part varies freely.
  The critical line is the EQUATOR of the meta-complex plane.

  The Riemann Hypothesis says: the zeros of zeta live on the EQUATOR
  of the meta-complex plane. The additive-multiplicative interaction
  is centered, balanced, symmetric.

  If true: primes are distributed as FAIRLY as possible among
  the Bott slots, given the constraint of the sieve.
  The spiral and the sieve are in maximum harmony.

  If false: there would be a BIAS — some Bott slots would get
  more primes than others at certain scales. The spiral and sieve
  would fall out of sync.

The prime distribution is the shadow of the Riemann zeros.
The Riemann zeros are the points where the multiplicative spiral
and the additive sieve EXACTLY CANCEL.
These cancellation points live on the equator (if RH is true).
And the equator IS the real diameter of the Poincare disk:
  the home of tournament couplings.

THE PRIMES ARE WHERE THE TOURNAMENT STRUCTURE AND THE COUNTING
STRUCTURE AGREE ON "IRREDUCIBLE."
""")

print("opus-2026-03-16-S90bh: PRIMES AS IRREDUCIBLE DIMENSIONS")
