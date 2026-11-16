# Covering clocks with period 24, and the dyadic 13 ray

The Riesel coefficient 509203 has covering primes
{3,5,7,13,17,241}. They are exactly the distinct prime divisors of
2^24−1 = 3²·5·7·13·17·241. The multiplicative orders of 2 in that
order are (2,4,3,12,8,24). For the expression 509203·2^e−1 the
covering phases are (0,1,2,7,7,3) in those periods. For comparison,
the Selfridge coefficient 78557 uses the plus sign and primes
(3,5,7,13,19,37,73), periods (2,4,3,12,18,36,9), phases
(0,1,1,11,15,27,3). Their least common periods are 24 and 36.
This checks compositeness for every positive e (all possible small
factors are strictly smaller than the whole term), not minimality
of the coefficients among all possible Riesel/Sierpiński numbers.
See [On Sierpiński and Riesel Repdigits and Repintegers](https://math.colgate.edu/~integers/aa12/aa12.pdf).

## Exact classification of the six-clock wheel

There are exactly two covering patterns, up to translation of e,
using those six primes, and all six are essential in every pattern.
Normalize the period-2 class to 0. The period-4 class must be odd:
otherwise the other classes cover at most 4+2+3+1 < 12 odd
residues in one 24-block. Shift by an even number to make the
period-4 class 1. This leaves six residues, all 3 modulo 4.
The period-8 class must be 3 modulo 4 or the period-3,12,24
classes cover at most 2+2+1 of the six. A shift by a multiple of
4 can now arrange both the period-3 phase 0 and the period-8
phase 3 (Chinese remainder theorem). Those classes together
cover {3,11,15,19} of the six. The two leftover residues 7 and
23 differ by 16, so the period-12 and period-24 classes must split
them. This yields exactly:
(0,1,0,7,3,23) and (0,1,0,11,3,7).
Each class is provably essential from either normalized row
(by inspecting its residues covered by no other class), and all
24 phase translations are distinct because of the period-24 entry.

For a fixed sign, the two patterns give precisely 48 CRT classes
of coefficients modulo 5592405, the product of the six primes.
For the minus sign their generating representatives before exponent
shift are 3968308 and 3177553; multiply either by any power of 2
modulo 5592405, taking the odd representative when needed.
The plus-sign classes are their negatives. If a prime cannot
divide a term, padding its empty phase with an arbitrary phase
could only make coverage easier; the essentiality thus also shows
that there are no dropped-prime covers missed by this enumeration.

**Brier obstruction for this wheel.** No coefficient can cover
both signs using only these six primes. Both covers would have to
use 7, yet 2 has odd order 3 modulo 7, so −1 is not in its power
subgroup: for a given coefficient, the same prime 7 cannot ever
occur on both signs. This is a limitation of the wheel, not an
obstruction to Brier numbers with different or larger wheels.
The period 24 here is a multiplicative-order clock; no Leech
lattice construction is needed to explain it.

## Why the 27→53→105→209 ray is different

Write A_e=13·2^e+1 for e≥1. In the 24-wheel its uncovered phases
are {2,4,8,10,16,20} mod 24. An uncovered phase is not a claim
of primality: A₂=53 and A₈=3329 are prime, but A₄=209=11·19.
Nor does every covered phase yield an automatic tournament factor:
prime 7 cannot be used as an ordinal-sum component. Indeed
7|A_e iff 3|e; 9|A_e iff e≡1 mod 6. Thus NO term in this ray
is divisible by 63, even though 63 itself and every 7·3^j for
j≥2 are realizable by the separate ×3 construction.
In the odd exponents divisible by 3 the 3-adic valuation is
exactly 1. At 7 the next two conditions are
49|A_e iff e≡6 mod 21 and 343|A_e iff e≡69 mod 147. Those phases
do allow a factor 49 or 343, but then one must still realize the
actual cofactor. The distinction is precisely between divisibility
on a finite clock and a constructive inverse at the absolute target.
