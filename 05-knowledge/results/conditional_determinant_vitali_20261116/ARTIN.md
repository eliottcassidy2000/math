# Artin and GRH really control the modular 67/69 problem, not the absolute spectrum

Use the two literal arms of the
[67/69 third-module construction](../chanpak_pair_67_69_20261114/THIRD_MODULE.md)
with the first two marked modules singleton and third module TT_r.
For every r>=1 their counts are

    A_r = 34·2^r−1,      B_r = 32·2^r+5.

In particular r=1 gives 67/69. At r=2 the ordering reverses:
A_r=135 and B_r=133=7·19. The modular question asks which odd
moduli M see *all* residue classes in the union of the two rays.

## Exact answer

For odd M>1 the union is all of Z/M precisely when

    M=p^k, p an odd prime other than 3 or 17,
    and ord_(p^k)(2)=φ(p^k).

For the forward direction when 2 is primitive, the sets are
−1+(Z/M)^× and 5+(Z/M)^×. Their complements consist respectively
of the −1 class and the 5 class modulo p, which are disjoint
since p does not divide 6. If M is a prime power but 2 is not
primitive, each ray contributes at most ord_M(2)≤φ(M)/2
residues, so the union misses something. Modulo 3 both rays
omit their common translate 2. Modulo 17, the order is 8 and
the first ray is constant, so no higher power of 17 can work.

For a modulus with at least two distinct odd prime factors, the
orders at the distinct prime-power factors divide their local
totients, which are even. Hence ord_M(2)≤φ(M)/2 by the least
common multiple bound. Even the union of both full orbits misses
some residue. This argument does not require their coefficients
to be invertible.

Thus the statement that infinitely many prime moduli are fully
covered by these particular *actual tournament counts* is
equivalent to Artin's primitive-root conjecture for base 2.
Under the Dedekind-zeta GRH used in [Hooley's theorem](https://doi.org/10.1515/crll.1967.225.209),
the prime moduli with this property have positive relative
density (Artin's constant). This is a genuine famous-conditional
theorem about a modular projection; it is **not** an implication
between GRH/Artin and the full Chan–Pak conjecture.

At 67 the order is 66 and 2^66=671 modulo 67², not 1, so
primitive-root lifting gives coverage at all powers of 67.
At 69=3·23 the order is 22 and the union has exactly 32 residues
out of 69 (22 in each arm with 12 in common). The 67/69 names
happen to distinguish entirely different modulus behavior here,
for explicit order reasons. Neither the Riesel wheel nor a
change of sign makes the arm values fill absolute intervals:
there are only O(log X) values across both arms through X.

Replay: `python check_artin.py` checks the order classification
through M=257, the two specific moduli, and the lifting witness.

## Twin-prime contrast: a different 67/69 arm has a permanent wheel

Instead inflate only the first module. The
[other inherited two-scale chart](../chanpak_pair_67_69_20261114/TWO_SCALE_AND_PORTS.md)
gives the consecutive physical pair

    (32·2^n+3, 32·2^n+5),   n≥1.

No pair on this ray consists of two primes. This is an exact
period-12 obstruction, using the composite member indicated:

| n condition | member | proper prime divisor |
|--|--|--:|
| odd | upper | 3 |
| 0 mod 4 | lower | 5 |
| 2 mod 12 | upper | 7 |
| 6 mod 12 | lower | 7 |
| 10 mod 12 | upper | 13 |

These conditions exhaust n. The divisibilities follow from the
orders of 2 modulo 3,5,7,13 (2,4,3,12); each entry is larger
than its prime divisor. This wheel preserves a physical +2
edit at all heights while forbidding a twin-prime pair. The
Twin Prime Conjecture, even if true, cannot supply its pairs
from this particular family. This is not a claim about other
physical +2 seams or general infinitude of prime gaps.

Two independent axes can escape that wheel. In the same inherited
table their transitive sizes (8,3) give counts 1,734,041 and
1,734,043. Both are prime by trial division (checked below), so
this is a literal twin-prime count pair outside the one-axis
obstruction. No infinitude assertion is made.
