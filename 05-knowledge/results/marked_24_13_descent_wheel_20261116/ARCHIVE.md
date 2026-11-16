# 24 and 13: a ledger of real transfers and false ones

| Roles of 24 and 13 | Concrete connection | What does not transfer |
|---|---|---|
| The wheel at ambient 24 | 13−1=12 divides 24, and 2 has order 12 mod 13 (2^6≡−1, with 2^4≢1). This one label gives one residue mod 12, two positions in Z/24. The three non-dyadic labels, 7,13,241, all have odd order part 3. | No primitive order-six label exists. The order-12 label cannot be assigned twice. |
| Bernoulli and pure-A Niemeier | The same 12 divides 24 in two other senses: p=13 contributes to the denominator of B_24, and A_12² has rank 24 and Coxeter number 13. The complete p−1 divisibility core of 24 is {2,3,5,7,13}. | The Leech lattice itself is rootless, and the agreement supplies no lattice action on a Riesel covering. |
| The 120 expansion | The p−1 divisibility core of 120 is {2,3,5,7,11,13,31,41,61}. Fermat gives their orders dividing 120 for odd p. | Their simultaneous phases and two-sign constraints do not follow from Bernoulli denominators. |
| Pythagorean doubles and triples | (5,12,13) doubles to (10,24,26), putting 24 beside twice 13. The primitive triples (5,12,13) and (7,24,25) are consecutive m=2,3 instances of (2m+1,2m(m+1),2m(m+1)+1). | 13 and 24 cannot be legs together, or leg and hypotenuse, of an integral right triangle: 745=13²+24² and 407=24²−13² are nonsquares. |
| Fibonacci boxes | The same neighboring triples have Euclid parameters (3,2), (4,3) and box quadruples (1,2,3,5), (1,3,4,7). The actual input is a marked sum/difference witness. | A width recursion is not an exact count recursion merely because the profiles sit on one ray. |
| Square-sum prefixes | The archived census puts 13 and 24 in the same set {13,14,19,20,21,22,24} of low double-connection vertices that do **not** themselves yield a Hamilton path. At 24 and 25 the twice-square graph displays Pythagorean component lengths 3,4,5, but extra components remain. | Enough incident edges or compatible component sizes need not connect the selected endpoints. |
| Monster work | In the archived Monster/Baby Monster ratio the exponent of 13 is 2; the level-6 Ogg lane gives Fi22 primes {2,3,5,7,11,13}; the separate involution lane has centralizer 2^{1+24}.Co1, whose Co1 factor also contains 13. | The archived prime-to-valuation cycle 4→2→5→3→7→13 is not induced by a congruence quotient modulo 9 or a tape evolution. |
| The actual marked tape | The simultaneous tests 24L_a=13L_b and 24R_a=13R_b give zero cross determinants upon gluing. Retaining the two opposite carries reads the gap exactly. | A tape’s width, residue and gap do not determine its starting Hamilton count. |

Archive pointers:

- [Wheel vs Bernoulli/rank-24 comparison](../two_three_wheels_synthesis_20261116/WHEELS.md)
  and [two-sheet follow-up](../unit_shields_wheel_followup_20261116/wheel/WHEEL_FOLLOWUP.md).
- [Leech deep holes and Coxeter frames](../leech_deep_hole_swap_theorem_20261022.md);
  [lattice/repunits note comparing {7,21} and {8,24}](../lattices_modular_repunits_20261021.md).
- [Square-sum and Monster ledger](../monster_levels_square_sum_7_21_20261022.md).
- [Fibonacci boxes and the K4 matching triple](../../../07-reflections/fibonacci-boxes-farey-and-the-k4-matching-triple-opus-20260728.md);
  [primitive-triple contraction and inverse guards](../ternary_triples_20260925.md);
  [Gaussian prime allocations as marked tournament endpoint fibers](../endpoint_pythagorean_fibers_20261019.md).
  The last link gives an actual marked bijection for one endpoint structure,
  not a rule translating arbitrary tournament counts to hypotenuses.

Riesel covering starts from a progression which already implies a
fixed divisor. A Chan–Pak wheel needs a stronger label: a progression
plus a physical inverse valid for **every sufficiently large absolute
target in that chart**. The ray series gives exact witnesses at every
sufficiently large even width for a chosen gap. Substituting “the
count is on the marked ray” for the absolute inverse would discard
precisely the missing coordinate.

## A Pell analogy with an exact failed transfer

The marked pair happens to satisfy 24²−3·13²=69.
The norm-69 Pell sequence under multiplication by 2+√3
contains

    (9,2), (24,13), (87,50), (324,187), (1209,698), …

with the larger coordinate first. The analogy is natural after
the repo's [norm-junction work](../chanpak_waves_20261113/README.md),
but the norm is blind to the binary carry guard.

Indeed for any marked ray (B/A=b/a) with (b>a>0), (b)
odd and (a) even, the generalized carry arithmetic in RAY uses
digits (0,-a,b,b-a), with forward carry strictly below (a).
When the reverse carry is zero, its parity requires an odd digit;
the only possible letters are 2 (digit (b), forbidden by the
forward bound for nonnegative forward carry) and 3 (digit (b-a)).
The latter leaves reverse carry zero and strictly increases the
already nonnegative forward carry. From the origin a nonempty
return is therefore impossible. In this Pell sequence parity
alternates because the unit update is (b,a) ↦ (2b+3a,b+2a).
In particular the ray at (87,50) fails, even though its norm
equals that of the successful (24,13).

Conversely, the following 86-letter tape on the ray (187,324)
satisfies both address equations and has zero signed gap, by
direct contraction:

    13221321333110313000202102220231121001012103212323321220131113213133302032200021021102

Its forward and reversed complement vectors are equal, so it
is reverse-complement invariant. This proves neither a nonzero
unit edit at that ray nor preservation of realizability on the
remaining parity subsequence.
