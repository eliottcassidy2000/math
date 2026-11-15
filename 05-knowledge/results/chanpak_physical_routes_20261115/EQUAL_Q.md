# A physical crossing that survives unbounded padding

There is a three-state route with two fixed two-arc moves which survives
unbounded padding outside its eight-vertex core. Unlike the one-apex
calibration, this route uses an exact fixed optional-slot interface. The
step size can be fixed or grown as a prescribed power of the absolute
count. This is an existence theorem for an explicit family of targets,
not a route for arbitrary bracketing witnesses.

## Constant-Q core and exact crossing

Write Q(D) for the two ordered optional path slots and H(D) for the
Hamilton path count. Use row-packed adjacency encoding on eight
vertices: bit 8i+j is 1 iff i→j. The following are literal cores:

| Core code | H | Q |
|---:|---:|---:|
| 1623487625516629146 | 645 | 8784 |
| 1695545228144487450 | 651 | 8784 |
| 1550304140166016026 | 657 | 8784 |

The first move reverses (0,7) and (1,4); the second reverses
(1,7) and (2,6). They have disjoint arc supports. These are physically
specified two-arc edits, not replacements by unrelated witnesses.
The verifier below computes H and Q by subset DP directly from
the encoded arcs. It does not rely on completeness of the
[existing order-eight inventory](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/n8_classes.csv)
for the resulting theorem.

Pad each core by ordinal sum with a t-vertex transitive tournament.
H is unchanged and Q becomes 8784·2^t: each of the two slots in
the ordinal suffix receives a subset of its vertices, and the
two fixed slot words concatenate independently. Feed the padded module
D into the plus side of the [older bouquet port](../chanpak_mixed_defects_20261115/ARC_PORT.md)
reading 6Q(D)+H(D). Only the two indicated arcs inside the fixed
eight-vertex core change on either step, regardless of t.
The ambient tournament is strong and the three counts are

    52704·2^t + 645,   52704·2^t + 651,   52704·2^t + 657.

They give an exact same-residue crossing at every t and fixed modulus 6.
For t=0 these are 53349, 53355, 53361. Q is a continuation invariant
with an integer reading, not merely a modular compatibility check;
it keeps the old step after arbitrary transitive padding.

## A polynomial-size moving modulus with the very same edits

The [old cyclic rank-two reader](../chanpak_bott_cayley_periods_20261115/CAYLEY.md)
says that the cycle of a module D, a module E and a singleton c
has Hamilton count

    (Q(D)−H(D))(Q(E)−H(E)) + 2 H(D)H(E).

In particular, if E is transitive of size s≥2, this is

    (2^s−1)Q(D) − (2^s−3)H(D).

So reverse the order of the same three cores and use this cyclic
reader in place of the bouquet. The route is strong because the
three nonempty blocks form a directed cycle. Every designated
count is congruent to the middle count N modulo M, and the steps are
exactly M:

    M = 6(2^s−3),
    N = (2^s−1)·8784·2^t − (2^s−3)·651,
    route: N−M → N → N+M.

For example, s=4,t=0 gives 123219 → 123297 → 123375, in one residue
modulo 78. No large edit hides inside the scale choice: each move
reverses the same two old arcs. Choosing t=159s gives M on the order
of N^(1/160); any prescribed fixed exponent between zero and one
can be obtained by analogous integer rounding. This is a
polynomial-modulus physical crossing with bounded edit size, in
contrast to a construction where the steps themselves must touch
more and more sites. The exact whole tournament, not just its
residue or an asymptotic approximation, is built at each vertex.

## What this does not supply

These targets belong to a particular exponential two-parameter family,
not all N or the old independent brackets. In this rank-two lane a
reader using any fixed second nonempty module E multiplies a constant-Q
difference by 3H(E)−Q(E), an odd integer: H is odd and Q even. Thus
rescaling cannot raise the 2-adic valuation of the input gap while Q
is frozen. For the explicit cores above v₂(M)=1, including the growing
modulus. This parity obstruction is intrinsic to this reader, not a
failure to select a large enough second module. To match the old
growing high-dyadic modulus, either produce source Q-fibers with growing
2-adic stride or vary more than one continuation coordinate.

There is a literal multiplicative-order question in the odd part of
the modulus: for p not dividing 6, divisibility p|M is equivalent to
2^s=3 modulo p, and when soluble the exponents form one residue class
modulo ord_p(2). Thus even a strong statement about 2's multiplicative
order could control which odd prime factors can be included. It would
not prescribe the absolute additive anchor N. The distinction is
between finding one orbit in a quotient and controlling the output
integer along the actual path.

## Complementary invariant: an actual growing high-dyadic route

There is also an explicit equal-H route by reversing the roles of the
two interface coordinates. The same row-packed convention gives:

| Core code | H | Q |
|---:|---:|---:|
| 1119081703992204326 | 589 | 8080 |
| 1119081703992269858 | 589 | 8104 |
| 974966515916576802 | 589 | 8128 |

Pad each by an ordinal transitive suffix of size t. Read the padded core
using the same cyclic reader with s=1, now giving Q(D)+H(D), since
Q(E)=2 and H(E)=1. The three counts are

    8080·2^t + 589 → 8104·2^t + 589 → 8128·2^t + 589.

The first move reverses arc (0,2); the second reverses (1,2) and (1,7). For example t=1 gives
16749 → 16797 → 16845 in one residue modulo 48. The modulus is now
exactly 3·2^(t+3), of the same arithmetic **form** as the earlier
high-dyadic bracket. But M is on the order of the absolute count in
this simple lane, not on the order of its 1/160-th power. So there
are two separate successes: a sublinear polynomial-size modulus
with dyadic valuation only one, and a high-dyadic modulus with a
sparse, much coarser pitch. Neither proves that one can join the
earlier independently chosen witnesses.

The underlying transport criterion is not “freeze the most
coordinates”. Freeze the coordinate whose readout forces the
wrong scale, and vary the coordinate whose response grows in the
desired direction. Keeping Q fixed forbids dyadic amplification
in this rank-two reader; keeping H fixed allows it via padding.
The next task is to find a reader that grows an absolute anchor
without dilating the pitch by the same amount.
