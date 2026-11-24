# Diagonal fugacity in fixed-point-free self-converse tournaments

Let T have 2m vertices paired by an involutive antiautomorphism
α with no fixed vertex. The cube holds every nonmate arc fixed
and independently toggles the m arcs within mate pairs. Every cube
member is self-converse with the same α. This applies to the
uncompressed Moore observers with sources u_i, observers w_i
and mates u_i,w_i. Toggling diagonals of those observers need not
retain the original Moore Gram.

Let H(z), C(z) be the integer multilinear polynomials on the cube
counting Hamilton words and directed Hamilton cycles up to rotation.
A literal weight realization puts weights z_i and 1−z_i on the
two directions of mate arc i, and weight 1 on present nonmate arcs.
Each word or cycle uses at most one of the two directions of a mate
arc (take m≥2 for cycles), so these weighted sums are multilinear.

## An exact determinant formula and a path-count congruence

Write M(z) for the weighted directed adjacency. In the commutative
squarefree ring on vertex variables x_v with x_v²=0, put D=diag(x_v).
Then

    H(z) = [∏ x_v] (∏(1+x_v)) det(I+D M(z))/det(I−D M(z)).

To prove it, take the trace logarithm of the determinant quotient.
Each term in tr((DM)^l) involving a repeated vertex vanishes.
The l cyclic starts cancel the denominator l; subtraction between
the signs leaves precisely twice the sum of weighted simple odd
directed cycles. Exponentiating gives the generating function of
vertex-disjoint collections of such cycles with weight 2 per cycle,
since overlapping cycles annihilate and the factorial cancels the
ordering of a disjoint collection. The factors (1+x_v) allow unused
vertices. At Boolean z, the exact formula for Hamilton words as
the sum of 2^(number of cycles) over odd directed cycle packings is
[Grinberg–Stanley, arXiv:2307.05569, equation (42)](https://arxiv.org/abs/2307.05569).
Both sides are multilinear also before Boolean specialization:
a disjoint collection of simple cycles of length at least three
cannot use two orientations of a given mate arc. Thus the identity
holds on the entire formal diagonal fugacity cube. This formula
uses an exponential-size squarefree algebra; it is not an
ordinary polynomial-time determinant algorithm for H.

The reversal induced by α has no fixed directed odd cycle: an
invariant odd vertex set cannot be a union of α-pairs. Let O count
its two-element orbits of odd cycles, and J count those orbits
whose two cycles are vertex-disjoint. Modulo 8 only singleton
and doubleton cycle packings contribute. There are 2O singleton
packings; a doubleton packing is fixed under α exactly when
its cycles are a disjoint α-pair, so the doubleton count is J mod 2.
Consequently

    H(z) ≡ 1 mod 4,
    H(z) ≡ 1 + 4(O−J) mod 8                 (for Boolean z).

An orbit contributes to O−J iff its cycle uses both vertices
of at least one mate pair, whether or not their mate arc appears.
The formula implies every nonconstant coefficient of the integer
multilinear H polynomial is divisible by four, by finite differences.
In particular no mate flip can yield a Chan–Pak unit step of two.

More locally, for a single mate i the parity of ΔH/4 equals the
number of anti-orbits of odd cycles that use either orientation
of the i-th mate arc, with all other arcs fixed. Such orbits are
intersecting automatically, since both cycles then use that mate
pair. Those avoiding the arc do not change. For this count one
may temporarily allow both directions of the flipped arc and
count odd simple cycles using either direction.

## Hamilton cycle parity is affine, with exact quotient slopes

The half-size quotient Q on mate labels is defined as follows.
For ordered distinct labels i,j, a block of arcs from pair i
to pair j is a two-by-two Boolean matrix of the form

    [[a,b],[1−b,1−a]].

Its 1 entries form either a full row (tail sign pinned) or
a full column (head sign pinned). Reversing the ordered pair
exchanges the two kinds. Put i→j in Q precisely when it is
tail-pinned. This makes Q a tournament.

Act on Hamilton cycles by applying α and reversing direction.
Nonfixed cycles cancel modulo two. In a fixed cycle, α acts as a
reflection with no fixed vertices, and therefore fixes two opposite
edges setwise. Those edges are mate arcs. If their labels are i,j,
cut just after the i-th mate arc and just before the j-th. The
resulting half-word runs from label i to label j and chooses one
sign per label. It is paired with its mirror half-word, so the fixed
cycle is recovered uniquely. Fixing either label as the first seam
does not introduce a multiplicity: each has one directed mate arc,
and the forward half between those two edges is unique.

For the coefficient of z_i z_j, sum over all signs in these
half-words, since the two orientations of each seam both have
coefficient 1 modulo two. For a fixed label order there are m−1
transition constraints on m signs and each transition fixes only
one sign. The consistent assignments therefore have even
cardinality, and inconsistent ones have zero. Thus every mixed
quadratic coefficient vanishes; there are no terms of larger degree
because fixed cycles have just two mate edges.

For the linear coefficient of z_i, consider each possible other
seam j, holding it at the baseline direction. Then the last sign
is fixed (tail of seam j), whereas the first sign is summed freely
(the derivative at seam i). To have an odd number of assignments
for this label order, the m−1 transition pins must address exactly
the m−1 labels other than j, once each. The last transition is
therefore tail-pinned, and working backward shows that every
transition is tail-pinned. Conversely an all-tail half-word gives
exactly one sign assignment for either fixed orientation at j.
Its order is precisely a Hamilton word of Q starting at i
and ending at j. Hence

    C(z) ≡ C(0) + Σ_i H_start(Q,i) z_i mod 2.

Since H(Q) is odd, the number of nonzero slopes is odd. Exactly
half the cube has odd cycle count. Complementing all mate bits
always changes cycle-count parity while keeping H in 1 mod 4.

For the Moore observer with increasing sources and decreasing
observers, Q has arc i→j for i<j exactly on Moore graph edges
(and the reverse on graph nonedges). Thus it is precisely the
all-nonedge-backward one-layer graph-to-tournament orientation,
with no native-spine requirement. The quotient slope is a
Hamilton start count on that one-layer orientation.

## Exact pentagon calibration

The checker runs all 32 mate assignments on each of the eight
decoded order charts, counting Hamilton words and directed
Hamilton cycles by independent endpoint subset recurrences.
It verifies the mod-eight cycle-orbit formula directly on all
32 assignments of the first chart, and the preceding path and
cycle polynomial congruences on all eight charts. In each
chart both residues 1 and 5 mod 8 occur for H, so the universal
mod-four statement cannot be strengthened uniformly to mod eight.
Among the 640 literal cube edges, exactly two have nonzero
height difference four. Both use the pentagon order 01342:

    mate mask 1 → 5:    4225 = 65²     → 4221 = 3²·7·67,
    mate mask 18 → 22:  8369 (prime)   → 8365 = 5·7·239.

In each case the flipped mate has index 2 (zero-based). The
minimum nonzero absolute change is therefore exactly four over
this calibration. The same flip changes cycle parity in both
examples. These examples give individual count transitions
and neither a general prime-to-7p controller nor an unrestricted
spectrum result.

Reproduce with `python3 moore-diagonal-fugacity-check.py`.
