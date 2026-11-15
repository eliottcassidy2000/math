# A logarithmic width barrier for localized absolute repairs

This obstruction covers single-site moves, two-site moves, and arbitrary
interleaving *inside a bounded interval*. It does not restrict how the
external addresses are chosen from the target.

**Theorem.** In the one-apex transitive-spine family, count endpoint
values ≤X that differ from another value by a nonzero even gap ≤2G
via a change of bits contained in an interval of at most w sites.
The number of such different endpoint values is at most

    O( G [ 4^w + sqrt(X+1) · 2^(3w/2) ] ).

The constant is absolute. In particular bounded-width unit
repairs have zero density, and even width at most
(1/3−ε)log₂X has zero density for every fixed ε>0. This is
a statement about incidence to an exact repair, not merely about
the size of the full one-apex count spectrum.

## Proof

For width m, fix a pair of distinct m-bit interior words.
Use the [exact two-cut seam](INTERLEAVED.md) with left address a,
right reversed-complement address b, and interior addresses u,v.
Once the two exterior parities are fixed, the change in H is

    constant + a Δv + b Δu.

**Both** coefficients are nonzero: the ordinary and the
reversed-complement integer readings are each injective on a
fixed-width word. This holds even when the first or last bit of
the interior changes. The parity terms are constants on each
exterior parity sheet. Thus each signed nonzero prescribed
gap and pair of words imposes an affine line on the actual
positive or zero integer addresses. For either address fixed
there is at most one possible value of the other address.

The same seam formula gives H≥2^m ab: all its other terms are
nonnegative after grouping H(a)+H(u)+H(b)−2, which is at least one.
If one endpoint is ≤X, at least one of a,b is at most
sqrt(X/2^m). Count the potential solutions by the smaller
address, including the two zero boundaries. The number is
O(1+sqrt(X/2^m)), uniformly in the sign and size of the
affine coefficients and with no need for injectivity of H
along that line. Padding representations of a,b with zeros
does not add any counts to the seam formula.

There are at most 4^m ordered pairs of words, four exterior
parities, and 2G signed nonzero gaps. Sum the two geometric
series for m≤w. This yields the claim. It also shows when allowing
small but nonunit steps makes a difference: if G=O(X^γ), then
width w≤(1/3−2γ/3−ε)log₂X still gives o(X) for 0≤γ<1/2.
The full interval-coverage strategy using only these short
moves cannot even make a dense set of vertices incident to an edge.

The [one-site inverse](SINGLE_SITE.md) gives a matching
Θ(sqrt(X)) fixed-gap estimate at w=1. At growing width this is
an upper bound, not an assertion of sharpness. Changing the
physical reader, changing bits at a span comparable to the
spine, or combining non-additive large moves is not excluded.
