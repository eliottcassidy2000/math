# Complete retained-source extension tests through three births

Status: **PROVED** Boolean polynomial degree bound; **FINITE-EXACT** complete
extension census for the stated labelled sources. No exclusion for four or
more births, or for arbitrary H55 parents, is asserted.

## B1. Retaining the actual source

The previous [inverse65/INVERSE I5](../inverse65_zeckendorf_20261103/INVERSE.md)
completely classified the bank's six H55 source presentations, all on seven
vertices, and excluded111 after one arbitrary new vertex. Literal isomorphism
maps now group them as {0,1,4} and{2,3,5}; the supplied processing ports remain
different data. The successful external H55 source is a third ordinary class.

We hold every old arc fixed and vary EVERY arc incident with the new vertices,
including the arcs between new vertices. With r births there are
`D=7r+binom(r,2)` binary variables. A bit1 means lower-index->higher-index.
No strong-connectivity, endpoint, chart, or bank-membership filter is imposed.

| Source universe | Births | Orientations | Number with H111 |
|---|---:|---:|---:|
| All six labelled bank presentations |2|6*32768=196608|0|
| Two ordinary bank representatives |3|2*16777216=33554432|0|
| Known external H55 positive control |2|32768|4|
| Same external positive control |3|16777216|18|

The recorded isomorphisms transport the complete three-birth exclusion to all
six presentations while fixing the names of the new vertices. The one-birth
failure is replayed independently. Consequently **any induced H111 extension
of any of these six sources needs at least four new vertices**, if it exists.
This is an extension barrier for supplied sources, not a global H111 gap.

## B2. A degree bound that makes the finite audit cheap

For any fixed n-vertex tournament and r new vertices, H as a function of the
new orientation bits is a multilinear polynomial of degree at most
`min(2r,n+r-1)`. For each candidate Hamilton word, check every adjacent old-old
pair against the retained source. A compatible word then prescribes t arcs
incident to new vertices; t<=2r because each new vertex has at most two
neighbors in the path. Its indicator is

    product_(i in on) x_i product_(j in off) (1-x_j).

Sum these indicators. This proves the degree bound directly and retains all
mixed new-port terms. It refines the use of the inherited full Hamilton path
indicator; no claim of first discovering Boolean path expansions is made.

`two_birth.cpp` computes the entire value table in two bases:

1. Expand the indicator in ordinary monomials. A subset zeta transform gives
   H(mask) for every mask from these coefficients.
2. Put chi_i=(-1)^x_i. The indicator is a product of `(1-chi_i)/2` for on and
   `(1+chi_i)/2` for off. Multiply every coefficient by2^(2r). The coefficient
   of chi_S is `2^(2r-t)(-1)^|S intersect on|` when S lies in the used arcs.
   A Walsh-Hadamard transform gives exactly2^(2r) H(mask).

The two resulting integer arrays are compared at **every** orientation, not
just by a hash or histogram. Their common path-constraint input is additionally
checked by an independent subset-DP implementation on585 literal extensions,
comparing the exact count at each sampled mask in both birth universes.
Python independently reconstructs every two-birth value from the published
sparse monomial coefficients and recovers the complete histograms.

For the largest case D=24 and n+r=10, signed64-bit arithmetic is safe:
even the crude bound `2^24 * 2^6 * 10!` is below2^63. The code also verifies
the degree bound, positive odd outputs, complete mask totals, and all named
positive witnesses. Permutation enumeration is factorial in n+r and the
transforms are exponential in D; this is an exact small-r tool, not a
polynomial-time algorithm for general tournaments.

## B3. Relation to the recurrence carrier

The source-replacement constructions in THEORY deliberately change a bounded
old source block. They are therefore not contradicted by this retained-edge
extension result. The distinction is now testable rather than verbal:
the six H55 representatives fail even all three-birth choices, while another
actual H55 source succeeds with one birth. Neither H nor a weighted mean of
its children can identify which source is available for the next operation.

The81/27 graphs are also replayed: identical67 labelled cyclic triples and
H462915, but different responses2390013 and1754919 after the same sink and
named-port operation. These checks preserve the hostile source, not merely
its automorphism-order labels.
