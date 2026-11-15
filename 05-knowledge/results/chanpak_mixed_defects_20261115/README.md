# Exact mixed defects, not a larger construction tree

This session asked what the previous repairs actually preserve. The
answer is an exact interaction law on the marked choices, plus a
continuation condition for moving the whole repair to an absolute
target. The distinction yields positive constructions and all-size
obstructions. Chan–Pak remains open.

## A positive nonlocal repair

Two disjoint, interleaved controls on one strong 12-vertex tournament
give **721, 723, 725, 727**. Their mixed Hamilton-count difference is
exactly zero over the integers. This escapes the last session's
contiguous obstruction. The exact exterior response transports
this block as two +2 pairs to unbounded size, but not as a full
width-six window: the second step grows under that move.
[Interleaving, the mirror cut graph, and the continuation test](INTERLEAVED.md).
A smaller nonlinear example gives 213, 215, 217, 219.
[Target-side inverse and controls](NONLINEAR_ESCAPE.md).

## The unifying theorem: the degree is a physical constraint

If controlled arcs touch r distinguished vertices, the exact
Hamilton count in the switch variables is a polynomial of degree
at most 2r, with even nonconstant coefficients. Thus an additive
repair on a legal family such as matchings requires the mixed
differences on precisely its feasible higher-order faces to vanish.
At one apex, pairwise cancellation suffices. With two apices it
does not: an actual 7-vertex example has **all three base pair
defects zero and a third defect of −2**.
[The degree theorem and witness](INTERACTION_DEGREE.md).

This realizes the abstract continuation-dividend idea and makes
the earlier graceful matching proposal falsifiable. In one apex,
two-site switches have a complete zero-interaction classification;
even with nonzero interactions they cannot produce four consecutive
odds. If a unit repair of any shape changes only a short interval
of the spine, the count of incident endpoint values is power-sparse:
width below (1/3−ε)log₂ X still gives o(X) distinct counts ≤X.
[Two-site classification](TWO_SITE.md) · [Width barrier](SUPPORT_WIDTH.md).

## Why 7 fails: local overlap forces missing interactions

Two triangles meeting at a vertex force a mixed triangle and a
mixed five-cycle, not just the two original cycles. The same
argument yields a weighted product lower bound for *any* number
of arms and all cross orientations. It gives a new non-enumerative
proof that 7 is impossible. In contrast, the 63 witness has 31
intersecting odd cycles across lengths 3, 5, and 7: its numerical
factorization does not describe its decomposition into components.
[Overlapping-cycle bound](BOUQUET.md). A fixed one-arc flip in the
bouquet also realizes the two signed forms 6Q−H and 6Q+H across
the same arbitrary module. [One-arc repair](ARC_PORT.md).

## Why another local clock would not close the proof

A power-sparse ordinal-sum-closed subset of the **actual** spectrum
can retain every previous local guarantee: literal +2 choices on
both sides, every fixed modulus, growing precision, and any finite
census. It can still miss almost all primes. The mass of local
information is Θ(XM/R), for window radius R and modulus M, not
Θ(X). [Local-to-global obstruction](LOCAL_GLOBAL.md).

What remains is an exact lift in one physical residue fiber. If
opposite-side witnesses for N can be connected by real moves
preserving the residue M and changing the count by at most M,
the first crossing hits N exactly. The priority is a target-selected
reader with that property, not another collection of isolated
branches or a larger finite census.

## Replay

Run `python verify_bouquet.py`, `python verify_interaction.py`,
`python verify_interleaved.py`, `python verify_two_site.py`, and
`python verify_escape.py` from this directory.
