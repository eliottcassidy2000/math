# Two compatibility tests, not a literal plus construction

The historical [stellohedron transport](../session285_five_web_20261112/holonomy/PROOF.md)
separates closed increments from a chosen absolute integration constant.
The [graceful star buffers](../entry_resources_20261113/GRACEFUL.md) separate
a legal exchange from adding the new crossing contact that makes it legal.
In the [marked Collatz Bott experiment](../collatz_bott_marked_periodicity_20261007d.md),
similar affine maps can act on disjoint native residue cylinders. These three
precedents give a useful common rule: **retain the witness to compatibility
when composing, not just the invariant on which it agrees.**

## The exact categorical test

Model a reservoir transition as a span of sets A <- E -> B, where E consists
of literal moves and A,B retain the data needed to use them. Composing it
with B <- F -> C requires the pullback E x_B F: the arrival and departure
must be the same state. If we first project B to absolute Hamiltonian count,
we instead form the larger E x_H F. The map from the former to the latter
is injective and need not be surjective. It is surjective on a pair exactly
when its endpoints match in B. More generally, if neutral repair moves
on B are allowed, replace E x_B F by triples (e,repair,f) with that
literal repair between the named endpoints. Merely having a nonempty
equal-count fiber does not guarantee a repair; even two endpoints of a
true strong exact switch need not have the same continuation menu.

The sheaf analogy has two *different* obligations. The first is to identify
local descriptions that represent the same state upon refinement; the
second is to glue matching local data into an actual global object. In the
[Stacks Project plus construction](https://stacks.math.columbia.edu/tag/00W1),
one plus makes a presheaf separated and a second makes it a sheaf. The
span observation specifies what would have to be checked for a similarly
structured tournament argument. No site or sheaf theorem about the
Hamiltonian spectrum has been constructed. Likewise taking a second
difference detects interaction of two cuts, but cannot supply the
missing anchor: all affine functions have the same zero second difference.

## Arithmetic coverage is weaker than reusable gluing

For every witnessed triple (a,a+2,a+4), with a odd, put a vertex at a
and join starts at distance at most 6. Two such triples at distance 6
do not intersect, but their union has no missing odd integer between them.
A component has an interval as its union on the odd lattice: each edge
joins overlapping or adjacent integer intervals after rescaling, and any
two vertices have a finite chain between them. In particular a component
starting at 23 and unbounded upward would prove eventual coverage from 23
(with the separately known smaller values). Conversely, if every odd
integer from 23 onward is covered by such triples with starts >=23,
all those triples are in a single component. A missing odd x is a cut:
no edge can join a triple entirely below x to one entirely above x.
The edge definition here certifies **arithmetic** coverage, not a reusable
move. A physical transition graph can be a strict subgraph. Any unbounded
component of that subgraph still suffices, so a useful strategy is to
prove the stronger property using interfaces, not to assume the two
graphs are the same.

The parallel with Collatz is exact only at the level of the obligation:
one needs native entry and a well-founded descent, not just an expanding
collection of admissible words. For accelerated Collatz, proving that
each odd n>1 reaches a strictly smaller positive odd integer in finitely
many legal steps suffices by strong induction. For Chan–Pak, it suffices
to give, for every proposed target above a finite certified base, a strictly
smaller target with a certified sufficient interface and a literal lift
to the required interface at the original target. These induction steps
would both be universal; they are not logically equivalent. In particular
the sparse exact switches in [Dyadic switch](DYADIC_SWITCH.md) are not a
proof of either universal entry.
