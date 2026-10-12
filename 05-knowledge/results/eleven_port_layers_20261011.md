# Eleventh-root generators: an exact port relation and the next layer

2026-10-11. **PROVED, scoped and self-audited:** the universal one-layer
boundary relation and the all-length generator colouring below. **FINITE-EXACT:** a complete8436-vertex mixed-generator
unit graph and an explicit three-colouring. **OPEN:** six-colour forcing in
Q(zeta55), independent ROOT positivity and universal Collatz. This does not
settle the chromatic number of the field or the plane. No peer audit or
formalization is claimed.

## Inheritance, and why the construction mechanism must change

Put `L=Q(zeta55)`, `alpha=zeta5`, `gamma=zeta11`, `beta=sqrt(-11)`, and
`K=Q(alpha,beta)`. Then `[K:Q]=8`, `[L:K]=5`. The original pentagon P and
rhombus R are fixed pointwise:

    P_j=sum_(i=0)^(j-1) alpha^i,  j=0,...,4,
    eta=(-3+sqrt(-55))/8,
    R=(0,-eta,3/4,-bar(eta)).

The nearest mechanisms are [G2 and the odd-degree barrier](rhombus_overlap_elliptic_20261011.md),
which give `chi(K)<=5` and prohibit a quadratic circle-intersection procedure
from leaving K inside L; [the unit component](cyclotomic_unit_component_20261011.md)
proves that the full L graph is triangle-free and identifies its exact
localization. The stronger [pentagon list classification](pentagon_context_number_fields_20261011.md)
already says that a five-colourable external context cannot supply the bad
six-colour boundary lists. These statements are retained, not replaced by
weaker finite tests.

The hostile is genuine field escape with a surviving three-colouring.
The corrected near miss is expecting a larger bank of eleventh-root unit
directions to force a new boundary relation. The least-used sidecar is now
*which preserved port can be a unit neighbour of each new point*.
The board and anchor/niche/wildcard transfer contract are in the
[charge-receipt companion](collatz_charge_receipts_20261011.md).

## P1. A new point can touch at most one old K-port

If `x in L\K` and distinct `p,q in K` both satisfy
`|x-p|=|x-q|=1`, put `t=(x-p)/(q-p)`. Subtracting the two norm equations gives

    t+bar(t)=1,   t*bar(t)=1/N(q-p),
    t^2-t+1/N(q-p)=0.

Thus `[K(t):K]<=2`, but it divides `[L:K]=5`. Therefore t and x lie in K,
a contradiction. This is a consequence of the previously proved circle
barrier, now used as an attachment restriction. It does not say that new
points cannot touch one another or cannot touch two ports that are themselves
outside K.

The full unit graph on L is triangle-free: an equilateral unit triangle
would put a primitive cube root of unity in L, whereas adjoining it doubles
the cyclotomic degree from40 to80. Consequently the unit neighbours of any
one point form an independent set.

## P2. Every one-layer pentagon context has exactly the bare six-colour relation

Let B be any set of at most5 distinct K-points. Let X be any set of points
in `L\K`, each having a unit neighbour in B. Retain **all** actual edges on
`B union X`, including every cross-edge and coordinate identification.
Then every proper six-colouring of the induced graph on B extends to the
whole context. This is valid for finite or infinite X.

**Proof.** P1 assigns every x a unique neighbour p in B. Partition X into
these port groups X_p. Each X_p is independent by triangle-freeness.
For each p choose a colour d_p different from its boundary colour, with all
d_p distinct. There are five allowable colours per port and at most five
ports; choosing successively always works. Colour all of X_p by d_p.
Edges to B are proper by construction; an edge inside one X_p is impossible;
all cross-group edges have distinct colours. No crossing edge is dropped.
The reverse inclusion of boundary relations is restriction of a colouring.
QED.

Taking B=P gives the **authenticated exact relation**

    R_6(context; P)=R_6(C5; P).

It has `5^5-5=3120` labelled assignments. The checker independently enumerates
all3120 and supplies the distinct satellite colours for each. The proof
allows arbitrary eleventh-root generators and arbitrarily many satellites;
this is not a bounded search result.

In particular these one-layer contexts admit a colouring in which all
neighbours attached to one pentagon vertex share one colour. They cannot
force the four distinct forbidden colours required by that vertex's bad
list. Merely changing a generator or adding more points on those five unit
circles cannot repair the problem.

The same proof works with k colours and at most k-1 ports. With exactly k
ports it works whenever their boundary colours are not all the same:
Hall's condition holds because each allowed list omits just one colour,
and the union of all k lists is all k colours in the nonconstant case.
This generalization is about the same attachment hypotheses, not arbitrary
list-colouring problems.

A potentially productive candidate must therefore add vertices **not
adjacent to any original port**, move additional ports into the extension,
or use a larger retained port interface. These possibilities are not claimed
to force six; they escape the specific theorem. The fixed pentagon and
rhombus together have eight distinct ports, so P2 alone does not classify
that larger boundary relation.

## P3. Mixed generators and an exact three-layer test

The first probe used norm-one directions

    u_(a,j)=(a-gamma^j)/(a-gamma^(-j)), a=2,3,

rather than only roots multiplied by elements of K. These remain in the
smaller field Q(gamma); field membership alone cannot identify a forcing
mechanism. One- and two-layer finite probes remained three-colourable.
Those exploratory sizes are not used as a proof or filed as an exhaustive
census.

The maintained test makes the two cyclotomic directions interact inside a
single unit:

    u=(1-alpha-gamma)/(1-alpha^(-1)-gamma^(-1)).            (1)

The denominator is nonzero: otherwise two unit vectors would sum to1,
creating a forbidden unit triangle. Equation(1) gives N(u)=1. Its exact
40-coordinate power-basis denominator is60002251, and the automorphism
fixing K and sending gamma to gamma^3 does not fix u. Thus u genuinely lies
outside K. The exact vector and this test are retained in the output.
There is the six-unit relation

    1+u*alpha^(-1)+u*gamma^(-1)=u+alpha+gamma.             (2)

This is an actual additive relation in L, with its norm-one factors kept;
it is not inferred from the common number six.

Let D be the symmetric set consisting of all alpha powers, all gamma powers,
and `u,u*alpha^(-1),u*gamma^(-1)`, together with their negatives. Form all
sums of at most three elements of D, then adjoin the original P and R
without moving them. Authenticate every pairwise unit distance.
The completed graph has

    8436 distinct vertices, 24238 unit edges, chi=3.

The lower bound is the retained pentagon. The upper bound is an explicit
colour vector checked on every actual edge. The graph has second- and
third-layer points, so it is a genuine test beyond P2, but it supplies no
six-colour forcing context.

The exact geometry uses integer vectors in `Z[z]/Phi55(z)` after clearing
a common denominator240009004. Evaluation modulo331 and661 rejects pairs
that cannot have unit norm. These are ring homomorphisms, hence they cannot
reject a true edge. Every survivor is checked by full power-basis norm
multiplication; there are24426 survivors and24238 true edges. A separate
small complete check tests every marked-port pair without the filters.
A search using only the declared generator edges would have counted24234,
so retaining all distances adds four edges. Those four edges are included
in the final colouring check and its graph hash.

[Reproduction data](charge_ports_eleven_20261011.json) stores the full colour
vector, all eight port indices, exact generator coordinates, edge and
coordinate hashes. The deterministic generator reconstructs every point
and edge. Normal and optimized runs must agree byte for byte.

## P4. A 40-trit certificate screens the entire declared generator graph

The finite test has a stronger explanation. Let `A_(3)` be the elements of
`Q[z]/Phi55(z)` whose40 rational power-basis coefficients have denominators
coprime to3. Reduction of those coefficients defines an additive map to
`F3^40`. Take the dot product with

    ell=(1,2,1,0,2,1,1,1,0,0,1,1,0,0,0,1,0,0,0,0,
         1,0,1,0,0,1,0,0,0,0,1,0,0,1,0,1,0,0,0,0).

Exact reduction gives nonzero values on every direction in

    D' = D union {eta,bar(eta),-eta,-bar(eta)}.

Therefore `x -> ell(x) mod3` properly colours **the whole Cayley graph on
A_(3) with connection set D'**, at every word length. This is a direct
finite arithmetic certificate for an infinite graph. It does not rely on
a completeness theorem for additive colourings.

All24238 actual edges in the completed finite graph have their directions
in D'. The checker authenticates that statement and independently verifies
this additive colouring on every edge, in addition to the SAT-produced
colour vector. Thus the failure is explained, not merely observed.

The distinction is essential: the unit graph on A_(3), or on the full field
L, can have directions outside D'. No upper bound of3 on either full unit
graph follows. Any construction defeating this particular colouring must
introduce an actual norm-one direction v with `ell(v)=0`, or leave A_(3).
Such a direction is necessary to defeat the displayed certificate, not
sufficient for four or six colours. This is a cheap next-generator test
before another large point search.

The cheap hostile probe succeeds immediately: `v=zeta55^3` is a unit root
with `ell(v)=0`. Yet the full integral unit graph already has a different
three-colouring, proved in the inherited unit-component note. Killing one
functional is therefore easy and gives no chromatic lower bound. The next
substantive test must exclude *every* extension of a specified joint port
assignment, with an actual geometric context; one failed colouring is only
a search hint.

## P5. What the other work contributes, and the next tests

The concurrent [multiport spindle analysis](spindle_tournament_interfaces_20261011g.md)
shows why pairwise projections can miss a joint boundary constraint. Here
P2 goes the other way: it gives a complete extension for every joint
assignment, so a hidden higher-order obstruction cannot be supplied by
this one-layer class. The [diagonal palette analysis](diagonal_palette_holonomy_20261011g.md)
keeps all parallel seam equations; the full-distance completion in P3 applies
the same discipline to actual geometric cross-edges. The
[finite-field escape family](unit_field_escape_recursion_20261011f.md)
remains a warning that field escape does not imply more colours.

We also read the pinned external working notes
[three colours over number fields](https://github.com/decalion89/chromatic-number-of-the-plane/blob/cb660f322dce47a773f347d39818ffc46e6ca72c/notes/three_colours_number_fields.md)
and its winding-lemma argument as exploratory inspiration. No general
classification or new literature claim from those notes is a dependency
here. A plane `F^2` and the norm graph on our CM field L are different
carriers; an additive three-colour test is not a completeness theorem for
arbitrary colourings. For example the earlier rational-fifteenth connector
already prevents an additive F3 colouring on its generated group while its
quadratic-field graph still has a three-colouring by reduction at its
ramified prime11 and a colouring of the resulting11-cycle.

The next geometric test should first defeat P4's exact colour functional,
then target a complete nontrivial boundary
relation involving unanchored interior points, and minimize that relation's
support. Growing a graph that still admits an explicit three-colouring
cannot supply the six-colour context. The next arithmetic test should expose
a new actual join past a proved inverse shield, as at39; increasing precision
alone does not pay the fixed source's odd charge.

| Connection | Exact map / preserved predicate | Destruction and required sidecar |
|---|---|---|
| Odd field degree -> attachment | Two norm equations -> degree<=2 | Does not constrain all new-new edges; keep all incidences |
| Attachment -> complete boundary relation | Assign one distinct colour per independent port group | Requires every new vertex be anchored and outside K |
| Mixed unit -> geometric graph | Six-unit relation and actual coordinates | A relation is not a forcing proof; keep complete edges and lifts |
| Inverse shield -> port-layer analogy | Fix the interface, classify every allowed continuation | No Collatz/geometry conjugacy; new terminal or new layer separately checked |

No change to the closed genus-two carrier is asserted. No rational
Hamiltonian, odd-zeta irrationality, or whole-plane conclusion follows from
these port and charge constructions.
