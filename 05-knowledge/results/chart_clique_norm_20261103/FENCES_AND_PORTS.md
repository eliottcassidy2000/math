# Shared boundaries force a change of graph class

**PROVED restricted-class bounds and exact constructions / CITED gallery
entries.** The two Friedman problems have different objectives and constraints.
No global optimum in either gallery is newly established here.

**Incoming precedence.** Public Session145 commit4c6f0ca135 contains an
independent proof of the same n/4 pseudoforest bound, a stronger exact
product theorem for all endpoint-only cacti, and a scoped optimum for
theta(5,5,2) with its two-edge separator and interior terminal chord.
See its `portfolio/fences/PROOF.md` and `theta_addendum/PROOF.md`, pinned
in SOURCE_PINS.json. Our proof and radical-coordinate replay below are
retained as independent deductions and checks, not claimed as a new
advance beyond that incoming packet. General product optimality remains
open; the stated theta subclass is no longer an open target.

## F1. What the two source pages ask

The [fence-area page](https://erich-friedman.github.io/packing/fence/) maximizes
the sum of enclosed areas, each at most1, using n unit fences without
crossings and with every fence endpoint on another fence. The
[product page](https://erich-friedman.github.io/packing/product/) instead
maximizes the product of polygon areas; edges meet only at endpoints and
there is no area-at-most1 condition. The pages describe best known
configurations. Their gallery values are not treated as general optimality
certificates.

For either problem, split all permitted T-junctions into vertices and take
the resulting embedded graph. We assume distinct segments have disjoint
interiors: coincident overlapping fences are not multiple graph edges.
The graph has minimum degree at least2. A connected pseudoforest has at
most one cycle; together these facts force every component to be a simple
polygon cycle. Indeed E<=V and minimum degree2 imply E=V and every degree2.
There can be no attached trees or degree-three shared-boundary vertices.
The splitting issue matters: original fence count is not generally the
same as the split graph's edge count, but in this pseudoforest case there
are no proper T-junctions and each cycle edge is an original unit fence.

## F2. A sharp12-fence pseudoforest obstruction for total area

Associate to each polygon cycle the face immediately inside it, subtracting
the interiors of its immediate nested child cycles if present. This pairs
the cycles with the bounded faces even when components are nested.
For a triangle of unit sides that face has area at most sqrt3/4<3/4.
For a cycle with m>=4 unit sides its associated face obeys the field-area
cap, hence has area<=1<=m/4. Consequently every
pseudoforest configuration satisfies

    total enclosed area <= n/4.                          (1)

The argument bounds the actual faces, not the uncapped gross area inside
an outer polygon; it therefore includes nested polygon components.
At n=12, three disjoint unit squares attain3, so the restricted maximum is
exactly3. A2-by2 grid of unit squares uses12 edges and encloses area4, each
face of area1. Thus attaining4 requires leaving the pseudoforest class.
The grid has9 vertices,12 edges and cycle rank4. This agrees with the
gallery's elementary12-fence construction without claiming global
optimality from the picture alone.

More generally, a k-by-l square grid has2kl+k+l edges and area kl; whenever
2kl>k+l it exceeds the pseudoforest bound at that same edge budget. This is
an infinite exact family of shared-boundary gains.

## F3. A sharp12-edge pseudoforest obstruction for product area

A simple equilateral m-gon has area at most

    a_m=(m/4) cot(pi/m),                                 (2)

attained by the regular polygon. One self-contained proof is discrete
Fourier decomposition of its complex vertices z_j: signed area is
(1/2) Im sum conjugate(z_j) z_(j+1), while the sum of squared edge lengths
is m. The ratio on Fourier mode k is(1/4)cot(pi k/m), whose maximum is
(1/4)cot(pi/m). Summing gives(2), and the first mode gives equality.

For one12-cycle the maximum is3(2+sqrt3). For k>=2 polygon components with
total12 edges, equation(2), cot x<1/x and AM-GM bound the product by

    (36/(pi k^2))^k.

There are at most four components. For k=2 this is<9 since pi>3; for k=3
it is<64/27, and for k=4 it is<(3/4)^4. All are below the regular12-gon
area, which is>9. Nested components only reduce the areas of bounded
regions relative to their constituent polygon areas. Therefore the exact
pseudoforest maximum at12 edges is

    3(2+sqrt3)=11.19615... .                              (3)

Now construct two reflected cyclic hexagons with side lengths2,1,1,1,1,1,
sharing the length2 side, which is split into two unit edges. The two
remaining boundaries each have five unit edges, for12 edges in all.
This is the theta graph with path lengths2,5,5: V=11,E=12,beta=2.

Let r=sqrt13, c=(r-1)/4, s=sqrt((1+r)/8), and h=(r-1)s/3. Start at(-1,0)
on the circle with centre(0,h). Rotate the radius clockwise five times by
the angle with cosine c and sine s. Each chord has length1 and the final
point is(1,0). Reflect this upper polygon in the horizontal axis for the
second face. Exact shoelace calculation gives the product

    ((89+26sqrt13)/16)=11.42152... .                      (4)

This independently realizes the formula credited to Erich Friedman at
n=12 on the product page. It exceeds(3) by0.2253683999..., so any design
achieving this value requires a multicyclic component. We do not claim
global optimality for(4). Its individual face areas exceed1, so it is
**not** a feasible configuration for the other, capped-area problem.

The exact coordinates, unit-distance checks and shoelace identity are in
verification.json. The numerical drawing is illustrative only; the proof
uses the radical identities.

## F4. Why pseudoforests, knots and the path carrier cannot be identified

The recovered functional-graph theorem says that a single-valued evolution,
after forgetting directions and collapsing duplicate edges, gives a
pseudoforest. Its simple cycles are vertex-disjoint. In particular the
two directed edges of a PRIMEGAME2-cycle collapse to one simple edge.
The separate numerator/denominator entry multigraph retains two parallel
edges for11/13 and13/11 and does have cycle rank1. These two graph models
must not be confused.

Enclosed-face counting uses the embedded graph's cycle rank E-V+c;
endpoint-path attachment means use a signed linear-forest kernel. Adding
shared faces changes the permitted incidence class, as(1)--(4) demonstrate.
No change of a scalar weight on a pseudoforest can recreate the missing
theta graph while preserving its embedding and shared edges.

The older knot work supplies a different rigorous continuation example:
THM-2176, `gordian-continuation-profile-and-interaction-cocycle`, retains
connected-sum contexts and distinguishes a relative chirality bit in its
specified mirror ambiguity. THM-2242,
`tournament-complement-transport-and-knot-kernel-green-rigidity`, makes
its full metric continuation kernel faithful, while warning that an
exterior histogram loses the labelled core kernel. These are inherited
results, not new knot invariants here.

A geometric fence needs a planar rotation system and lengths; a knot
diagram additionally needs over/under crossing data; a path carrier needs
its selected endpoint incidence and allowed operation language. The shared
principle is to retain the data tested by the next gluing. An isomorphism
between the resulting mathematical objects has not been supplied.
