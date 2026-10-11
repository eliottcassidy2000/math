# Pair bits, pants seams, and a multiscale six-versus-seven target

2026-10-11. **PROVED / FINITE-EXACT:** the scoped interfaces, identities,
packing bound and enumerations linked below. **CITED:** high-girth preservation.
**CONDITIONAL:** six/seven consequences under the accepted plane lower-bound
premise H6. **OPEN:** deciding between six and seven, and universal Collatz.
No algebraic six- or seven-chromatic unit-distance witness is supplied here.

## The structural count, and the distinction it exposes

The literal expression 10^3+10^2-10 equals1090. The nearby identity is

\[
 1099=10^3+10^2-1=1+2+8+64+1024
      =\sum_{n=1}^{5}2^{\binom n2}.
\]

It counts all labelled simple graphs on one through five vertices. Palette
threshold three is an additional question about those graphs, not a filter
implicit in this sum. The meaningful triangular number is binom(n,2), the
number of independent unordered-pair bits. A chosen vertex order converts
the same bitstrings to tournament orientations, but changes the relabelling
action and does not preserve the graph-colouring predicate.

The other tournament/even-graph identity concerns **unlabelled tournaments**
and **automorphism-even graphs**, not arbitrary labelled graphs or even-degree
graphs. The current inherited mechanism is
[THM-4531, twisted Mallows--Sloane and equivariant-map obstructions](../../01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md).
Equal finite cardinalities always give some bijection. What can fail is a
bijection respecting the specified relabelling or switching structure.
[The typed transfer packet](tournament_even_transfer_20261011b.md) preserves
those qualifications and the earlier cycle-basis correction.

## The connection board

Anchor: the plane six/seven alternative. Niche: tournament characters and
composite-modulus flows. Wildcard: pants decompositions and triangular seam
counts. The useful common object is a **marked constraint system with a lift**,
not the recurrence of a numeral. This table records the actual transfers.

| Source -> target | Map and preserved predicate | Lost information / retained sidecar | Cheapest hostile |
|---|---|---|---|
| Labelled graph -> pair-bit cube | Record edges on binom(n,2) pairs; preserves exact graph with labels | Vertex order is needed for orientation; graph colouring is not a tournament invariant | Same bits can encode an empty graph and a transitive orientation under the chosen convention |
| Tournament orbit counts -> even graphs | Signed fixed-point/Fourier identity; preserves averaged orbit count | A count gives no vertexwise equivariant decoder | Empty graph has full symmetric stabilizer; no tournament does |
| Six-colouring -> ternary labels plus bits | c=2t+b; proper iff same-label edges admit alternating bits | Labels alone lose within-label obstructions; retain parity paths | Monochromatic C5 passes every triangle test |
| Graph fragment -> boundary parity | Keep named components and relative bit parity for fixed t; preserves extendibility | Retain paths, internal obstruction state, and full colour lift | An even path and an odd path have incompatible common endpoints |
| Pants -> cuff bits | Restrict a binary character; preserves local cuff holonomy | Genus g leaves g independent seam bits; retain marked dual graph and reference transport | Genus-two theta gluing with all cuff bits zero and a nonzero seam class |
| High-girth graph -> unit realization | Preserve graph adjacency and chromatic number through cited construction | Minimum point spacing can shrink; retain exact geometry and scale | Fixed spacing plus sufficiently high girth forces three colours |

The last two rows explain the same kind of failure in different categories:
valid local restrictions need not determine the global object. They are not
an identification of surface cohomology with a plane colouring invariant.
Unit-distance graphs drawn in the plane may have crossing edges; they are
not automatically planar graphs. Planar duality and the four-colour theorem
therefore cannot be imported without another proved map.

## A different-looking exact statement deciding six or seven

For any finite graph G, form the hypergraph whose hyperedges are its odd-cycle
vertex sets. Then

\[
 G\text{ admits six colours}
 \quad\Longleftrightarrow\quad
 \text{its odd-cycle hypergraph admits three labels with no monochromatic edge}.
\]

This follows by partitioning the six colours into three pairs. Each pair
must induce a bipartite graph, and any such partition lifts to six colours.
The [channel packet](plane_odd_cycle_channels_20261011b.md) proves the general
2r version, derives its exact counting/Fourier formula and implements a
terminating finite learner. Each failure produces an actual odd cycle, which
forbids precisely the corresponding monochromatic label assignments.

A cycle basis of the **original** graph is not enough for that learner:
K4 with labels (1,0,0,0) passes every triangle through vertex0 but has a
monochromatic triangle123. The residual same-label graph must be analysed
anew, or its complete constraints retained. This is the nonlinear loss
behind a tempting reuse of a linear cycle-space compression.

Under H6 and the usual seven-colour upper bound, the plane is six-chromatic
exactly when every finite algebraic unit-distance graph passes this ternary
test; it is seven-chromatic exactly when one fails. The geometry and the
complete refutation of all ternary assignments remain indispensable.

The [signed boundary compiler](plane_six_reframes_parity_20261011b.md)
exports exact relations and path witnesses for a fixed ternary assignment.
Its gluing failures expand back into real odd cycles. It therefore supplies
a compositional oracle for the channel learner, rather than a scalar score
that can declare success without an actual colouring lift.

## A genuine geometric restriction on every high-girth route

Tomon's cited high-girth preservation theorem leads to the conditional
reformulation: proving a six-colour upper bound for **all** unit-distance
graphs above even one fixed girth threshold proves it for all unit-distance
graphs. Conversely, a seven-chromatic witness would yield such witnesses
at every prescribed girth.

The [crowding packet](plane_girth_crowding_20261011b.md) proves that a
k-critical unit graph of girth>=2r+2 has minimum point separation at most

\[
 \frac{2r}{\sqrt{1+(k-1)((k-2)^r-1)/(k-3)}-1}\qquad(k\ge4).
\]

The mechanism is exponential breadth-first tree growth versus quadratic disk
area. A related proved upper bound is stronger at fixed resolution:
separation>=1/1000 and girth>=64 force 2-degeneracy, hence three colours.
The six/seven alternative therefore lives beyond every one fixed spacing
when the high-girth construction is used. This gives a concrete missing
coordinate rather than an appeal to vaguely difficult global structure.

## Pants and triangular numbers, with an explicit decoder

In a fixed three-band seam system on a pair of pants, let x,y,z be arc counts
between the three pairs of distinct cuffs. Their endpoint counts are

\[
 (a,b,c)=(x+y,x+z,y+z).
\]

The determinant is -2. Nonnegative integer recovery is possible exactly when
the triangle inequalities hold and a+b+c is even. At x+y+z=N, the counts
are in bijection with two cut positions among N+2 places, hence number
binom(N+2,2). This is a second concrete triangular-number mechanism.
It describes the fixed cross-cuff bands, not all arcs or all twists.

The [pants packet](pants_boundary_holonomy_20261011b.md) then proves, for a
marked genus-g closed surface, that compatible cuff bits have dimension g
while complete binary characters have dimension2g. Each cuff assignment
has a torsor of 2^g seam gluings. Tree normalization retains the missing
g bits and gives a decoder. This is the useful information-storage lesson:
compression can be exact, but only relative to the observable and the
marking whose reconstruction it promises.

## What incoming work changed, and the next precise targets

The incoming [source-frame carrier work](source_frames_eit_20261011.md)
preserves the actual source tests rather than their count; the incoming
[Collatz weight transport](collatz_bridge_weight_transport_20261011.md)
preserves source weights and affine offsets under refinement. These support
the same research move: identify and carry the data that future composition
still observes. Their theorems do not themselves prove plane colourability,
and the colouring interface does not supply missing ROOT certificates.

The recovered knight-torus comparison,
[THM-4552, six as a two-place product and seven as a Paley torus](../../01-canon/theorems/THM-4552-exceptional-knight-tori-six-is-a-two-place-product-seven-is-a-paley-torus.md),
is an instructive failed shortcut: its numeral six indexes a bipartite
finite quotient, not a six-chromatic real unit graph. Categorical, Cartesian
and lexicographic products preserve different relations.

The next targets are now precise:

1. **Completed prototype:** [adaptive learning through boundary fragments](plane_cycle_boundary_integration_20261011b.md)
   now retains a full lift or a complete finite refutation. The integrated
   and direct oracles agree on all 75 graphs through order4, all 6,594 label
   cases, and the geometric pentagon/spindle controls. Different witnessed
   cycles can legitimately give different search traces.
2. Extend geometric assembly through progressively close clusters while
   retaining exact unit-edge equations and named port transport. Nearness
   alone never permits identifying two vertices.
3. Seek a uniform six-colour gluing invariant across those scales, or an
   algebraic unit realization of a finite odd-cycle hypergraph that needs
   four weak labels. Every local rule must be checked against the long-cycle,
   original-basis and seam-twist hostiles above.

The results give exact proof obligations and reusable checked operations;
they do not establish either global alternative.

## Reproducibility close-out

All six companion scripts have explicit finite universes and hostile controls;
normal and optimized runs agree with their saved outputs. The proofs and
implementation contracts received independent agent review. The checks below
are explicit assertions within those declared universes, not independent
theorems or a plane-wide enumeration.

| Packet | Exact checks |
|---|---:|
| [Odd-cycle channels](plane_odd_cycle_channels_20261011b.out) | 8,382 |
| [Signed boundary compiler](plane_six_reframes_parity_20261011b.out) | 2,040,717 |
| [Integrated adaptive learner](plane_cycle_boundary_integration_20261011b.out) | 13,811 |
| [Tournament and flow transfer](tournament_even_transfer_20261011b.out) | 19,067 |
| [Pants and seam holonomy](pants_boundary_holonomy_20261011b.out) | 14,857 |
| [Girth and crowding](plane_girth_crowding_20261011b.out) | 2,309 |
| Total | 2,099,143 |
