# Full colour-boundary tomography by actual unit-distance ears

**PROVED in the scopes below; FINITE-EXACT placement and counting controls
in `receipt.json`.** No new chromatic lower bound is claimed. In this packet,
"a six-chromatic witness" requires an actual unit-distance graph with a
certificate of no five-colouring. "No six-colouring" would be a different,
stronger seven-chromatic claim. The graphs constructed here are three-colourable.

## Inheritance and the connection being tested

The [full boundary interface](../../boundary_palette_compression_20261011.md)
already retains colour-equality partitions, exact gluing and witness lifts.
The [six-colour staircase](../../hn_metagraph_information_20261026/six_color/RESULT.md)
separates prescribed frames from ordinary colourings; the
[metric packet](../../rank_geometry_metagraph_20261026/unit_distance/RESULT.md)
retains cross-unit contacts and has a positive two-terminal diamond transfer.
The [five-interior theorem](../../paired_budget_selectors_matrix_20261028/FIVE_INTERIORS.md)
excludes a particular small paired-palette selector. These are inherited,
not conclusions newly inferred from the present experiment.

Closest geometric mechanism: [Session140's generic circle avoidance](../../session140_collatz_hexagonal_geometry_20261102/hadwiger_nelson_refinement/PROPER_UNIT_DISTANCE_REFINEMENT.md).
Canonical hostile: unrecorded cross-unit edges change a colour tensor.
Corrected near miss: complete numerical observations do not make a forcing
gadget; the [rooted-cost discussion](../../rooted_graph_costs_20261020/README.md)
already distinguishes these and does not claim q>=3 star completeness.
Least-used sidecar here: the geometric availability of *joint products* of
tests, rather than merely their individual pair marginals.

Board: full equality partition; actual planar placement; two private probe
vertices; simultaneous probe products; positive versus signed operations.
Anchor = an exact geometric counterpart of endpoint-response tomography;
niche = the smallest probe's metric boundary; wildcard = testing whether
complete observation can yield a chromatic obstruction.

## G1. A two-private-vertex path measures the complete two-terminal table

Let G be a finite simple graph with two distinct named vertices u,v. Write
Z_G(q) for its number of proper colourings with q labelled colours, q>=2, and

    E_G(u,v;q)=#{proper q-colourings c of G with c(u)=c(v)}.

Adjoin two distinct private vertices x,y and exactly the three new edges
u-x, x-y, y-v. Call this a three-edge **ear**. Retain any old edge u-v.
For an old colouring c, literal counting of the two new colours gives

    number of extensions = K_q - [c(u)=c(v)],
    K_q = q^2 - 3q + 3.                                      (1)

If the endpoint colours agree, the count is(q-1)(q-2)=K_q-1; if they differ,
it is(q-1)^2-(q-2)=K_q. Thus

    Z_(G+ear)(q)=K_q Z_G(q)-E_G(u,v;q).                       (2)

For any specified equal endpoint colour the full boundary weight is E_G/q;
for any specified unequal pair it is(Z_G-E_G)/(q(q-1)). Global palette
permutations justify these divisions. Consequently these two scalar counts
recover the entire two-terminal q-colour tensor, including its zero entries.
This is a direct counting identity, not a new chromatic-polynomial theorem.

## G2. Exact planar realization and its sharp metric boundary

Let P be any finite set of distinct planar points and u,v two of them.
If 0<|u-v|<3, two new points x,y can be chosen so that the complete induced
unit-distance graph gains **exactly** the ear in G1. There are no collisions
and no other new unit contacts to P.

**Proof.** Vary x on the unit circle about u within an open arc where
0<|x-v|<2. Such an arc exists exactly in the indicated distance range:
the distance from v to that circle has infimum ||u-v|-1|<2.
The two unit circles centred at x and v intersect transversely, giving
two local continuous choices for y. If |u-v|=1, one choice is identically
y=u; discard it and use the other. Otherwise neither choice is identically u.
Shrink the open arc if necessary to avoid the finitely many degeneracies.

There are only finitely many forbidden x positions. A collision with an
old point, or an unintended unit contact from x, excludes finitely many
points on the circle about u. For y to hit an old point w, x must also be
a unit neighbour of w: the circles about u and w meet in at most two points,
except the already discarded case w=u. For y to have an unwanted unit
contact with old w!=v, y belongs to the at most two intersections of the
circles about v and w. Each such y!=u permits at most two positions of x.
Hence only finitely many x values are excluded, and an open arc remains.
The required x-y contact itself precludes x=y. QED.

This can be repeated for finitely many ears, including repeated ears on a
pair, treating all earlier private vertices as old points. Therefore a
whole finite family can be placed simultaneously with no cross-ear edges.
It is not necessary to assume that independently drawn fragments happen
to avoid one another.

**Boundary.** If |u-v|>3, the triangle inequality forbids such a unit ear.
If |u-v|=3, equality forces the unique straight positions
x=u+(v-u)/3 and y=u+2(v-u)/3. These can work when the old set is just u,v,
but may fail injectivity or gain extra contacts. For u=(0,0),v=(3,0) and
an old point w=(1,0), x=w is forced and no private realization exists.
Thus <3 is a uniform guarantee with the old placement fixed, and the
distance3 boundary cannot silently be included. An arbitrary abstract
neighbor-set probe is not licensed by this much smaller construction.

## G3. Bell-many realizable contexts recover the full boundary, not marginals

Let B=(b_1,...,b_m) be an ordered boundary in a fixed injective finite planar
placement, and let G be its full induced unit-distance graph, with every
pair of boundary points at distance<3.
Keep q>=2 and **all** boundary equality partitions. For a partition beta
with k<=q blocks, let w(beta) be the number of extensions of any particular
labelled colour assignment of that pattern. Put

    mu(beta)=(q)_k w(beta),

the total number of colourings having that pattern. For k>q set mu(beta)=0.
The current graph is not assumed colourable, and no boundary state is removed
in advance because it seems incompatible with a proposed palette pairing.

For each set partition pi choose the canonical forest F_pi: in each block,
join its smallest index to every other index in that block. Form the actual
unit-distance context C_pi by adding one ear for each edge of F_pi.
G2 realizes all these ears simultaneously in a common point set using at
most2*binom(m,2) new vertices. Each individual query uses2(m-|pi|) vertices.
Subsets of that common set realize each C_pi as an induced graph.

By(1), the observed total is

    Z_pi = sum_beta mu(beta) product_(ij in F_pi)(K_q-[beta_i=beta_j]). (3)

There are exactly Bell(m) such queries, including the empty-forest query
Z_G. These suffice to recover every mu(beta), by the following explicit
integer inverse; no rank assumption about realized graphs is needed.

Let M_pi count base colourings whose equality partition coarsens pi,
so M_pi=sum_(beta coarsens pi)mu(beta). For S subset F_pi, let pi(S) be
the component partition of the subforest S. Because every block uses its
smallest vertex as star centre, F_(pi(S))=S. Expanding delta=K_q-(K_q-delta)
on each forest edge gives

    M_pi = sum_(S subset F_pi) (-1)^|S| K_q^(|F_pi|-|S|) Z_(pi(S)). (4)

Finally process partitions from fewer to more blocks and use

    mu(pi)=M_pi-sum_(beta strictly coarsens pi)mu(beta).            (5)

Divide by(q)_|pi| when |pi|<=q to recover the exact labelled boundary
weights. Equations(4) and(5) give the promised complete tensor and support.
For impossible >q-block states the recovered masses must be zero.

The joint contexts in(3) are essential. Merely observing each individual
pair ear supplies pairwise equality marginals; that is a different and
generally smaller measurement family. No minimal-query theorem for actual
unit-distance graphs is claimed. Bell(m) is a sufficient explicit protocol.

## G4. A sharp obstruction to turning these measurements into forcing

For every q>=3,

    K_q-1=(q-1)(q-2)>0.

Therefore every proper colouring of G extends across each ear, regardless
of its endpoint pattern. For any finite collection of ears with disjoint
private vertices and exactly their specified edges, the complete boundary
support is **unchanged**. In particular the assembly is q-colourable if
and only if G is q-colourable. This holds for the simultaneous induced
unit-distance construction in G2, not just for an abstract drawing.

At q=6 the two weights are20 and21; at q=5 they are12 and13. Exact arithmetic
can extract a zero boundary weight from signed combinations of these
measurements, but a physical graph union multiplies nonnegative extension
counts. Inclusion-exclusion subtraction in(4) is not a graph gadget.
Hence this entire ear-only construction language cannot turn a five-colourable
base into a six-chromatic witness, or a six-colourable base into a graph
with no six-colouring. Complete observability does not supply forcing.

The palette threshold is sharp. At q=2 the equal-endpoint factor is zero.
Take an induced two-edge unit path with endpoints(-1,0),(1,0) and centre(0,0).
It has two proper two-colourings, both with equal endpoint colours. G2 adds
an induced three-edge ear, giving C5 and zero two-colourings. The supplied
exact coordinates realize this hostile; the graph has30 three-colourings.

This result is scoped to ears with no additional private interactions.
Extra edges, identified private vertices, or gadgets with a genuinely zero
boundary entry require a new analysis and remain possible research lanes.
The inherited paired-selector no-go is not being generalized to all gadgets.

## G5. One exact37-point experiment recovers the six-colour wheel table

Use six boundary points at the vertices of a regular unit hexagon and its
centre as the base graph. Its induced graph is the wheel on seven vertices.
Coordinates are in Q(sqrt3). Its full six-colour boundary weight is

    w(beta)=6-|beta|, if consecutive hexagon ports have different colours;
           =0 otherwise.                                        (6)

There are203 equality partitions. Exactly40 have positive weight, accounting
for14910 of the46656 labelled boundary assignments. There are24600 full
proper six-colourings. The all-distinct orbit is a genuine zero: all six
colours around the centre leave none for it. The graph is nevertheless
three-colourable. A zero orbit is not an empty full relation.

The exact construction in `coordinates.json` adds one ear for each of the15
boundary pairs. It has37 distinct points,57 unit edges, and609 other pairs
certified nonunit. The203 forest contexts are literal induced subgraphs of
this same placement. Their totals range from24600 to97219384920, and the
integer inverse(4)–(5) recovers every entry of(6), including all163 zeros.

For the three possible endpoint distance types, one-ear counts are:

| Ports | Squared distance | Total six-colourings | Recovered equal-endpoint count |
|---|---:|---:|---:|
|0,1|1|516600|0|
|0,2|3|510360|6240|
|0,3|4|512280|4320|

The full37-point graph comes with a literal three-colouring, checked on
every unit edge. Thus this experiment tests a full boundary recovery on
actual geometry while supplying its own decisive hostile to an inferred
six-chromatic obstruction.

## Exactness, reproduction and remaining target

Run `python3 verify.py` here; `python3 -O verify.py` must give the same receipt.
`--write` is an explicit producer mode for the coordinate and table artifacts;
ordinary verification compares them without overwriting either.

Each new x is given by the rational unit-circle parameter
((1-t^2)/(1+t^2),2t/(1+t^2)) about u. With Delta=v-x and d2=Delta dot Delta,
the new y is

    (x+v)/2 + sign*sqrt(1/d2-1/4)*perpendicular(Delta).

All coefficients are in Q(sqrt3). The checker verifies the required unit
identities exactly in this quadratic field. Every collision and nonedge is
excluded by rational interval arithmetic for these explicit real radicals,
using integer square roots and160 binary fractional bits. No floating-point
distance tolerance, probabilistic test, or coincidence assumption is used.
The largest search used13 candidate placements; this is a finite construction
receipt, while G2's unbounded existence proof is separate.

Counting has two implementations: partition-weighted equations and a literal
loop over all46656 boundary assignments with every pair of probe colours.
Literal graph-colouring backtracking independently checks the three one-ear
totals and the C5 hostile. All203 reconstructed masses match the entire
original tensor. Frozen inherited source hashes are in `SOURCE_PINS.json`.

Source-to-target map: actual boundary colour tensor -> realized ear-context
counts -> equality moments -> full tensor. Preserved: every labelled boundary
partition weight and its support. Lost by individual marginals: higher joint
patterns, restored by actual simultaneous forest probes. Lost by counting
alone: the ability to enforce a zero or attach arbitrary geometric neighbors.
The sidecars are the complete point configuration, cross-unit incidence,
and the distinction between signed reconstruction and positive graph gluing.

The next forcing target must leave this positive-kernel language, while
retaining a realizable metric and the full five- or six-colour boundary
relation appropriate to the claimed chromatic bound. No UNSAT certificate
for such a new target is supplied by this experiment.
