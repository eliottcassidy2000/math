# What the geometry preserves, and what it must pay

## Direct progress: the hexagonal picture has an exact survivor

For ordinary `C(n)=n/2` on even n and `3n+1` on odd n, put
`S={n>0:3 does not divide n}`. Choose k by residues modulo nine:

```
n mod9: 1 2 4 5 7 8
k:      2 1 2 3 4 1
F0(n)=32n; F1(n)=(16n-2^(4-k))/3.
```

**PROVED.** Both children are distinct, larger than n, remain in S, and
map to n under exactly five ordinary steps. Equal-clock fibres of distinct
parents are disjoint, yielding `2^r` endpoints. Five is the sharp uniform
clock for two unit-surviving predecessors: at four steps, root7 has the
complete fibre `{18,112}`, only one of which is a 3-unit. The shortcut map
has clock four and the corresponding three-step hostile `{18,56}`.
Sharpness concerns this fixed-clock requirement, not optimal growth rate.

Every component contains a 3-unit at finite forward distance. Therefore its
undirected balls have exponential growth, without assuming convergence or
classifying possible unknown components. If selected depth-r predecessor
walks map injectively to the honeycomb or triangular lattice, some ordinary
edge has endpoint distance at least `(2^(r/2)-1)/(10r)`. Any fixed finite
fibre multiplicity still fails. Fixed-dimensional polynomial-growth targets
have the same capacity obstruction. Uniform Euclidean separation gives a
packing version; arbitrarily small scales evade that hypothesis.

**PROVED with a CITED input.** The ordinary graph is a countable locally
finite subcubic pseudoforest. Its finite-level strip drawing has no vertex
accumulation point and no component has infinitely many nested cycles.
Kühn's square-grid minor theorem applies with its actual proviso. Contracting
the honeycomb vertical perfect matching gives the square grid; a median in
each subcubic minor branch set turns the model into a subdivision. Directed
paths preserve every transition, but their lengths necessarily grow without
a uniform bound. The result is a subgraph, not induced or spanning, and not
a canonical effective layout. Branch markers and clocks must be retained.

## Niche progress: unit distance is not lattice locality

**PROVED.** Every countable locally finite simple pseudoforest has an induced
unit-distance realization: start each cyclic component with a regular unit
polygon, then place leaves on unit circles while avoiding finitely many
collision or unwanted-distance locations. Fair interleaving reaches all
vertices. Crossings and arbitrarily close vertices are allowed; a placement
need not give an algorithm for discovering Collatz's unknown components.

Functional pseudoforests are at most three-colorable, and the ordinary
Collatz triangle forces exactly three. Hence a faithful representation of
this graph alone cannot force four, six, or seven colors. The ambient plane
unit-distance graph and a chosen induced subgraph are different objects.

**PROVED.** For specified oriented edge-displacement labels, exact cycle
sums are the obstruction to a position lift. Modular cycle closure alone
is insufficient: of all30 rooted triangles in the seven-label quotient,
12 close exactly and18 only modulo the ideal. Injective positions and
unwanted unit distances remain separate checks after a position lift.

## Wildcard recovery: seven historical HN mechanisms, not seven new claims

The source-pinned recovery distinguishes logical compactness, corrected
all-lattice coloring, arithmetic residue/spectral obstructions, named-field
configuration exclusion, moment/elliptic correspondences, perfect-code
quotients, and the metric hexagonal-cell lift. Corrections and original
status travel with each item. In particular, a quotient requiring seven
colors is not seven unit-distance points in the plane; the classical plane
seven-color construction needs explicit cell-diameter and separation bounds.
Inherited descent and modularity are not silently promoted to fresh proofs.

## Deeper lesson and honest remaining work

The useful common principle is **preserve the predicate, not the picture**.
Degree-three compatibility does not preserve metric capacity. A graph
quotient does not preserve exact displacement. A subdivision preserves
incidence and reachability only with original vertices and clocks identified.
The literal arithmetic hexagon is doubling on `(Z/9Z)^*`; division by three
consumes precision, so mod9 is not a closed dynamical model by itself.
Abstract graph rigidity, even when valid, does not prescribe unique geometry.

No single forward Collatz orbit visits every positive integer. A geometric
space-filling traversal would require an extra rule; the outdegree-one
condition does not produce one. The independent audit validates the scoped
claims, not convergence, a new HN bound, or equivalence of major conjectures.
The first-return and exact-clock addressing continuation is independently
accepted below. No new global method-card promotion is made here.


## Checkpoint2: three faithful models, with different costs

The new synthesis makes the retained predicate explicit rather than choosing
one picture as universally correct:

| Model | Retained exactly | Necessary cost / boundary |
|---|---|---|
| Ring-order site conjugacy | All lattice sites and original one-step arithmetic dynamics | Nonlocal jumps; jump lengths grow without bound; not a single spanning trajectory |
| Honeycomb subdivision | Directed incidence and marked-state reachability | Auxiliary routing sites, branch markers and variable clocks; unbounded path lengths |
| Proper induced Euclidean unit distance | Every simple edge/nonedge at distance one; compact local finiteness | Arbitrarily close nonadjacent pairs at infinity; straight-edge crossings allowed |

The first-return theorem is the exact clock repair. Let J be injective and
let each tau(x) be the first positive return to J(X), with return endpoint
J(F(x)). Successive return times are sums of these positive clocks. Thus
marked reachability and first hits are equivalent, least cycle periods are
clock sums, and eventual periodicity is preserved in both directions. The
extra condition that every target state enters J(X) is needed before ruling
out all off-marker cycles. Merely prescribing `R^tau J=JF` is insufficient:
a two-state swap with proposed clock 2 realizes the endpoints of an identity
map but creates a spurious marked hit at time 1.

First-return exactness does **not** require different routes to have disjoint
interiors. Shared interiors are permitted when their successor agrees and
no extra marked state is crossed. A canonical disjoint suspension and an
actual native expansion can induce the same marked-state dynamics while
having different shared-tail incidence. Pairwise meeting times, distances
and vertex counts are different consumers and may need that incidence.
In particular, subdividing original one-step Collatz edges must not be
confused with independently inflating accelerated odd edges: native halving
tails can be shared. This is consumer-relative exactness, not a defect in
the first-return theorem.

The proper Euclidean refinement uses nested horizontal intervals and upward
unit arcs. Every child rises by more than `sqrt(15)/4`; finite depth levels
and separated component strips imply local finiteness in compact sets.
Finite forbidden circle intersections preserve induced adjacency. This is
not a Collatz descent coordinate: the chosen root need not lie forward of
every vertex in a hypothetical nonconvergent component. Exponential growth
still excludes uniform separation. Properness means the increasingly close
pairs can move out to infinity rather than accumulating at a finite point.

## Curvature is a capacity control, not a Collatz solution

Place the abstract depth-k binary-tree vertex j at
`((j+1/2)2^(-k),2^(-k))` in the hyperbolic upper half-plane. The standard
metric formula gives constant edge cosh-distance 21/16. Exact parity of the
horizontal numerator shows that only tree edges attain the sharp minimum
L. The geodesic circle equations keep edges inside disjoint height bands
and horizontal intervals, proving noncrossing and compact local finiteness.
Consequently exponential branching, uniform separation and bounded routes
can coexist in a space topologically equivalent to the plane; Euclidean
lattice capacity, not plane topology, caused the earlier obstruction.

The indispensable hostile is retained: depth-k neighbors across the root
split have intrinsic tree distance 2k but constant hyperbolic distance
`arcosh(3/2)`. This placement is therefore not quasi-isometric. Its boundary
crowding is on y=0, outside hyperbolic space. It does not route all ordinary
Collatz edges or construct a hyperbolic hexagonal tessellation.

## Consumer-relative sufficiency across the recovered history

An observable factors through a quotient exactly when it is constant on
fibres. An operation also needs fibre-constant legality and output. For
linear boundary responses, the erased kernel must annihilate the span of
**legal** consumers; hypothetical delta pins do not establish realizable
geometric gadgets. Joint boundary assignments must be multiplied before
summation, not replaced by proper marginals. Arithmetic repair depends on
the coefficient group: the triangle Laplacian defect is repairable over R
but not Z. A complete derivative jet can retain a residue norm floor when
a scalar value vanishes, but the actual coefficient congruences and bounded
parameter cube remain hypotheses. These are precise successful mechanisms
and hostiles, not shared-word analogies or a recertification of all history.


## Concurrent native-port bridge, with a rank boundary repaired

The [narrow independent joint-port audit](independent/joint_ports/AUDIT.md)
recovers the concurrent J1–J4 identities, not its whole packet. Odd children
of the same parent share halving tails: merging those exact integer states
recovers the native graph, whereas independent macro-edge subdivisions do
not. The full exponents and odd lowest common ancestor recover native first
meetings and their two asynchronous clocks. Sources 13 and 53 first meet at 40
at times (1,3); sources 853 and 53 first meet at 160 at times (5,1), even though
the changed exponents 3 and 9 agree modulo 6. A periodic exponent summary
therefore loses this fixed-port response while a first-return clock can
still preserve the odd-state map.

**Explicit boundary correction:** the leaf-path Gram rank `2^h` holds for
`h>=1`. At `h=0`, the sole leaf is the stopped root, its edge-indicator is
zero, and the rank is **0**, not 1. The auditor checked this exception and
the positive-depth private-edge mechanism. The upstream packet and its
PNG-dependent seal remain byte-frozen; the included SOURCE.md is labelled
historical evidence, not an endorsement of its unqualified h=0 sentence.
The correction and exact source identities travel in this addendum.
