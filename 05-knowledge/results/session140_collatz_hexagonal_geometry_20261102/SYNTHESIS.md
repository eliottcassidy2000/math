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
First-return and exact-clock addressing results are reserved for a separately
reviewed continuation. No new global method-card promotion is made here.
