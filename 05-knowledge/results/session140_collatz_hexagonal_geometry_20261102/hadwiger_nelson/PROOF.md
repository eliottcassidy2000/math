# Session 140: functional dynamics, hexagonal carriers, and exact lift certificates

**Status.** Elementary theorems proved below; finite controls are supplementary.
These are new-to-this-session connections, not literature-wide priority claims.
They neither prove Collatz convergence nor improve the chromatic number of the plane.
The recovered Hadwiger–Nelson results and their correction history are in
`HISTORICAL_RECOVERY.md`; exact Git source identities are in `SOURCE_PINS.json`.

## 1. Objects that must not be conflated

Use standard Collatz `C(n)=n/2` for even positive n and `3n+1` for odd n.
Its directed graph has one outgoing edge per vertex. Its incoming vertices at m
are `2m`, and additionally `(m-1)/3` exactly when `m=4 mod6`.
The shortcut map `T(n)=n/2` or `(3n+1)/2`, and the odd map
`U(n)=(3n+1)/2^{v_2(3n+1)}`, are different directed graphs.
In particular C has the triangle `1→4→2→1`; U has a loop at 1.
Throughout, an *underlying simple graph* forgets orientations, removes loops,
and collapses parallel/opposite edges. A loop is not a unit-distance edge.

The honeycomb graph joins **vertices** of regular hexagons and is bipartite,
with degree three. The hexagon-**center** neighbor graph is the triangular
lattice T, degree six and chromatic number three. The full unit-distance graph
of R² joins every pair at distance one: its straight drawn edges may cross,
and it is not the class of planar graphs. Arbitrary subdivisions preserve
neither edge length nor chromatic number nor one-step dynamical time.

## 2. Functional graphs cannot themselves force four colors

**Theorem 1 (functional pseudoforest).** If a countable directed graph has
outdegree at most one, its underlying simple graph is a pseudoforest:
each component contains at most one finite cycle. It is three-colorable.
Consequently the underlying standard Collatz graph has chromatic number
**exactly three**, without assuming the Collatz conjecture.

**Proof.** A finite vertex set spans at most as many undirected edges as
vertices: charge each retained edge to one of its directed sources; a source
is charged at most once. A connected finite graph containing two distinct
cycles has more edges than vertices. If a component contained two cycles,
join them by a finite path and apply the count to their finite union, a
contradiction. Every finite subgraph therefore has a vertex of degree at
most two, and recursive deletion/reinsertion gives a three-coloring.
For a countable graph, the finitely branching tree of partial three-colorings
has arbitrarily deep nodes. Choosing the least child with arbitrarily deep
descendants gives a coloring, exactly the choice-free argument recovered
in THM-4622. The triangle on 1,4,2 gives Collatz's matching lower bound. ∎

This does not say that an arbitrary geometric placement has no additional unit edges.
Those additional edges are extra constraints, not consequences of outdegree
one. Conversely, selecting one outgoing edge at each point of a difficult
unit-distance graph usually discards the chromatic obstruction.

## 3. Some unit-distance placement is cheap; a lattice placement is not

**Theorem 2 (induced Euclidean realization).** Every countable locally finite
simple pseudoforest has an injective placement in R² with distance one
**if and only if** two vertices are adjacent. Crossings are allowed.

**Proof.** In each component with a cycle, first put its cycle on a regular
polygon of side one. It has no other unit chords. In an acyclic component
start at one root. The remaining vertices form rooted trees. Add a tree
vertex on the unit circle about its already placed parent, avoiding all
old vertices and all unit circles about other old vertices. At each finite
stage the forbidden set on the parent circle is finite: two distinct
unit circles intersect in at most two points. Thus a permitted point exists.
For a new component choose its initial root/polygon by a generic translation,
avoiding finitely many collisions and unit distances to earlier vertices.
Local finiteness permits a countable fair interleaving of the rooted-tree
construction over components, with finitely many vertices at each stage.
Every pair's adjacency or nonadjacency is secured when its later endpoint
is placed. ∎

Therefore a Euclidean unit-edge drawing of the abstract Collatz graph does
not require solving Collatz. It may be induced but it does not, by this
construction, give a planar straight-line drawing or a lattice drawing.
This is an **existence proof**, not an effective algorithm deciding the
components or cycles of the infinite Collatz graph. Such structural knowledge
is used to order the construction; it is not supplied by the construction.
The obstruction below forces arbitrarily close nonadjacent points whenever
all edge lengths are one in such a realization.

**Nonuniqueness.** Fix the labeled equilateral placement of 1,4,2. The
vertex 8 is attached to 4, and may be chosen at any of continuum many points
on its unit circle after finitely many exclusions. Each choice extends by
the same existence construction. Any isometry between two resulting labeled
placements must fix the three noncollinear anchors and hence be the identity.
The placements with different positions of 8 are therefore noncongruent.
The dynamical degree rules do not select a unique geometric curve.

## 4. Independent slower check of the lattice-growth obstruction

The primary Session140 Collatz lane proves the **stronger equal-five-step**
fiber bound `|C^{-5r}(m)|≥2^r` for every positive m with 3∤m. This section
retains an independently derived odd-preimage proof as a redundant control,
not as a separate primary advance or the sharp session constant.

Let `S={m≥1 : m odd and 3∤m}`. For `m∈S`, define

```
A(m)={a∈{1,2,3,4,5,6} : 2^a*m mod9 is 4 or 7},
P(m)={(2^a*m-1)/3 : a∈A(m)}.
```

**Lemma 3.** `|A(m)|=|P(m)|=2`; every n in P(m) belongs to S and U(n)=m.

**Proof.** Multiplication by m permutes the six units modulo nine, and two
has order six modulo nine. Exactly two exponents hit 4 or 7. Divisibility
by three holds, the quotient is odd, and residues 4,7 make the quotient
nonzero modulo three. Moreover `3n+1=2^a*m` with m odd, so its exact 2-adic
valuation is a. The two exponents give different positive integers. ∎

Set `S_0={1}` and `S_{k+1}=⋃_{m∈S_k}P(m)`.

**Lemma 4.** `|S_k|=2^k`, and each n in S_k reaches 1 in at most 7k
iterations of C. Thus the undirected Collatz ball of radius 7k around 1
contains at least `2^k` vertices.

**Proof.** Predecessor sets of distinct parents are disjoint because U is a
function. Hence each level doubles. A U step with valuation a is one odd
C step followed by a halving steps, at most seven C steps. Concatenate.
The sets at different depths need not be disjoint: 1∈P(1), reflecting the
known C cycle. Only disjointness of sibling-parent preimages at the same
depth is used. ∎

**Theorem 5 (no uniformly bounded lattice routing).** There is no injective
map φ from the standard Collatz graph into a triangular or honeycomb
lattice for which a constant L bounds the lattice-path distance between
φ(n) and φ(C(n)) for every n. The same conclusion holds for any target
graph whose balls have polynomial volume growth. It holds already on
the known component containing 1.

**Proof.** The triangular ball of integer radius R has
`1+3R(R+1)` vertices. A honeycomb lattice is a subgraph of a suitably aligned
triangular lattice, so the same upper bound suffices. For every n∈S_k,
concatenating the routes gives target distance at most `7Lk` from φ(1).
Injectivity would imply

```
2^k ≤ 1+3(7Lk)(7Lk+1)  for every k,
```

impossible since an exponential eventually exceeds every polynomial.
For a polynomial-growth target replace the right side by its ball bound. ∎

L can be taken integral by rounding up. No noncrossing assumption, edge-
disjoint routing condition, or surjectivity onto the lattice was used.
The result therefore also prohibits any stronger space-filling lattice
model with uniformly bounded subdivision. The honeycomb's absence of a
triangle is a smaller separate obstruction for unsplit C edges.
For the shortcut T and odd U, replace 7k by 6k and k respectively, after
ignoring the odd map's loop for the underlying graph metric.

**Quantitative finite version.** For the finite union of the chosen inverse
paths through level k, let L_k be its largest routed-edge distance. Then

```
L_k ≥ [sqrt(12*2^k-3)-3]/(42k).
```

Indeed the least possible target radius satisfies
`R≥(sqrt(12*2^k-3)-3)/6`, and `R≤7kL_k`.
This forces stretch of order `2^{k/2}/k` along a growing inverse family.
For L=1, already k=16 gives `65536 > 1+3*112*113 = 37969`.
The k=16 calculation is an explicit obstruction threshold for this proof,
not a claim to the smallest forbidden Collatz subgraph.

**Euclidean separation version.** If φ takes values in R², all Collatz
edges have Euclidean length at most L, and distinct vertices are separated
by at least δ>0, disjoint radius-δ/2 disks about φ(S_k) fit into a radius
`7Lk+δ/2` disk. Hence
`2^k≤(14Lk/δ+1)^2`, again impossible. In particular Theorem 2's unit-distance
realizations cannot have uniform positive vertex separation.

This is a growth/metric obstruction, not a convergence theorem. It uses
explicit inverse families already known to reach 1 and says nothing about
whether every positive integer belongs to that component.

## 5. Exact cycle-lift certificates: what modular coloring forgets

Let A be an additive abelian group, `q:A→B` a surjective homomorphism,
and S a symmetric set of allowable nonzero steps with `0∉q(S)`.
For a connected finite graph H, suppose each vertex v has a distinct label
`y_v∈B`, and each oriented edge uv has a selected lift `s_uv∈S` with
`q(s_uv)=y_v-y_u`, `s_vu=-s_uv`.

**Theorem 6 (cycle displacement criterion).** These edge choices admit
positions `x_v∈A` with `q(x_v)=y_v` and `x_v-x_u=s_uv` if and only if the
signed sum of s along every graph cycle is zero in A. It is enough to
check a fundamental cycle basis. The resulting positions are injective.

**Proof.** Necessity is telescoping. Conversely choose a spanning tree and
one lift of its root label; define positions by path sums along the tree.
Every nontree edge has the correct difference exactly when its fundamental
cycle sum vanishes. The q-labels telescope correctly, and distinct q-labels
ensure distinct positions. ∎

For number-field unit-distance problems take A to be the additive group
of the relevant local ring and S its exact norm-one vectors. Reduction
supplies a graph homomorphism; a finite-field closed polygon only has
sum zero **after reduction**. A Euclidean unit-distance lift needs zero
cycle displacement in characteristic zero, simultaneously for every cycle.
Kloosterman moments count the former constraint and are not certificates
of the latter. Fixed edge lifts may fail even when all local incidences
are correct; choosing compatible exact lifts is the remaining problem.

### A seven-color quotient with only three-color source

Use `ω=e^{2πi/3}`, `Λ=Z[ω]`, and unit steps
`±1,±ω,±(1+ω)`. The triangular lattice T is three-colorable by
`a+b mod3` on `a+bω`. The residue map

```
c(a+bω)=a+2b mod7
```

has kernel `(2-ω)Λ`; its six unit steps map bijectively to the six nonzero
residues, so the quotient graph is K7. The quotient triangle
`0→1→2→0` has uniquely prescribed unit lifts `1,1,-ω`, whose sum is
`2-ω≠0`. It is closed modulo the ideal but not in the lattice.
This explicitly exhibits the lost cycle displacement.

Thus `χ(T)=3` and `χ(T/(2-ω))=7` are compatible. A proper target coloring
pulls back, so `χ(source)≤χ(target)`; a lower bound on an arbitrary quotient
need not transfer back. The recovered Session61 graph square T² also has
chromatic number seven, but its edges have lengths 1,√3,2. A separate
Voronoi-cell diameter/separation argument is what gives the actual
seven-coloring of every point of the Euclidean plane.

## 6. The common preservation principle, stated as a test

For dynamics, a map `q:X→Y` factors f into an exact map `g:Y→Y` precisely
when `q(x)=q(x')` implies `q(f(x))=q(f(x'))`. This is necessary by applying
g and sufficient by defining `g(q(x))=q(f(x))`. Even an exact factor may
create periodic orbits: `f(n)=n+1` on Z factors modulo m to a cycle although
f has none. Branch information and a lifting condition are necessary to
recover the original trajectory.

For coloring, a graph map is enough to pull back upper bounds, but lower
bounds require a witness mapping the other way with its exact edges.
For arithmetic configurations, modular cycle sums must be lifted to exact
cycle sums. For lattice geometry, cardinality alone is irrelevant: an
injection can exist while bounded edge stretch cannot.

The reusable question is therefore not whether two models have the same
hexagons, seven labels, or binary branching, but which **predicate** the
map preserves, what information its fibers discard, and what finite or
infinite lifting obligation remains. This connects the corrected residue
prime, descent, lattice, and Collatz stories without equating their open
problems.
