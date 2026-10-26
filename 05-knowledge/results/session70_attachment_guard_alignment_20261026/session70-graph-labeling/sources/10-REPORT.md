# Square-sum islands, shared ancestry, and an authenticated three-vertex move

**PROVED:** the founder-reachability tail starts exactly at71; a pinned
three-vertex recoloring repairs an actual attachment although every allowed
one-vertex recoloring fails; reflection/residue and square-scaling transports.
**FINITE-EXACT:** all founder-copy colorings through35 and an explicitly
defined family of95 arithmetic coloring candidates. No claim about the
chromatic number of the infinite square-sum graph, universal Hamiltonicity,
Collatz, or an unrestricted tournament wiring algorithm is made.

## 1. Inheritance and the meaning of the proposed islands

Let `Q_n` have vertices `1,...,n` and edge `a--b` iff `a != b` and `a+b` is a
square. Edges are undirected. When studying ancestry, orient an edge from its
larger endpoint to its smaller endpoint; that extra birth order is retained.

Closest proved mechanism: the reflection `x -> k²-x` and its directed arrival
neighbors. Canonical hostile: a graph-degree-two vertex can be a Hamiltonian
path endpoint. Corrected near miss: three startup components do not recur as
three connected components at every scale. Least-used sidecar: actual
descending paths to each founder, rather than a single inherited color.

Relevant current sources:

- [THM-4510, reflection orbits and Pythagorean zigzags](../../../../01-canon/theorems/THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md).
- [THM-3357, Berggren three-branch Walsh level collapse and parent circuit](../../../../01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md).
- [Session40 square-sum packet](../../session40_marked_ports_square_sum_20261023/square-sum/README.md): actual endpoints, shared path incidences, edge-history colors, and boundary-conditioned birth formulas.
- [Summand closure and square filter](../../collatz_mod6_20260922_w6_summand_closure_square_filter.md): founders `{1,2,4}`, connectivity from14, and the distinction from the strict-summand closure holes `{1,4,6}`.

The live board is: founder sets / exact descent addresses / reflection phases /
proper versus ancestry colors / pinned coupled moves / labelled Pythagorean
chains. The anchor is shared attachment compatibility; the niche is square-sum
ancestry; the wildcard is the three-branch Pythagorean transport.

At12 the three actual components are the paths

```
2--7--9                   (3 vertices),
12--4--5--11               (4 vertices),
10--6--3--1--8             (5 vertices).
```

They are founded by2,4,1, respectively. These sizes really are `3,4,5`.
At13 two components merge;14 merges the remaining two. The full graph stays
connected thereafter.

The user's later three sets are exactly the components of the induced graph
`Q_18 minus {13,14,15}`:

```
{12,4,5,11}, {2,18,7,9,16}, {10,6,3,1,8,17}.
```

This is a valid deleted-vertex view, with component sizes4,5,6. Its sidecar is
the three omitted vertices and all their incidences. In particular15 has only
founder1 ancestry: removing15 cuts a cycle but does not remove a mixed-founder
arrival. The three sets are not components of the full `Q_18`, and are not a
proper vertex coloring because consecutive vertices within each set are
adjacent.

## 2. The exact least tail where all three ancestries coexist

Write

```
D(n) = {k²-n : n < k² < 2n},
A(1)={1}, A(2)={2}, A(4)={4},
A(n) = union_(u in D(n)) A(u) otherwise.
```

These are exactly the founders reachable by strictly descending square-sum
paths. No quotient by undirected connectivity is taken. The recurrence uses
sets, not a chosen parent, and `A(n)` need not increase with n.

**Theorem.** The first `n` with `A(n)={1,2,4}` is22. The least integer `N` such
that every `n>=N` has this property is **71**.

**Proof.** The elementary recurrence gives the finite base table through99,
stored in [receipt.json](receipt.json). A separately written descending DFS
checks every entry through120. At22 the union of the ancestries of3 and14 is
all three; all earlier unions omit a founder. At70 the predecessors are
`11,30,51`, with ancestry sets `{4}`, `{1}`, `{1,4}`, respectively. Thus70
does not reach2. Every entry71 through99 reaches all three. For example:

```
71 -> 10 -> 6 -> 3 -> 1,
71 -> 29 -> 7 -> 2,
71 -> 29 -> 20 -> 5 -> 4.
```

For the infinite step, let `n>=100` and `k=ceil(sqrt(n+71))`. Since

```
(n-73)² - 4(n+71) = n²-150n+5045 > 0,
```

and the right side is45 at100 and increases thereafter,
`2 sqrt(n+71) < n-73`. Hence

```
n+71 <= k² < n+72+2 sqrt(n+71) < 2n-1.
```

The actual older neighbor `u=k²-n` satisfies `71<=u<n`. Strong induction
gives `A(u)={1,2,4}`, and therefore the same holds at n. The failure at70
proves minimality. QED.

Thus exclusive founder allegiance cannot be intrinsic descending reachability
beyond70: every later vertex has all three possibilities. This does not forbid
a *chosen* three-color ancestry forest. It says that such a choice is extra
information, and must remain compatible at shared predecessors.

## 3. Unary possibility sets lose the first coupling at22

Fix colors0,1,2 on founders1,2,4. A **founder-copy coloring** satisfies

```
c(n) in {c(u):u in D(n)}
```

at each nonfounder. Equivalently choose one older parent per nonfounder, and
copy that parent's color. This is a coloring of the birth ancestry; it is
usually not a proper coloring of `Q_n`. The set of possible colors of one
vertex, marginalized over all coherent choices, is exactly `A(n)`.
The parent-to-color map is not injective:15 may choose parent1 or10, both
colored0. Thus the counts below count distinct colorings, not parent forests.

**FINITE-EXACT:** the first prefix for which the product of all unary
possibility counts overcounts coherent colorings is22. There are32 coherent
colorings, versus48 unary combinations. The obstruction is explicit:

```
c(14)=1 and c(22)=2
```

is impossible, even though each requirement is separately realizable.
Vertex22 can only copy3, whose color is0, or14. If14 has color1, color2 is
unavailable at22. This is the same shared-owner failure as independently
choosing compatible attachment charts at two incidences and forgetting that
the middle owner is shared.

Complete counts from orders20 through35 are

```
8,16,32,64,128,192,336,588,882,1323,1323,2646,5292,13230,26460,66780.
```

The complete enumeration preserves the whole prefix, not just the color
histogram. At35 there are66,780 colorings. A global product of unary sets is
already invalid at22; the all-three tail theorem does not make later choices
independent.

For lossless descent information use the three-component path polynomial

```
F_seed(seed)=1, F_other(seed)=0,
F_a(n)=sum_(u in D(n)) x_(n,u) F_a(u).
```

Here edge variables retain actual ordered endpoints. Its monomials are the
strictly descending paths to a. Setting all variables to1 retains their
multiplicities; replacing coefficients by support gives `A(n)` and loses
both multiplicity and address. This is the direct analogue of the Fibonacci
zero-slot warning: discarding a length-zero or already-summed coordinate does
not license identifying its two adjacent owners. No particular Fibonacci
recurrence is asserted for these square-sum polynomials.

## 4. An actual pinned coupled move repairs the first greedy three-color failure

Now use **proper** colors, so adjacent vertices have different colors.
Ascending first-fit chooses the smallest available color from `{0,1,2}`.
It colors `Q_34`, then fails at35: the older neighbors

```
1,14,29 have colors 0,1,2.
```

Failure of this algorithm is not failure of three-colorability.

**Proposition.** Freeze every old vertex except `K={1,15,24}` at its
first-fit color. Every allowed recoloring of a single vertex of K is invalid.
Nevertheless a simultaneous swap of colors0 and2 on K is valid, preserves
every frozen boundary owner, and extends to a proper coloring of `Q_35` by
giving35 color0.

**Proof.** The induced graph on K is the star with edges `1--15,1--24`.
Its colors are `(0,2,2)`. Its actual outside boundary is

```
{3,8,10,12,21,25,34},
```

all colored1. Every vertex of K has an internal neighbor of the other color
and an outside neighbor of color1, so neither alternative single-site color
is legal. The simultaneous0/2 swap retains every internal inequality and all
outside inequalities. The colors of the new vertex's neighbors become
`2,1,2`, so color0 is available. QED.

This is a component switch with authenticated boundary pins, not an analogy
based only on the number3. Without these pins, recoloring14 alone from1 to2
also repairs the arrival; the frozen outside data is essential to the claim
that simultaneous movement is necessary. The state space of K with that
boundary has precisely two states, `(0,2,2)` and `(2,0,0)`, although the move
touches three vertices. The measured number of color0 vertices in K changes
from1 to2; the one-site feasible-move graph nevertheless has no edges. This
is an explicit observation that a test using only legal one-site moves would
incorrectly declare invariant. A longer connected bipartite block with all boundary
vertices fixed at the third color has the same two-state versus large-support
phenomenon whenever each internal vertex has both types of neighbor.

Connection contract: source = proper coloring with pinned boundary; target =
compatible attachment fiber; map = restrict to actual block and boundary
owners; preserved predicate = every edge inequality and every pin; lost data
under one-site projection = the phase correlation across the block; restoring
sidecar = one common component-flip bit; decisive test = all one-site moves
fail while the displayed simultaneous move passes.

## 5. Two and three colors that carry exact arithmetic information

### Parity is a valid edge cocycle, not a proper vertex coloring

For an edge with sum `k²`,

```
(a mod2) XOR (b mod2) = k mod2.
```

Thus the edge color "odd or even square root" is the binary coboundary of
vertex parity. Every closed walk has an even number of odd-root edges. This
is an all-size conservation law. It allows even-root edges to join same-parity
vertices; it does not assert a bipartition. Indeed the inherited five-cycle
`1--3--6--10--15--1` already excludes proper two-coloring at15.

### Modulo three gives two edge symbols acting on three actual states

Put `r_t(x)=t-x (mod3)`. A square sum is0 or1 modulo3, so the two edge symbols
act by

```
r_0 = (1 2), fixing0;    r_1 = (0 1), fixing2.
```

They generate `S_3`. Following an edge word and an initial vertex residue
recovers the terminal residue exactly. Two consecutive reflections translate
by the difference of their square targets. The initial residue and edge order
are necessary sidecars; symbol totals alone do not determine their product.
Fixed points of these transpositions show why the residue labels are not
proper colors.

This is the precise sine/cosine connection. On a formal Fourier mode pair
`cos(theta*x),sin(theta*x)`, reflection about `t/2` acts by

```
[ cos(theta*t)  sin(theta*t) ]
[ sin(theta*t) -cos(theta*t) ].
```

Composing reflections yields the rotation/translation phase of `t2-t1`.
This calculation lives on the full arithmetic carrier; on `Q_n` each
reflection is a partial matching and must retain its actual domain, endpoints,
and no-loop condition. A mode calculation does not restore a deleted port.

### Doubling toggles the target carrier; scaling by four returns

Let `R_n` join positive integers whose sum is twice a square. Multiplying all
vertex labels by2 identifies `Q_n` with the induced even-vertex subgraph of
`R_(2n)`, and identifies `R_n` with the induced even-vertex subgraph of
`Q_(2n)`. Multiplication by4 is therefore an induced embedding of `Q_n` in
`Q_(4n)`. The parity/target carrier is essential: doubling alone is not an
endomorphism of the square-sum graph (`1+3=4`, but `2+6=8`).

### Cheap hostile probes for procedural coloring searches

No periodic finite vertex coloring can properly color the infinite square-sum
graph. If M is a period, `M²` and `3M²` have the same color and their sum is
`(2M)²`. Thus any coloring based on finitely many modular registers has an
explicit hostile after taking their common period.

For a finite procedural screen, test every nonzero coefficient vector in

```
c(n)=a*n+b*v_2(n)+c*popcount(n)+d*floor(sqrt(n)) (mod q), q=2 or3.
```

There are15 binary and80 ternary coefficient vectors. All fail proper
coloring by order10 and11, respectively. The receipt gives the first bad edge
for every vector. Constants are immaterial to properness and were omitted;
distinct vectors need not define distinct functions. This is a scoped
95-candidate failure, not a classification of arithmetic colorings.

## 6. The genuine Berggren connection, and the coordinate it loses

THM-4510 already supplies a stronger bridge than matching three color names:
every primitive triple `s²+t²=u²`, with `s<t`, gives a labelled Hamiltonian
chain on `[1,t²-1]` using the three target squares `s²,t²,u²`. Its edge
reflections retain the arithmetic. In THM-3357's branch convention the root
and its three children give:

| branch | triple (odd leg, even leg, hypotenuse) | chain order | ends |
|---|---|---:|---|
|root|(3,4,5)|15|8,9|
|L|(5,12,13)|143|25,72|
|M|(21,20,29)|440|200,400|
|R|(15,8,17)|224|32,64|

The checker builds each entire three-square union and verifies that it is a
connected path, not merely that its size matches a formula. Keeping the three
square labels recovers the triple and its Euclid parameters, hence its
Berggren address. The naked path shape loses this information.

Even chain order does not recover the first branch. The primitive triples
`(29,420,421)` and `(341,420,541)` both give a path on176399 vertices. Their
parameter pairs are `(14,15)` and `(10,21)`, with branch words `L^13` and
`R^4 M`. Thus their unlabelled paths and vertex counts coincide although the
first branches differ. This is an explicit obstruction to recovering a
Berggren allegiance from the unlabeled square-sum chain alone.

The startup islands have the lengths3,4,5, but no branch-preserving map from
all founder-copy colorings to the Berggren tree has been supplied. The exact
route is through labelled Pythagorean target matchings, with content, leg
ordering, and branch address retained. The exceptional Hamiltonian orders
`15,16,17,23` remain inherited finite-startup behavior; the numbers7 and21
acquire no additional identification from this packet.

## 7. Reproduction and productive continuation

Run from the repository root:

```
python 05-knowledge/results/attachment_lattices_square_sum_20261023/square_sum/check.py
python -O 05-knowledge/results/attachment_lattices_square_sum_20261023/square_sum/check.py
```

The frozen [receipt](receipt.json) contains all1..100 ancestry masks, the
finite induction base, explicit paths, founder-copy counts, all95 coloring
hostiles, both pinned colorings, actual boundary vertices, and the Berggren
branch witnesses. The standard and optimized runs execute the same checks;
no correctness condition relies on a Python `assert`. Independent small
adjacency and DFS implementations check the generated relations.

The next useful object is a boundary-conditioned transfer polynomial for
ancestry/color choices with actual owner variables, not another list of unary
colors. Contract a region only after identifying its shared owners; retain its
joint relation and multiplicity. The three-vertex pinned example supplies a
positive test for coupled moves, while the22 hostile rejects a compiler that
multiplies marginal supports. The71 theorem states exactly when a coarse
three-founder support ceases to distinguish new vertices.
