# From the four-state microcosm to factorization, arithmetic and plane interfaces

2026-10-11. **PROVED:** the elementary structural and interface statements
below, with the all-order type bijection using the cited tournament/even-graph
count theorem. **CONDITIONAL:** the explicitly marked ABC radical limit.
**FINITE-EXACT:** the declared bounded implementations. **OPEN:** a rank-free
natural atom correspondence and the plane six/seven decision. No arithmetic
or colouring consequence is inferred merely from equinumerosity.

## Inheritance, portfolio and live board

Anchor: six-versus-seven plane colouring. Niche: an explicit tournament/even
type bijection. Wildcard: ABC difference families and golden factorizations.
The closest proved mechanisms are the
[four-state SCC codec](tournament_recursive_four_state_20261004.md),
[THM-4531, twisted equinumerosity and natural-map obstructions](../../01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md),
and [the prior channel/boundary integration](plane_cycle_boundary_integration_20261011b.md).
The hostile is the same four-state count with different allowed operations.
The corrected near miss is using a quotient after it has discarded a port,
norm sign or conditional hypothesis. The least-used sidecar is irreducible
factorization with its order, primitive scale and actual geometric coordinates.

The board is **marked two-bit section / ordered SCC word / connected graph
atoms / primitive arithmetic packet / unit-rhombus relation / geometric carrier**.
The session applies the existing META-PATTERNS cards "Type every analogy and
every implication" and "Find the hidden second coordinate in a nearly true
theorem". No new universal meta-pattern is promoted from this one session.

## 1. The annotated four-arc object is recovered exactly

The partial outdegrees {0,1,1,2}, with the specified missing-pair incidence,
force a unique frame

\[
 a\to b\to c\to d,\qquad a\to d,
 \qquad\text{missing }ac,bd.
\]

Its two free directions give the four tournament types once each. The
[complete recognition packet](partial_four_tournament_20261011c.md) also
classifies the other oriented squares. In particular, opposite source and
sink have the same degree list but yield T,T,S,S, so the missing-pair
incidence is an essential part of the information.

The previous K4 cycle-basis hostile shares the square as an actual binary
cycle-space element. It does not share the orientation predicate: inserting
the undirected edge bd creates the same monochromatic triangle in either
direction. The orientation frame is a valid two-bit comparison gadget;
unit-edge presence and colour equality require different interfaces.

## 2. A specific reversible tournament/even-graph construction

Here "even" means **automorphism-even**, as in the cited counting theorem,
not even degree at every vertex. A tournament decomposes uniquely into an
ordered word of strong tournament types. An automorphism-even graph
decomposes into a multiset of connected automorphism-even types.

The [new type codec](tournament_even_bijection_20261011c.md) bridges those
operations by the unique nonincreasing Lyndon factorization of the SCC word.
Weighted Lyndon factors and connected even-graph types have equal numbers
at each vertex weight: the count theorem and uniqueness of Euler-product
exponents prove this. Pair those **irreducible** lists by explicit canonical
rank, then replace the factors by their graph atoms and take disjoint union.
The inverse reads components, recovers their factors, sorts them and rejoins.

This is an all-order, terminating recipe on isomorphism types; the shipped
implementation caps its exhaustive atlas at order six. Canonical ranking is
still used on atoms. The construction is not a new proof of the count theorem,
not rank-free, and not relabelling-equivariant on labelled vertices.

At order four the resulting explicit inverse table is:

| Free reversal bits (x,y) | Tournament | Automorphism-even graph |
|---|---|---|
| 00 | Transitive | Four isolated vertices |
| 10 | Directed triangle over a sink | P3 plus an isolated vertex |
| 01 | Source over a directed triangle | K1,3 |
| 11 | Strong tournament | C4 |

The target graph edge parity is y+xy over F2. The mixed term is unavoidable:
these four graph types have edge counts0,2,3,4, so exactly one has odd parity.
On an affine Boolean square a linear parity observable has0,2,or4 ones.
Therefore no choice of labelled representatives turns this correspondence
into two fixed edge-toggle operations. The graph updates must retain context.

### The small interaction recurs with unbounded size

Let a denote the one-vertex strong type and b the directed triangle, a<b.
For every r>=1 the word a^r b is Lyndon: any proper suffix beginning with a
has fewer leading copies before its b, and a suffix beginning with b is
larger immediately. Thus its image under the type codec is connected.
In contrast b a^r factors as b|a|...|a, so

\[
 \Phi(a^r b)\text{ is connected},\qquad
 \Phi(b a^r)=P_3\sqcup rK_1.
\]

The order-four example is r=1; the same factorization mechanism persists at
every r. In the expanded graph, passing from Phi(a^r) disjoint-union Phi(b)
to Phi(a^r b) must add at least r edges, whatever vertex labelling or
additional deletions are used, because r+1 components must become one.
This is not a time lower bound for symbolic or compressed representations.
It is a precise reason a local operation before factorization can have an
unbounded response after the change of representation.

## 3. Prime microcosms have exact boundaries, not a finite prediction rule

The earlier twin-centre work already explains the exceptional chain3--5--7:
one of any three consecutive odd candidates is divisible by3. It is the
only chain of overlapping twin-prime pairs. The pair2,3 is uniquely adjacent
as integers. All other twin centres are divisible by6; that condition is
necessary, not sufficient.

The [ABC/difference-family packet](abc_difference_families_20261011c.md)
goes beyond that observation in three concrete directions:

* Every finite prime wheel has an explicit infinite arithmetic progression
  of survivors whose two neighbours are both composite. For primes3,5,7
  the progression9450+30030t has forced proper neighbour factors11 and13.
* For odd n>=9, with a=(n+3)/2 and b=(n-3)/2,
  F_a L_b=F_n+2(-1)^b and F_b L_a=F_n-2(-1)^b. Both products are nontrivial.
  The prime Fibonacci values with a twin-prime neighbour are exactly3,5,13.
* **Assuming ABC**, primitive normalization of the signed Pell relation
  L_n^2-5F_n^2=4(-1)^n gives
  log(rad(F_n L_n))/log(F_n L_n)->1. It does not imply prime values or a
  density of simultaneous primes. Twin-centre additive packets themselves
  have radical larger than height, so are not exceptional small-radical cases.

The historical Fibonacci/Lucas exploration had a wrong general norm sign,
a parity-based gcd rule and an unsupported independence inference. The
producer and saved prose are repaired with explicit historical banners;
the old numerical output is retained and is not claimed to be a fresh run.
The correct gcd is2 exactly when3 divides the index. Current conditional
proofs retain this primitive normalization.

## 4. Actual progress and obstructions for the plane

The [rhombus packet](unit_rhombus_interfaces_20261011c.md) proves that a
four-vertex unit C4 has diagonal squared lengths summing to4. Thus both
missing pairs cannot simultaneously be unit edges. All four orientation
completions exist abstractly, but only three geometric edge-presence states
exist. The required geometric map fails already in the microcosm.

A unit diamond forces its nonadjacent apices to have equal colours at
palette three. At palette six they can differ: its exact colour count is
k(k-1)(k-2)^2. Thus reusing the familiar three-colour equality gadget at six
colours is another false extension. This identifies what a stronger gadget
would actually have to replace.

There is a positive restricted family: two marked isometric unit triangular
lattices give an explicit six-colouring, with colour3j+(a-b mod3) on the
selected carrier address(j,a,b). The exact Moser spindle is a checked example.
Every edge-sharing connected triangle patch lies in one such lattice and
has exactly six labelled three-colourings. The coordinate and palette
hypotheses remain part of the certificate.

This cannot be a universal unit-scale covering theorem: three pairwise
closer-than-one points already require three unit carriers, despite their
unit graph being empty. Separately, the arithmetic packet supplies a complete
two-colouring of every rational-coordinate unit graph, so putting prime
filters on Pythagorean rational parameters cannot create a six/seven witness.

The surviving target is therefore exact irrational, multiscale geometry
with a complete palette interface. The prior adaptive cycle/boundary solver
can authenticate a candidate or a finite refutation; the new partial-square
and carrier rules can supply tested subroutines. Neither the bijection nor
ABC transports a plane-colouring predicate without an additional proved map.

## Incoming work: keep the operation as well as the state

Incoming commit85ba11616 was integrated before publication. Its
[Vitali recovery](vitali_atom_recovery_20261011.md) sharpens the same boundary:
a target survives an observation precisely when it is constant on every
observation fibre. In particular, finite tournament score fibres need not
be binary affine spaces. Our four-state section really does have two marked
bits, but its graph image does not inherit their fixed-toggle operation.
Keeping four distinguishable states and keeping their algebra are different
requirements. The codec preserves the former and explicitly computes the
context-dependent replacement for the latter.

The subsequent incoming commit60ea42e66 supplies a contrasting four-element
object: [M5 of the affine-source packet](mobius_affine_source_return_20261011.md)
classifies the rational symplectic stabilizer of its fixed generic invariant,
using a cited Legendre section-group calculation. Its four maps really form
a Klein group under composition. This is inherited work, not a new audit of
that input. The concrete comparison is at the level of operations: its group
law and the marked tournament bit flips are explicit, whereas cardinality
alone would not give either law or make graph edge toggles implement it.

The incoming adjoint repair also distinguishes conjugacy from equality;
self-converse tournaments need not have Hermitian adjacency. This agrees
with retaining the actual converse action on the marked ports, rather than
identifying a state with its transformed presentation. Canonical finite
catalogue choice here is explicit and does not invoke a nonmeasurable
transversal or transfer a dense-translation obstruction to a lattice.

The [composed Collatz weight calculation](collatz_composed_weight_floor_20261011.md)
offers a separate reusable mechanism: retaining affine offsets puts all
composite anchors in one fixed interval, so distortion need not accumulate
with composition depth. Its native guards and independently grounded ROOT
certificates remain additional coordinates. It does not establish universal
coverage, and neither the graph codec nor the arithmetic factorization
supplies those missing certificates. The next useful transport must preserve
the actual terminal or colouring predicate, not just the number of states.

## Verification

Every companion declares its finite universe, positive controls, hostile
examples and reproduction commands. Independent agents audited the proofs
and APIs; normal and optimized runs agree with saved outputs. The all-order
mathematical statements are supported by the included proofs, rather than
extrapolated from those finite checks. No plane chromatic decision follows.

| Packet | Explicit exact checks |
|---|---:|
| [Partial four-tournament](partial_four_tournament_20261011c.out) | 5,194 |
| [Tournament/even type codec](tournament_even_bijection_20261011c.out) | 84,251 |
| [ABC and difference families](abc_difference_families_20261011c.out) | 5,685 |
| [Unit rhombi and carrier certificates](unit_rhombus_interfaces_20261011c.out) | 243 |
| Total | 95,373 |
