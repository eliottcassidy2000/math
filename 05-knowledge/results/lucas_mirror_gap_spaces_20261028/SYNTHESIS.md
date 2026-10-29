# What these bridges preserve: oddness, dimension, and a missing matching

**PROVED elementary observations / INHERITED / CITED / OPEN**, separated below.
This is a source-and-operation atlas, not a claimed equivalence of major open
problems. Source baseline90e704a8f356; concurrent work incorporated through
b07ec51b20. No universal Collatz, LRC(14), planar Jacobian, or six-colour result.

## 1. The concept board and inheritance pass

| Concept | Closest proved mechanism | Canonical hostile / missing data |
|---|---|---|
| Golden words | THM-4592's all-period reader, Session74 proof | even-period zero fibre needs three source tags |
| Geometric symmetry | actual Lucas toggle graph; Session91 golden metric | arithmetic translation need not be a plane isometry |
| Incidence and tournaments | Paley blocks and the signed-even Walsh basis | point and block actions can differ |
| LRC gap states | THM-370's two-neighbour criterion | all-runner failure is not marked-runner failure |
| Supersingular primes | corrected level filtration1\|2\|6 | Baby Monster is a quotient of2.B, not a subgroup ofM |
| Dimension | induced cube realizations and retained square equations | unit edges do not retain cubical angles or all distances |

Anchor: give the eleven-state reader its actual moves and reflection.
Niche: identify the continuous gap-space behind the Lucas states.
Wildcard: determine why the same symmetry closure at29 states is much larger.
The least-used relevant data were the point-to-block matching and the
positive-excess coordinates of a gap, rather than another total count.

The literature check recovered Brown's exact pentagonal biplane reflection,
so that permutation and PSL(2,11) are credited as classical. The new assembled
tests are the twelve matching frames, their fifteen-arc observation ambiguity,
the explicit A29 hostile, and the threshold-space homotopy construction.
The inherited generic plane hypercube realization was also found in a stronger
repo theorem, so it is cited below rather than reproved as new.

## 2. Four genuinely different occurrences of oddness

**(a) Odd cyclic period.** The golden word reader is injective for odd n.
This comes from the conjugate-root/carry analysis; it does not say its state
set has odd cardinality. For example n=3 has L3=4 states. The signed word
count is six-periodic and gives a topological Euler characteristic, as proved
in [GAP_TOPOLOGY](GAP_TOPOLOGY.md).

**(b) Tournament automorphism order.** Every finite tournament has an
automorphism group of odd order. An even-order group would contain an
involution; a nonidentity involution swaps some vertex pair and would have
to reverse the arc between them. For a self-converse tournament, adjoining
all anti-automorphisms doubles the automorphism-group order, so Cauchy's
theorem supplies an anti-involution. It has at most one fixed vertex, since
an arc between two fixed vertices could not reverse. At odd order it has
exactly one fixed vertex.

The Lucas mirror fixes three of its eleven vertices. It therefore cannot
be either an automorphism or an anti-involution of any completed tournament
on those vertices. That explains a symmetry obstruction before enumerating
any completions; the biplane repair retains it by allowing different actions
on points and blocks.

**(c) A sign-even graph need not have odd automorphism order.** Lambda5 has
Aut=D5 of order10. Its unique degree-five vertex is fixed; the degree-three
vertices and degree-two vertices identify the two outer layers. Suppressing
degree-two vertices yields a five-cycle, proving the full automorphism group.
Every automorphism preserves the bit-addition orientation, so no edge orbit
can swap its endpoints. Its Royle endpoint-reversal sign is therefore trivial:
Lambda5 is a connected sign-even graph with15 edges.

Here “even” is **not** even degree at each vertex and is not an even number
of edges. The canonical change of basis between tournament invariants and
signed copies of sign-even graphs is inherited
[THM-4721, twisted Burnside and Walsh basis](../../../01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md),
with its linked corrections. Equinumerosity does not produce a symmetry-
preserving graph-to-tournament bijection on the same vertex set.

**(d) Self-converse classes are a trace, not half of a smaller enumeration.**
Converse acts on the tournament-invariant function space. In the class basis
its trace is SC_n. In the sign-even Walsh basis the eigenvalue of the X row
is(-1)^|E(X)|. Hence the inherited basis theorem gives

    SC_n = (# sign-even graph types with even edge count)
           - (# sign-even graph types with odd edge count).

This is a character identity, not a set bijection preserving individual
objects. MISTAKE-081 already refuted SC_n=A000568(n-1).

### A concrete cancellation certificate supplied by the new graph

For a tournament T on eleven labelled vertices, sum over distinct embedded
copies of Lambda5 the product of its15 edge signs, measured against bit
addition. This is well-defined because all ten graph automorphisms preserve
that reference orientation. There are11!/10=3,991,680 such copies.

If T is self-converse, an anti-automorphism pairs positive and negative
copies: all15 factors reverse, so the product changes sign. No copy can be
fixed by that pairing action, because its sign would equal its negative.
Even if the chosen anti-automorphism has higher order, it still bijects the
positive and negative sets; an anti-involution supplies a literal pairing.
Thus exactly1,995,840 copies have each sign.

For Paley11 the given Lucas embedding is positive, and x->-x supplies an
explicit negative partner. So a genuine positive local geometric embedding
can coexist with exactly zero global signed invariant. This is relevant to
the ROOT research as a sharp analogy: a positive summand alone cannot prove
a positive aggregate. No map to actual Collatz source weights is supplied.

## 3. The Monster15 split has its own established mechanism

The corrected inherited
[Monster-level note](../monster_levels_square_sum_7_21_20261022.md), including
MISTAKE-668, already identifies the user's4/5/6 partition:

| Prime support layer | Primes | Size |
|---|---|---:|
| Monster minus Baby Monster |29,41,59,71|4|
| Baby Monster minus Fi22 |17,19,23,31,47|5|
| Fi22 |2,3,5,7,11,13|6|

The modular-curve levels are1,2,6, with primes dividing the level included
separately. This is the established rationality/centralizer dictionary,
not a new partition discovered this session. Baby Monster's relevant
Monster centralizer is the double cover2.B. Mere divisibility already gives
some nesting of the genus lists; that nesting alone is not new moonshine.

There is an exact representation-theoretic comparison to converse. For any
finite set with involution tau, if d is its cardinality and f its fixed count,
the rational permutation module has trace f and anti-invariant dimension
(d-f)/2. For supersingular j-invariants Frobenius is such an involution;
absolute rationality means this anti-invariant space vanishes. At higher
levels Aricheta's theorem identifies the analogous defect with a genus
difference. See [Aricheta, Theorem1.1](https://www.emis.de/journals/SIGMA/2019/007/sigma19-007.pdf).
For tournament classes the involution is converse and f=SC_n. These are
applications of the same involution decomposition to **different modules**;
no identification of the underlying sets or their dimensions is known here.

The exact PSL(2,11) biplane route reaches the classical Mathieu construction;
Brown's paper describes the associated code/design connections. It does not
turn the fifteen Lucas edges into the fifteen Monster primes. A direct map
would need to transport the involution or incidence data, not just the count.

## 4. Fifteen objects can fail the cheapest symmetry test

There are15 Lucas edges and15 pairs of six labels. Let D5 act naturally on
five cyclic labels with a sixth fixed label. On K6 pairs its orbits have
sizes5,5,5; on Lucas edges the orbit sizes are5,10. Under a reflection the
fixed-edge counts are respectively3 and1. Hence no D5-equivariant bijection
exists between these two fifteen-object sets with those actions.

Restriction to C5 removes this obstruction: both split into three free
five-orbits, so an equivariant bijection can be chosen after choosing how
to pair and phase those orbits. The lost information is reflection and the
choice of orbit alignment. This is a concrete, small example of why an
equinumerosity may become usable only after its symmetry contract is stated.

Similarly10=2*5 is intrinsic here: the two nonempty weight layers are two
rotation orbits. The arithmetic identity10=7+3 does not pick out a preserved
partition in this construction. One can impose such a split, but it has not
been explained by the Lucas or biplane action. The Monster4/5/6 split instead
comes from a separately specified level filtration.

## 5. Cuboids, simplices and spheres depend on the retained metric

The stronger inherited
[cubelike geometry theorem](../ports_regular_actions_hn_20261025/geometry/README.md)
already proves: a Cayley graph of an elementary binary group with chosen
directions admits an injective plane unit-edge realization iff the directions
are independent. In particular every finite hypercube, and hence each induced
Lucas cube, has a faithful plane unit-distance realization. High cube rank
does not force high unit-distance embedding dimension. This fact is not new
here; it is also compatible with Kovač–Predojević's work on flexible hypercube
patterns in the plane.

**Elementary spherical rigidity.** Let all2^d additive vertices
sum b_i u_i lie on one sphere, with nonzero u_i in a real Euclidean space.
Write the squared-distance equation for the four vertices0,u_i,u_j,u_i+u_j
and subtract:2<u_i,u_j>=0. Hence the u_i are mutually orthogonal and the
ambient dimension is at least d. The sphere centre projects to half their
sum, giving radius at least sqrt(sum |u_i|²)/2, with equality for a centre
in their span. Conversely an orthogonal cuboid has that circumsphere.

Thus “all cube edges have unit length” and “the full additive cube lies on
one sphere” retain different information. A regular simplex with m vertices
likewise requires dimension m-1 when all pairwise distances are retained:
its centred Gram matrix has rank m-1. Dimension becomes a natural-number
invariant after specifying independent directions, affine span, simplicial
dimension, or homology degree. Cardinality by itself does not specify which
dimension is intended.

For independent supports, the coordinate cubical complex is a cone; the
simplex section of that cone is an independence complex; that link can be a
sphere. This is the exact cuboid/simplex/sphere relation used by G1. It does
not require regarding a triangle and a circle as metrically interchangeable.

## 6. Concurrent synthesis and three next problems

The final incoming operator packet also shares the15-arc Lucas/Paley reader
derived here. Its contribution is credited explicitly in LUCAS_MIRRORS.
Combining its singleton Gram with our two cyclotomic plane embeddings gives
the precise eigen-coordinate equations Bz=phi²z and Bz'=phi^(-2)z'. A separate
star family in GAP_TOPOLOGY now proves that constant threshold-space homotopy
type can coexist with return-Gram rank2^m-1. The final joint-port/germ packet's
phase hostile strengthens OPEN B below: a positive local germ alone need not
survive adding a new velocity.

The incoming Session91 work sharpened three boundaries used here:

1. The golden reader may preserve rotation while failing local-toggle
   lumpability. Its three-word even-period fibre repair is distinct from
   our twelve matching frames and the final two-frame observation ambiguity.
2. The split golden norm controls exact linear-observer memory. Its roots4
   and8 modulo11 explain why the geometric real golden ratio and the cyclic
   reader select different scalar channels here.
3. The geometric Euler repair preserves a whole unit-distance core while
   forcing E=V-1 by disconnected padding. Therefore E=V-1 alone cannot
   identify trijectivity, packing optimality, or geometric difficulty.
   See [the exact padding theorem](../session91_packing_golden_carriers_20261029/geometry/GEOMETRIC_EULER_REPAIR.md).

The earlier [tree compression](../session90_even_graphs_odd_descent_20261029/trijection/TREE_COMPRESSION.md)
requires its oriented frame. In our biplane, the missing frame is an explicit
matching, whose orbit is P1(F11). Both examples suggest retaining a separate
space of frames and its action before taking a quotient. They do not justify
a common unframed codec.

**OPEN A — geometric selector.** Can a finite unit-distance attachment with
named ports distinguish the two surviving biplane frames by a specified
boundary-colouring extension relation? The nonunit probe1->2 identifies the
exact bit to observe; realizing that bit geometrically is still unproved.
The base graph is two-colourable, so no six-colour lower bound follows.

**OPEN B — arithmetic gap transport.** Construct a marked, collision-compatible
chain model for an integer-velocity orbit that retains the labelled runner.
Test it on AP13 and the mixed-parity controls before asserting a nonzero
invariant. G1 supplies the ambient topology, but not an escape theorem.

**OPEN C — closure complexity.** For odd periods with prime L_n, classify
<translation, transported word reflection>. Period5 gives PSL(2,11),
period7 gives A29 with an explicit certificate. Which digit languages or
restricted translation menus admit a uniformly small frame orbit? A positive
answer must retain source labels and legal moves; counting residue classes
alone is insufficient.

The unifying research move is concrete: replace a numerical coincidence by
a named relation and action, then compute the smallest additional frame that
keeps the desired operations meaningful. Here this produces exact maps and
small certificates, together with a large next-period obstruction.
