# Black detours around the blue periphery

**PROVED:** all three neighbors of every trivalent pure-blue node belong to
one black connected component; each of its blue edges lies on a cycle in
the simple support graph. The contracted blue quotient therefore has a
triple parallel attachment at each such node. Also proved are the
nonstrong-source rotation lemma, its boundary-conditioned black hub, and a
mixed-center family whose blue edges admit three-black-line detours.

**FINITE-EXACT:** explicit five-vertex certificates, small-census replay
through order seven, all normal-form words of weight at most twelve, and a
minimal six-vertex failure of the fixed-boundary method with an actual
shortest black detour. **OPEN:** connectivity of all black-active nodes and
the existence or exclusion of mixed-to-mixed blue bridges in general.
No actual non-leaf multigraph bridge has been found here.

## 1. Current inheritance and conventions

The authoritative pure-blue alphabet is now the independently proved
prefix-rerouting classification in
`../pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md`.
Its branch theorem and counts are in that folder's
`PURE_BLUE_COUNTS_AND_BRANCH_GEOMETRY.md`.
The old affine-only proof gap is not used. The short rigid-strong lemma and
dependency audit in `../boundary_support_20261020_incoming.md`, section 6,
give an alternative route to the needed rigidity input.

The closest proved mechanism is **THM-4731**,
`../../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`:
black lines have distinct parallel sigma twins. The complete bridge
criterion is `../boundary_support_20261020_metagraph.md`, Proposition 1.
The least-used sidecar is a Hamiltonian path's *first two and last*
vertices, retained across two complement operations. The corrected near
miss is compressing those concrete path choices to a node degree.
The inherited abstract bridge hostile was insufficiently restricted by the
new pure-blue classification; section 7 below supplies a stronger one.

We use the converse-merged **line multigraph**, with each tiling pair an
individual edge and loops contributing twice to degree. Its simple support
graph discards loops and merges parallel edges. An individual black line
can never be a multigraph bridge; a black support edge can be a bridge.
These meanings are kept separate below.

Let `K_P(T)` complement every arc outside a directed Hamiltonian path P,
retaining the arcs of P, while keeping the original vertex labels.
This is the tiling involution tau in P-coordinates. Writing `E(P)` for the
unordered path edges,

    K_P(T) = T^op triangle E(P),
    K_Q(K_P(T)) = T triangle E(P) triangle E(Q).       (1)

Here triangle means reversing the specified edges. A path presentation is
blue exactly when its positional reversal is an anti-automorphism.
The board is: **cycle cuts / endpoint scores / path choices / black hubs /
pure-blue branches / mixed quotient**.

## 2. A black rotation at every nonstrong self-converse source

**Lemma 1 (PROVED).** Let U be a nonstrong self-converse tournament of order
`n>=4`, and let `P=(p_0,...,p_(n-1))` be any Hamiltonian path of U, blue or
black. Put

    S=K_P(U),   Q=(p_1,...,p_(n-1),p_0).

Then Q is a directed Hamiltonian path of S, its presentation is black, and

    K_Q(S)=U triangle {p_0 p_1, p_0 p_(n-1)}.          (2)

**Proof.** Hamiltonian paths visit the strong components in their unique
ordered sequence. The first component dominates the last, so
`p_0 -> p_(n-1)` in U. This is a nonpath arc; its reversal in S closes the
Hamiltonian cycle `p_0,...,p_(n-1),p_0`. Hence Q is directed. Its edge set
replaces precisely `p_0 p_1` by `p_(n-1) p_0`, proving (2).

Write d for outdegree in U. We claim

    d(p_0)+d(p_1)>n.                                  (3)

If the source component has size `s>=3`, self-converseness gives a distinct
sink component of the same size, so `2s<=n`. Both first vertices dominate
all `n-s` later vertices, and the first internal path arc contributes at
least one to their combined internal degree. Their sum is at least
`2(n-s)+1>n`. A tournament has no strong component of size two.

If the source component is a singleton, so is the sink component.
Then `d(p_0)=n-1`. The vertex `p_1` dominates both `p_2` and the last
singleton; these vertices are distinct because `n>=4`. Thus `d(p_1)>=2`,
again proving (3).

In the target V of (2), the degrees of the first and last Q-vertices are

    d_V(p_1)=d(p_1)+1,  d_V(p_0)=d(p_0)-2.

Their sum is at least n, whereas a blue Q-presentation would pair these
endpoints under an anti-automorphism and force their sum to be `n-1`.
Therefore Q is black in V, and also in S because tau preserves blue/black.
No parity assumption on n is needed. ∎

**Boundary-conditioned hub.** For one fixed U satisfying the lemma, every
Hamiltonian path with the same ordered triple `(first,second,last)=(a,b,z)`
has its complement partner joined by a black line to the same labelled
tournament

    V=U triangle {ab,az}.                             (4)

This is an exact target, not merely the assertion that the targets have
equal degree statistics.

There is another useful version of the mechanism. If S has a blue
Hamiltonian-cycle presentation of odd order n, and its one-step cyclic
rotation is also blue, the two positional reversals compose to the shift
by two. That is a transitive cyclic automorphism because n is odd.
Therefore S must be regular. A nonregular S has a black rotation.
The assumption is necessary: every cyclic cut of the standard Hamiltonian
cycle in a regular cyclic tournament can remain blue. Lemma 1 supplies a
stronger endpoint-score argument for the nonstrong-source situation.

## 3. All trivalent pure-blue branches have a common black neighbor

**Theorem 2 (PROVED).** At every order `n>=5`, let U be a trivalent pure-blue
node. Its three distinct mixed neighbors have a common black neighbor.
Consequently each of the three blue lines at U is nonbridging, even in the
simple support graph.

The inherited classification gives

    U=A -> C5 -> A^op,

where A is an ordered word of singleton and directed-triangle components.
Write `r=|A|`; the central C5 is the regular five-vertex tournament.

First suppose `r>=1`. In C5 on residues modulo five, use the three path
orbits represented by

    P0=(0,1,2,3,4),
    P1=(0,2,3,4,1),
    P2=(0,2,4,1,3).                                  (5)

All three begin at vertex 0. Fix the outer paths once and extend each
central path by them. These three full Hamiltonian paths have the same
first two and last vertices: for `r>=2` the first two lie in A; for `r=1`
the second is the shared first central vertex. The last lies in `A^op`.

U is nonstrong and self-converse, so Lemma 1 applies to all three paths.
Their three complement neighbors therefore have the same black hub (4).
The inherited intrinsic central score-shell argument distinguishes these
neighbors even after converse merging: the central five vertices are
exactly those of score `r+2`, and their induced Hamiltonian counts are
respectively 11, 15, and 13. Equivalently, the full branch geometry theorem
already supplies this distinction.

For `r=0`, order five, the following finite certificate supplies the same
conclusion. After expressing each partner in the coordinates of its source
path from (5), use the indicated black Hamiltonian path Q:

| source path | partner adjacency rows | Q | target relabelling |
|---|---|---|---|
| P0 | `(26,20,9,18,4)` | `(3,1,4,2,0)` | `(0,1,2,3,4)` |
| P1 | `(14,20,24,18,1)` | `(1,2,3,4,0)` | `(0,2,3,1,4)` |
| P2 | `(22,12,24,17,2)` | `(2,4,1,3,0)` | `(1,3,0,2,4)` |

An adjacency row is its outgoing-neighbor bitmask, with labels zero through
four. Form `tau(partner read along Q)` and then read its vertices in the
last column's order. In each case the result is exactly
`(26,20,25,18,0)`. This common hub has H=5 and is a strong four-vertex
tournament followed by a sink. All three Q-presentations are black;
the script checks their arcs and failed blue tests directly.

To conclude nonbridging, choose two distinct neighbors `S_i,S_j`.
Their common black hub gives a black support path of length at most two
between them. The hub cannot be U because U has no black incidence.
Together with `S_i-U-S_j`, this gives a cycle containing the blue edge
`U-S_i`. If the hub coincides with one of the neighbors, shorten the black
walk to its path; the argument still gives a cycle. This uses actual
different support edges, not the existence of parallel black twins.
Hence it applies both to individual line bridges and support bridges. ∎

### Exact consequence for the remaining blue quotient

Contract the black support components, retaining all blue edge instances,
to obtain the quotient Q of the boundary-support bridge criterion.
Theorem 2 says that each trivalent pure-blue node attaches to **one**
black-active quotient vertex by **three parallel blue edges**.
Previously the three neighbors were known to be distinct before
contraction; they were not known to contract to one component.

Thus all pure-blue nodes can be removed from the remaining bridge search:
the mass-one nodes contribute the known leaf bridges, and the mass-three
nodes contribute no bridges. Every additional multigraph bridge, if one
exists, must be a **mixed-to-mixed blue line**, and remains a bridge after
discarding the pure-blue periphery and contracting black components.
This does not prove that there is only one black-active component.

## 4. A positive mixed-to-mixed family

**Proposition 3 (PROVED).** Suppose

    U=A -> C -> A^op,   |A|>=2,

where C is self-converse and mixed. Every blue line incident at U has a
black-only detour of length at most three between its endpoints.

**Proof.** Start with the blue Hamiltonian path P giving that line.
The ordered components force each block to be visited consecutively.
The positional reversal preserves the central block, so its central path
is blue. Because C is mixed, choose a black Hamiltonian path of C and
replace just that middle path, obtaining P'. The outer paths are unchanged.
The full P' is black because its reversal fails already on the central
block.

Since `|A|>=2`, P and P' have the same first two and last vertices.
Their partners `S=K_P(U)` and `S'=K_(P')(U)` both join the common hub V
by black lines, by Lemma 1. The line `U-S'` is itself black. Thus

    U --black-- S' --black-- V --black-- S             (6)

is a black walk between the endpoints of the original blue line.
Removing repeated vertices or loops gives a path of length at most three.
It cannot use the original blue edge. ∎

The source U here is mixed. This is therefore an infinite family of actual
mixed-to-mixed blue edges proved nonbridging, not only a statement about the
pure-blue periphery. The condition on the padding is sufficient; no
necessity claim is made.

## 5. The first boundary obstruction, and an actual repair

The tempting universal assumption

> Every mixed nonstrong self-converse tournament has a black Hamiltonian
> path sharing the first two and last vertices of every blue path

is false. The first mixed nonstrong self-converse order is six. Let

    U=1 -> R4 -> 1,

where `R4` is the inherited tribonacci partner, with forward adjacent
arcs and backward longer arcs. Its adjacency rows are
`(62,36,40,50,38,0)`. Its five Hamiltonian paths are

    (0,1,2,3,4,5)   blue,
    (0,2,3,4,1,5)   black,
    (0,3,4,1,2,5)   blue,
    (0,4,1,2,3,5)   black,
    (0,4,2,3,1,5)   blue.

Fixing `(first,second,last)=(0,1,5)` leaves only the standard blue path.
The first failed implication is the assumed existence of a black path
inside that boundary fibre, not the hub formula or the rotation lemma.

There nevertheless is a shortest black-only detour of length three from U
to its standard blue partner. The following are literal successive
`tau(T read along path)` steps:

| source rows | black path | target rows |
|---|---|---|
| `(62,36,40,50,38,0)` | `(0,2,3,4,1,5)` | `(2,28,9,17,37,15)` |
| `(2,28,9,17,37,15)` | `(0,1,3,4,5,2)` | `(62,20,56,18,32,10)` |
| `(62,20,56,18,32,10)` | `(0,1,2,4,5,3)` | `(2,52,9,51,37,5)` |

Reading the last tournament in vertex order `(0,1,4,5,2,3)` gives
`tau(U)=(2,28,25,17,33,15)` exactly. A complete black-neighbor BFS on merged
classes certifies that there is no shorter detour for this edge.
These four nodes are distinct; with the original blue edge they form an
actual support four-cycle.

The repair frees the last endpoint in the second and third reroutings.
All displayed paths continue to begin at vertex 0, but the final endpoint
is no longer pinned. This is a concrete replacement for the failed
boundary-preserving move, not an all-order endpoint-release theorem.

### An internal-cut replacement for the failed boundary condition

There is also a structural explanation of a second shortest detour.
Let P and P' be Hamiltonian paths of a nonstrong tournament U of order
`n>=3`, with the same first and last vertices a,z, and sharing a directed
path edge `u -> v`.
Rotate each partner's Hamiltonian cycle to start just after that edge.
Both twice-complement targets are then exactly

    V=U triangle {az,uv}.                              (7)

The two rotated presentations both have endpoints v,u. Their endpoint
degree sum in V is

    d_U(u)+d_U(v)-[u=a]+[v=z].                         (8)

Reversing `uv` cancels in the sum; reversing `az` gives the displayed
correction. If (8) differs from `n-1`, both rotations are black. Hence a
blue P and black P' satisfying this shared-edge imbalance condition give
the same three-black-line detour as (6). This is a sufficient certificate;
equality in (8) does not imply that a rotation is blue.

For the six-vertex hostile, use the black path
`P'=(0,4,1,2,3,5)`. It shares endpoints 0,5 and the internal edge `1 -> 2`
with the standard blue path. The partner rotations are respectively

    Q=(2,3,4,5,0,1),  Q'=(2,3,5,0,4,1).

Both reach `V=U triangle {05,12}=(30,32,42,50,38,1)`.
The endpoint degree sum is `d_U(1)+d_U(2)=4`, different from `n-1=5`.
Thus the obstruction is repaired by moving the *cut* to an internal
shared edge, while keeping the original source paths' endpoints fixed.
The script checks the literal matrices and both black paths.

The internal score test can be uninformative when the central block is
regular: all central vertices then have the middle score, so their
internal edges have the balanced sum `n-1`. The boundary cut can still
work; no impossibility theorem follows from this test's silence.
The generalized shared-edge condition succeeds on all 14 mixed nonstrong
blue presentations through order seven. It is **not universal**.
A cheap exhaustive extension generates only nonstrong blue sources:
choose their first strong component S, an arbitrary blue middle, and the
forced sink presentation `S^op`. Uniqueness of the ordered SCC list makes
this a complete presentation generator, without traversing the full tiling
cube. It yields:

| order | nonstrong blue presentations | mixed | fixed-triple failures | shared-edge-test failures |
|---|---:|---:|---:|---:|
| 8 | 70 | 67 | 5 | 3 |
| 9 | 519 | 511 | 9 | 3 |

The first obstruction is `U=R4 -> R4^op`, of order eight, with rows
`(242,244,249,243,32,64,144,48)`. Fixing the standard first and last vertices
leaves only the standard Hamiltonian path, which is blue. Hence there is
no black replacement path with those two endpoints, irrespective of the
choice of cut or the degree test.

This is an infinite structural hostile family:

    U=R4 -> Z -> R4^op,   Z any pure-blue tournament, including empty.

In R4, a Hamiltonian path starting at its standard first vertex is forced
to be its standard path: the first vertex's only outneighbor is the
second, and the remaining unused path is then forced. The dual statement
holds for the last vertex. Thus pinning both global endpoints forces
both outer paths, while every central path is blue by purity of Z.
Every resulting full path is blue. U is nevertheless mixed: choose
unmatched paths in the two outer R4 components. More explicitly, R4 has
trivial automorphism group and five Hamiltonian paths. Once one outer path
is fixed, only its prescribed reversed counterpart can make the pairing
of those two components arc-reversing; choosing any other outer path
produces a black presentation.

The missing operation is therefore genuinely an endpoint change.
Enriching a state while requiring both source endpoints to stay fixed
cannot repair this family. A general endpoint-release detour remains
open; the six-vertex certificate above proves only its stated instance.

A related cheap hostile explains why the trivalent proof required a global
cycle cut: at order five the three central partners have respectively
`1,3,1` Hamiltonian paths with both standard endpoints fixed. The first and
third have no local rerouting at all under those two fixed ports.

## 6. What the small censuses actually say

The checker reuses the committed THM-4731 class and line tables; it does not
rerun the order-eight cube. At orders three through seven the numbers are:

| n | individual multigraph bridges | simple-support bridges |
|---|---:|---:|
| 3 | 1 | 1 |
| 4 | 1 | 2 |
| 5 | 2 | 3 |
| 6 | 2 | 4 |
| 7 | 3 | 6 |

The extra support bridges cannot be used as counterexamples to a claim
about individual line bridges. In particular a doubled black edge can be a
support bridge while neither of its individual black lines is a bridge.

The trivalent node at order five has mixed neighbors labelled `(2,3,5)` in
the inherited table and common black hub 4. At order seven its neighbors
are `(145,199,296)` and common black hubs include
`187,201,227,238,241`. These labels identify that committed census only;
the all-order proof does not extrapolate from them.

## 7. A stronger unrealized mixed-quotient hostile

The theorem does not force black-core connectivity from colored degrees.
Here is an explicit **unrealized abstract rewiring**, using the actual
order-five merged node labels, actual fibre masses, and actual blue/black
degrees from the inherited table.

Keep every actual pure-blue incidence:

    blue (0,11), (2,9), (2,6), (3,6), (5,6).

Take the other blue edges to be `(0,2)`, `(0,1)`, and the loop `(3,3)`.
Set black multiplicities to

    (0,1):6,  (1,7):2,
    (2,4):4, (3,4):2, (4,5):4, (2,3):6.

All omitted multiplicities are zero. Every node retains its original blue
degree and black degree, so all actual merged fibre masses are retained.
There are eight blue and 24 black line instances, as in the actual
order-five staircase. Every black multiplicity is even. The pure-blue
branch 6 retains precisely its three actual neighbors and their common
black hub 4. The known strong/nonstrong separation is also respected.

Nevertheless `(0,2)` is the sole line between the two black-active sides,
and is a multigraph bridge with two mixed endpoints. The other bridges
are the two known pure-blue leaf edges `(0,11)` and `(2,9)`.

This is not a tournament counterexample. It is a witness that even the
actual fibre inventory, the full colored degree inventory, the exact
pure-blue attachment data, and the new common-hub condition do not imply
the desired conclusion. Further realization constraints, such as the
matrix-valued distance conditions in
`../merged_metagraph_distance_moment_constraints_20261020/PROOF.md`, have not
been certified for this rewiring. They must not be silently assumed.

In fact the subsequently read endpoint-current theorem already rejects
this fake. In
`../merged_metagraph_endpoint_transport_20261020/endpoint-transport-and-complement-leakage.md`,
equation (2.1), each colored row's Laurent polynomial in the change of
cyclic-triangle count is fixed by the source's Hamiltonian endpoint-score
inventory. For our actual census labels the blue rows change as follows:

    node 0: z^-3+2z      becomes z^-3+z+1;
    node 3: 2z^-1+z      becomes 2+z.

The new constant terms are already incompatible with blue-current parity
at odd order. The code computes these discrepancies from the actual score
multisets, counting the fake loop twice. Thus this hostile shows a failure
of the specifically listed degree/periphery data; it does **not** survive
the incoming endpoint-current gate, much less all known realization
constraints.

The remaining precise question is whether every mixed-to-mixed blue line
admits a black-only detour, or at least a detour avoiding that line after
black contraction. The first assertion would prove that blue never joins
different black-active components; it is stronger than merely excluding
bridges. Neither assertion is proved in this packet.

## Reproduction

    python 05-knowledge/results/defect_transitions_20261020/blue_check.py
    python -O 05-knowledge/results/defect_transitions_20261020/blue_check.py

The script uses exceptions rather than removable assertions. It checks the
five-vertex literal certificates; all 187 singleton/triangle words of
weights 1 through 12 and all their 561 three-branch routes, reaching order
29; regular rotation hostiles; all Hamiltonian paths of the small
nonstrong blue presentations through order seven; the mixed-center
bypasses and the shortest six-vertex detour; all 589 nonstrong blue
presentations at orders eight and nine using the canonical SCC generator;
and the full inherited small line tables. The generator is independently
compared with the full blue cubes through order seven. The abstract hostile
is checked against both the constraints it preserves and the incoming
endpoint-current gate that rejects it. The infinite conclusions are
supplied by the proofs above.
