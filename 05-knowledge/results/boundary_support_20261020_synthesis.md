# Boundary-support reductions, Petersen, and tournament metagraph wiring

**20 October 2026. PROVED scoped reductions and identities; FINITE-EXACT
specified censuses; CITED Hamiltonian-path bound; general snark reduction,
all-order degree extrema, and compact tournament realization remain OPEN.**

The main advance is an actual reduction: a four-claw strip can be replaced
by a two-claw strip with exactly the same labelled colouring support, and
the replacement is a graph minor. It reduces every obstructed cubic
claw-ring wiring to Petersen. The tournament connection is also explicit:
both settings retain labelled states, pair them, and then contract by a
fibre map. Their common mechanism is compatibility, not graph isomorphism.

## 1. Inheritance that changed the research

The [board](boundary_support_20261020_board.md) records the portfolio,
hostiles, sidecars, and preserved predicates. Three parallel investigations
read existing work before deriving, and then cross-audited their results.

- [THM-4116, boundary-state gluing](../../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md)
  supplied nonnegative extension vectors and the Petersen four-pole.
- The [earlier Pythagorean synthesis](graph_pythagorean_20261018_synthesis.md)
  separated the literal Seidel 3,4,5 identity from the colouring geometry.
- The incoming [105 matching atlas](hamming_matching_atlases_20261018e.md)
  identified the user's exact seven / forty-two / fifty-six decomposition.
- The [snark transfer packet](snark_boundary_transfer_20261018a.md)
  supplied authenticated Petersen subdivisions in the Blanusa examples;
  those certificates were reused and checked, not rediscovered.
- [THM-796, three-sorted incidence](../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md)
  supplied the tournament presentation/line/node span.
- [THM-4731, self-lines and sigma-twins](../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md)
  supplied doubled black multiplicities, the committed censuses, and the
  corrected tribonacci partner of the transitive class.
- [THM-840, continuation congruence](../../01-canon/theorems/THM-840-hamming-five-continuation-congruence-boundary.md)
  gave the exact criterion for when a quotient retains an operation.
- The [weighted-core inheritance audit](weighted_core_transport_20261020/inheritance.md)
  prevents use of the provisional all-order pure-blue alphabet. None of
  the proofs below depends on that classification.

The older bucket reflections correctly anticipate degree-constrained
wiring, but their proposed partner law `2^(n-2)+1` is false. The current
tribonacci theorem takes precedence.

## 2. The obstruction-preserving Petersen reduction

A claw module has three incoming and three outgoing edges. Its transfer
matrix T is indexed by the 27 ordered colour triples. The missing colour
at each leaf is its spoke colour; requiring the three spokes to differ
gives every matrix entry directly.

The complete state graph is

    2 K3 + 3 C6 + 3 isolated vertices.

Consequently T squared consists of eight active three-state blocks, each
equal to `I+J`, with three zero rows/columns. On those same blocks,

    T^2=I+J,       T^4=I+5J.

They have identical support. Delete the centres of the middle two claws
in a four-claw strip and suppress their six degree-two leaves. This gives
the two-claw strip as a minor, retaining all six outer ports.

For **every exterior**, a boundary word extends across one strip iff it
extends across the other. This proves preservation of colourability and
uncolourability, not merely equality of a statistic.

For the flower graph J_k, repeat this move and contract the final triangle:

    J_k -> J_(k-2) -> ... -> J_3 -> Petersen       (odd k>=5).

The last triangle and a cubic vertex have the same three-port relation:
three distinct boundary colours, with one extension. The terminal
Petersen isomorphism is explicit.

The result covers all cubic rings of k claws with arbitrary permutations
of their three strands at the junctions. Independent leaf relabellings
move all twists into a single closing permutation pi. The graph is
uncolourable precisely when **k is odd and pi is odd**. In S3 every odd
permutation is a transposition, so every obstructed ring is isomorphic
to an odd flower and has this reduction.

This is a complete answer for that class. It is not a proof that every
snark has such a strip.

The [full proof and edge formulas](boundary_support_20261020_snarks.md)
also establish the four adjacent-edge multiplicities, for odd k:

| Edge | m_e |
|---|---:|
| b-spoke | `(2^k+4)/3` |
| c/d-spoke | `(2^k-2)/3` |
| b-cycle edge | `(2^(k-1)-4)/3` |
| c/d-cycle edge | `(2^(k-1)+2)/3` |

Each complementary four-pole has signature `(m_e,m_e,0,0)`.
All three elementary two-vertex smoothings become colourable when
`m_e>0`, with counts `9m_e,3m_e,3m_e`. This explains why the larger
six-port move is needed.

## 3. The arithmetic now comes from actual extension multiplicities

On each active block,

    T^(2h)=I+alpha_h J,       alpha_h=(4^h-1)/3.

Its distinct eigenvalues are 1 and `1+3alpha_h=4^h`. Applying Euclid's
map to `(2^h,1)` gives the primitive triples

    (3alpha_h, 2^(h+1), 3alpha_h+2)
      = (4^h-1, 2^(h+1), 4^h+1).

They begin `(3,4,5)`, `(15,8,17)`, `(63,16,65)`. The arithmetic parameter
is a proved extension count in the very transfer matrix used by the
minor reduction.

The isosceles connection is a different exact chart. For nonzero boundary
vectors f,g put U=f tensor g and V=g tensor f. Their sum and difference
are orthogonal, with squared lengths

    2(M+C^2), 2(M-C^2),    M=||f||^2||g||^2, C=<f,g>.

Thus their possibly degenerate right triangle is nondegenerate and
isosceles iff the glued graph has no colouring. The reduction preserves
this zero predicate. It changes
positive multiplicities, so it need not preserve lengths. The integer
triples above flatten with h; they do not approach the isosceles limit.

## 4. What seven / forty-two / fifty-six actually does

The eight addresses are F2^3. Their 105 perfect matchings split into
7 parallel, 42 two-direction, and 56 four-direction matchings. For an
eight-pole and matching M, define

    B_M(s)=product_({x,y} in M) [s_x=s_y].

Then its closure count is exactly `<f,B_M>`. This turns the inherited
atlas into actual wiring operations.

The [complete matching analysis](boundary_support_20261020_seven.md) finds:

- The 105 ternary kernels have rank 91. Their 14 linear dependencies
  are generated by four-by-four determinant identities.
- The 98 nonparallel kernels already span all 105. A signed formula
  recovers the seven parallel measurements, but its negative terms
  cannot be read as individual obstruction certificates.
- The seven parallel probes miss 1,176 of the 1,641 parity-allowed words.
- All 105 closures of each named Petersen/prism eight-pole were tested.
  In each bank, two closures are Petersen and nine further uncolourable
  closures have a bridge. No other closure is uncolourable.

At tournament order four, the tiling cube itself is F2^3. The very same
105 pairings push forward through its ordinary fibres `(1,1,1,5)`, or
its converse-merged fibres `(1,2,5)`. Six distinct pairings satisfy the
staircase involution constraints and produce the same coloured
metagraph. This is a literal common carrier and an explicit example of
the information lost on passing to the quotient.

A late incoming defect-grammar proof supplies another exact dictionary:
these 105 matchings are the 105 leading feedback-budget-four tournament
templates, contributing `105*x^53/(1-x)^9`. The
[incoming integration](boundary_support_20261020_incoming.md) audits the
normalization and shows why the rank-91 colouring quotient cannot replace
105 in the tournament class count. Sharing a defect endpoint fails the
naive equality-kernel translation; one new cubic vertex repairs its
seven-port parity, giving a concrete next operation to test.

## 5. Tournament bridges and complete wiring

Throughout this section the graph is the **converse-merged complement-line
multigraph**, retaining loops and parallel lines, at order n>=3.

Every black line has a parallel twin. Therefore no individual black line
is a bridge. Contract each connected component of the black support graph,
including black-isolated vertices, and keep every blue line separately.

**A blue line is a bridge of the original graph iff its image is a bridge
of this contracted blue multigraph.**

This necessary-and-sufficient criterion covers non-leaf bridges. It does
not assume that all bridges touch leaves. Exact replay gives bridge
counts `1,1,2,2,3` at n=3,...,7, all leaf-incident. The inherited order-eight
census implies three leaf bridges. An explicit six-node abstract wiring
satisfies all the parity and twin laws but has a non-leaf bridge, showing
that those laws alone cannot prove a universal leaf classification.

The incoming strong-core work sharpens the even-order conclusion:
pure-blue nodes are exactly leaves, attached to mixed strong nodes.
Delete these leaves from the contracted blue graph Q. Every remaining
bridge corresponds precisely to a non-leaf bridge of the original
metagraph. The audited proof and its cited Lu-theorem dependency are
recorded in the [incoming integration](boundary_support_20261020_incoming.md).
No connectivity of the remaining black-active quotient is assumed.

There is a constructive abstract realization theorem. A blue multigraph
supported on self-converse nodes together with an even-multiplicity black
multigraph can be realized by commuting involutions: two flags per blue
line and four flags per black twin pair. Specified fibre sizes impose
exact degree conditions. Matching the staircase orbit counts gives an
isomorphic Klein-four action.

Actual tournament realization requires more. An order-five edge swap
preserves every actual fibre size and coloured degree but sends the
transitive node to the wrong partner. A sufficient actual certificate
retains a partition of the full tiling cube, tournament isomorphism
witnesses for its fibres, and every complement pair. The compatible-face
recursion of THM-796 supplies an inductive version.

Incoming work gives useful additional necessary screens: nonstrong
nodes have only strong neighbors, and blue lines obey a triangle-count
current law. Our false order-five wiring survives the first screen but
fails the second, changing allowed jumps 3 and 1 into forbidden 4 and 2.

If "complete metagraph" means a complete underlying simple graph, it is
impossible at n>=4: the transitive node has degree one, and there are at
least three nodes.

There is also an exact measure of quotient loss. Write J for the matrix
of fibre indicators, P for complement, F=diag(fibre sizes), and
A=J^T P J for line incidence. Then

    F-A F^(-1) A = R^T R >= 0,       R=PJ-J F^(-1) A.

The defect vanishes exactly when complement descends to an involution
of whole fibres. More generally its rank is `r-c-b`, where r is the
node count, c the number of support components, and b the number of
bipartite components. The kernel consists of functions constant on
even-walk classes. This measures precisely how much of the complement
operation is lost by retaining only the node.

Proofs and witnesses are in the
[metagraph note](boundary_support_20261020_metagraph.md).

## 6. D_n(z), M_n and the factorial normalization

Because the original notation was not uniquely recoverable, this session
explicitly adopts

    D_n(z)=sum_v z^deg(v),       M_n=deg D_n.

Write h(T)=H(T)/|Aut(T)| for an ordinary tournament class. Then

    D_n(z)=sum_(self-converse T) z^h(T)
            +(1/2)sum_(non-self-converse T) z^(2h(T)).

Odd coefficients count self-converse fibres; exponents 2 mod 4 count
paired nodes; positive multiples of four never occur. Exact formulas
for every coefficient, an all-order finite enumeration algorithm, and
the full coefficient dictionaries through n=8 are saved.

For example,

    D_3=2z,
    D_4=z+z^2+z^5,
    D_5=2z+z^2+z^3+z^5+2z^9+z^10+z^11+z^13.

| n | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|
| M_n | 1 | 5 | 13 | 86 | 306 | 1322 |
| Maximum H(T), for comparison | 3 | 5 | 15 | 45 | 189 | 661 |

The automorphism divisor and converse pairing matter: these are different
extremal sequences.

The [asymptotic proof](boundary_support_20261020_asymptotics.md) gives

    (1+o(1)) n!/2^(n-2) <= M_n <= C n^(3/2) n!/2^(n-2).

The lower bound uses an elementary Burnside derivation of the class-count
asymptotics and total degree; the upper bound uses
[Alon's Theorem 1.1](https://www.tau.ac.il/~nogaa/PDFS/hamilton.pdf).
Therefore

    M_n/n! = 2^(-n+O(log n)),
    (M_n/n!)^(1/n) -> 1/2,       M_n/n! -> 0.

An exact all-order extremal formula and an asymptotic equivalent with a
leading constant remain open in this work.

## 7. Verification and remaining frontier

| Computation | Explicit universe / control |
|---|---|
| Independent strip audit | Two constructions of all 729 entries; exact matrix powers and determinant |
| Flower colours | Independent perfect-matching counts for J3 through J15, up to 60 vertices; even/untwisted positive controls |
| Minor certificates | Explicit chains through J101, 404 vertices; authenticated rooted Blanusa models |
| Eight-port closures | All 105 matchings, all 6,561 ternary words, Petersen and prism controls; independent binary-flow counts |
| Tournament wiring | Full inherited line tables through n=7; independent small tiling construction; degree coefficients through n=8 |
| Proof review | Root and agents independently checked support, asymptotics, quotient types, and scope |

The five scripts and receipts are linked in the component notes. Normal
and optimized runs agree; no claim rests on assertions removed by `-O`.

The next serious obligation is to find an unavoidable, support-preserving
minor reduction outside the claw-ring class. Identical four-pole support
does not suffice: two subdivided K4s joined by a bridge have Petersen's
support but no Petersen minor. Likewise, the tournament degree polynomial
does not determine wiring. The new reductions and certificates make both
missing steps explicit; they do not settle the general Kuratowski/Tutte
program.
