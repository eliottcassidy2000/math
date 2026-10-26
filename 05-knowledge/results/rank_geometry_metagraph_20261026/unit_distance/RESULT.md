# Four-state compatibility with the complete planar metric retained

2026-10-26. **PROVED below:** a sharp rank obstruction for the inherited
two-loop transport; an explicit planar fixed-free transport after changing
the pairing; exact full boundary counts for a unit wheel and an all-length
induced diamond strip. **FINITE-EXACT:** the fifteen-pairing classification
and the declared independent arithmetic/coloring controls. **OPEN:** an
unprescribed six-color obstruction from these gadgets. Here that phrase
means a graph with no proper six-coloring, hence a seven-chromatic witness.
No such witness is claimed. The [independent parent review](AUDIT_PARENT.md)
accepts the proofs and declared finite controls within these scopes.

The main positive finding is that three of the fifteen pairings of an actual
regular hexagon admit a parity-preserving, fixed-point-free four-state
transport by reflection. The diameter pairing does not. However, complete
coloring tables expose the remaining gap: an isometry of a port set is not
a forced relation between color choices.

## 1. Inheritance and the precise source problem

Closest proved mechanism: the `C2 x S4` reduction in
[Fano/eight, §§3–5](../../hadwiger_six_seven_eight_20261025/fano_eight/RESULT.md).
The canonical hostile there consists of two affine maps on `A=F2^3`,

    g1(x)=pi(x),       g2(x)=pi(x)+101,
    pi: bit positions 0 -> 1 -> 2 -> 0.

Their fixed sets are `{000,111}` and `{011,100}`. Each loop admits a
frame of either completion parity, but the two loops have no common frame.
The even half `E={000,011,101,110}` has four states; after completion
parity closes, the residual action is S4. This is an exact abstract seam
system, not yet a graph-coloring relation.

The [Session61 palette-deck proof](../../session61_hadwiger_six_seven_20261025/palette-deck/RESULT.md)
retains the existential ternary assignment and full port relation.
The [cube boundary proof](../../hadwiger_six_seven_eight_20261025/cube_boundary/RESULT.md)
gives exact table composition and distinguishes feasibility from counting
marginals. The [six-color stress gate](../../strong_blue_games_20261020f/six_colour.md)
separates prescribed failure from ordinary chromatic number, and uses
circle intersections to exclude an abstract minimum blocker bank.
These are inherited mechanisms, not new results here.

Corrected near miss in this investigation: the common-frame property of
the **diameter pairing** cannot be generalized to all pairings of the same
hexagon; §4 gives a counterexample. Least-used sidecar: the full boundary
Gram matrix, including nonedge distances. Anchor: realize an actual seam;
niche: quantify rank obstruction; wildcard: change the pairing without
changing the points. The board is full color relation / full metric /
pairing / four-state holonomy / color quantifier / gluing operation.

For current external status, Session61's
[source record](../../session61_hadwiger_six_seven_20261025/SOURCE_STATUS.md)
cites the September 2026 six-color lower bound. This session reopened the
[pinned primary Lean source](https://raw.githubusercontent.com/openai/math/adc7f1241b42e322a6451854ab7e4b4c146bf78a/lean/OAI/Geometry/PlaneColoring/Five.lean)
and its paper link; no build or transitive axiom audit was performed.
That lower bound remains **CITED**, and none of the deductions below
depends on it. Older five-to-seven status prose in the stress note is
historical. A new finite seven-color witness is a separate target.

## 2. Full boundary composition requires a common geometric realization

For a finite graph F with named boundary B and palette `[q]`, retain

    w_F(beta) = #{proper q-colorings of F extending beta},
    beta in [q]^B.

Every assignment is retained, including zero entries and all color-label
patterns. If F and G have exactly B in common, disjoint private vertices,
and no undecomposed cross edges, their literal union has
`w_(F union G)(beta)=w_F(beta)w_G(beta)`. Forgetting an internal port is
summation over its color. For two-terminal tables, series composition is
matrix multiplication and parallel composition is entrywise product.
Additional edges contribute indicator factors only when both endpoint
colors are retained; otherwise the interface must be enlarged.

For the fixed four-state transport problem the seam matrix is
`R_g(x,y)=[y=g(x)]`. Closing a loop reads its diagonal; sharing a root
between two loops gives `sum_x R_g(x,x)R_h(x,x)`, not the product of
their separate traces. The inherited g1,g2 give zero although both traces
are positive. The checker retains all sixteen entries of each relation
and checks all 24^2 residual permutation compositions.

To realize geometric fragments jointly, their boundary distance matrices
must also agree under the proposed identification. After placing them,
check injectivity and every cross unit contact. Separate feasibility of
two boundary metrics is insufficient, as §3 proves. Separate absence of
cross edges inside each fragment is also insufficient: unit edges from
0 to 1 and from 0 to alpha, where `alpha=(1+i sqrt(3))/2`, create the
third unit edge. At a fixed center color, the independent product count
is 25 for six colors; the actual induced triangle count is 20.
The product formula remains correct for the noninduced two-edge graph;
it is the change to the complete unit-distance graph that adds the factor.

## 3. The inherited two-loop candidate needs dimension at least three

**Theorem U1.** Let six distinct Euclidean points be named `p_(j,b)`,
`j=0,1,2`, `b=0,1`. If both g1 and g2 above act on these labels by
isometries of the same point set, its affine span has dimension at least
three. This conclusion concerns the complete metric, not merely the
unit-edge graph or permutations of color names.

**Proof.** The composition `g2 g1^-1` flips endpoint pairs 0 and 2.
Conjugating by g1 gives the other two double flips. Write
`d_j=p_(j,1)-p_(j,0)`. Each double flip negates the two corresponding
difference vectors and fixes the third. For distinct i,j choose a
double flip which negates d_i and fixes d_j. Isometries preserve inner
products of differences, so `d_i dot d_j=-(d_i dot d_j)=0`.
All three d_j are nonzero by injectivity. Thus their Gram matrix is
positive diagonal and has rank three. The differences lie in the affine
span of the six points, proving the claim. QED.

No central symmetry, equal radii, unit edges, or generic position is
assumed in this proof. It already follows from the three double flips.
The cyclic g1 also makes the three squared lengths equal. In six
symmetric Gram variables, each loop alone gives four independent linear
conditions; together they give five, leaving precisely scalar multiples
of `I3`. Positivity and distinct paired points then force rank three.

**Sharpness and two separately valid planar certificates.** The points
`p_(j,0)=e_j`, `p_(j,1)=-e_j` realize both maps in R3. These are unit
vectors; their nonantipodal mutual distances are sqrt(2), not one. Scaling
by `1/sqrt(2)` gives unit edges of the octahedron graph if that graph is
desired. The no-go itself is about complete-metric isometries.

For g1 alone, take the planar hexagon `p_(j,b)=(-1)^b u_j`, with
`u_j=exp(2pi i j/3)`. For g2 alone, interchange the two endpoints of
pair 2 in that same coordinate list. Each is an injective unit-radius
planar metric. Their difference-vector Grams are respectively

    [ 4 -2 -2 ]             [ 4 -2  2 ]
    [-2  4 -2 ]     and     [-2  4  2 ].
    [-2 -2  4 ]             [ 2  2  4 ]

Both have rank two and are positive semidefinite. Neither preserves both
loop actions. This is an explicit local-to-global failure of compatible
metric realization, not an assertion that either loop is individually
nonplanar. Arbitrary graph automorphisms are outside U1: an induced unit
star has every permutation of its leaves as a graph automorphism, even
though most do not preserve its leaf-to-leaf distances.

## 4. Changing the pairing creates a genuine planar hidden transport

**Diameter pairing.** In the hexagon just used, its twelve isometries act
on paired endpoints as `S3` pair permutations and simultaneous reversal
of all three pairs. Thus they form `C2 x S3` inside the 48 palette maps
`C2 x S4`. After completion parity is removed, all fix the even frame
which selects one alternating triangle. Every family of such residual
holonomies has a common frame. The proposed two-loop obstruction cannot
come from these rigid diameter-palette seams.

**Theorem U2 (a realized positive alternative).** Number the same hexagon
vertices clockwise by `z_j=exp(pi i j/3)`. Instead choose

    P0={z0,z3},    P1={z1,z5},    P2={z2,z4}.

Complex conjugation fixes P0 pointwise and swaps both P1 and P2. Its
endpoint-frame action is `x -> x+110` (integer mask 6 under the named
bit convention). It preserves completion parity, but has no fixed point
on either four-state parity half. It is an actual complete-metric
isometry of six distinct planar points. The pair squared lengths are
`(4,3,3)`. This is not the diameter pairing, whose lengths are `(4,4,4)`.

The proof is the displayed coordinates and their conjugates. Thus a
planar hidden four-state obstruction is geometrically possible; U1
excludes its particular two-order-three-loop candidate, not all possible
four-state transports.

**Complete small classification (FINITE-EXACT).** For each of the fifteen
perfect matchings on this fixed hexagon, retain its full dihedral
stabilizer, the subgroup preserving completion parity, and the common
fixed frames in E. All quantities use one matching throughout.

| Squared pair lengths | Matchings | Stabilizer size | Parity-preserving size | Common E frames |
|---|---:|---:|---:|---:|
| 1,1,1 | 2 | 6 | 3 | 1 |
| 1,1,4 | 3 | 4 | 2 | 2 |
| 1,3,3 | 6 | 2 | 1 | 4 |
| 3,3,4 | 3 | 4 | 2 | 0 |
| 4,4,4 | 1 | 12 | 6 | 1 |

The three zero-common-frame cases are rotations of the mirror-pair
construction. This changes the actual geometric gate while preserving
the six points. Pairing is therefore an essential variable, not cosmetic
relabeling after a fixed pairing has been selected.

The same reflection has permutation type `2^2 1^2` on the six physical
ports. The [information lane's literal Frobenius boundary](../information/RESULT.md#5-a-literal-frobenius-boundary-for-the-new-planar-reflection)
proves that every tournament orientation on those ports changes at least
two arcs under it, sharply: its complete defect polynomial is
`512 t^2 (1+t^2)^6`, independently checked on all 32,768 orientations.
This retains the actual vertices and permutation but changes the predicate
from metric symmetry to invariance of a strict complete binary relation.
It supplies neither a forced color relation nor a chromatic obstruction.

## 5. Full unit-wheel relation exposes the remaining color gap

Add the center to the six hexagon points. The complete induced unit graph
is the wheel with rim C6: rim neighbours are at distance one, second
neighbours at sqrt(3), and opposites at two. Its exact unprescribed
six-color boundary table, with the rim as six named ports, is

    w(beta) = [beta is proper on C6] * (6 - |image(beta)|).

There are 14,910 boundary assignments of positive weight and 24,600 total
six-colorings. These are complete counts, not pair projections or one
prescribed palette test. The wheel is three-chromatic: a rim two-coloring
and a third center color give an explicit coloring, and a triangle gives
the lower bound.

In the mirror pairing of §4, a frame can be used to assign the two named
colors `2j,2j+1` bijectively to the two points of P_j. Every one of these
eight assignments is rainbow, so its wheel weight is zero. Nevertheless
the graph has the three-coloring above. Eliminating every frame in this
rainbow chart does not eliminate every proper six-coloring. Conversely,
on the rim alone every one of those eight assignments has weight one.
Neither full relation imposes the claimed reflection transport as a
coloring constraint. U2 supplies the metric carrier; the seam gadget
forcing that relation is still missing.

An independent rational unit-star control makes the same quantifier
visible without rim edges. Its ports are
`((1-j^2)/(1+j^2),2j/(1+j^2))`, j=1,...,6. The exact induced graph is
K_(1,6), its boundary weight is `6-|image(beta)|`, and only the 720
rainbow assignments have zero weight. Total six-colorings are 93,750,
although the graph is bipartite. For a three-port star, choosing one
named color from each of three disjoint pairs always leaves exactly
three center colors: all eight frames, hence both four-state halves,
extend. A frame transport needs actual constraints beyond this carrier.

## 6. An all-length realized family with full boundary transfer

**Theorem U3.** For every r>=1 place terminals
`U_j=(sqrt(3)j,0)`, j=0,...,r, and two private vertices per interval,

    A_j=(sqrt(3)(j+1/2), 1/2),
    B_j=(sqrt(3)(j+1/2),-1/2),   j=0,...,r-1.

The complete induced unit graph is a serial chain of r unit diamonds,
with `3r+1` distinct vertices and exactly `5r` edges. For q>=3, its full
two-terminal coloring table is `D_q^r`, where

    D_q(a,b) = (q-1)(q-2)  if a=b,
               (q-2)(q-3)  if a!=b.

**Proof of geometry.** Represent coordinates as `(sqrt(3)a/2,b/2)`.
The integer pairs are `(2j,0)` and `(2j+1,±1)`. Unit distance is exactly
`3(Delta a)^2+(Delta b)^2=4`. Integer solutions among these points have
`Delta a=0, |Delta b|=2`, or `|Delta a|=1, |Delta b|=1`. These are
precisely the private edge A_j B_j and its four incident terminal edges.
No other contact or collision occurs, for any r.

**Proof of color table.** In one diamond the two private vertices are
adjacent and adjacent to both tips. Equal tip colors leave q-1 choices
then q-2; unequal tip colors leave q-2 then q-3. Series gluing sums over
the exact shared terminal color, giving matrix multiplication. No pair
shadow or prescribed leaf color is used.

At q=3, `D_3=2I`: this is the familiar three-color equality wire. At q=6,

    D_6=8I+12J,
    D_6^r=8^r I + (80^r-8^r)J/6.

All entries are strictly positive for every r>=1. The r=2 diagonal and
off-diagonal weights are 1120 and 1056. Thus arbitrary extension of
this actually realized diamond chain cannot impose any six-color
two-port restriction. Positive tables also remain positive under series
and parallel composition; that statement is algebraic and does not
assert every resulting abstract parallel assembly has a planar unit
realization. This is a precise stopping rule for this grammar, not a
no-go for all spindle-based or high-chromatic gadgets.

## 7. Next constructive target and its validation interface

The mirror pairing is the stronger next geometric candidate: it supplies
a real parity-invisible transport that the diameter model could not.
For a proposed seam fragment F, keep **both** six-port boundaries,
compute its entire `[6]^(B_in union B_out)` relation (or a lossless
color-permutation quotient), and identify which full colorings induce
the proposed four-state channel. The desired transport must follow from
that relation, with every alternate channel/chart accounted for.
Attaching two algebraic permutation matrices without such a derivation
repeats the inherited gap.

A geometric certificate must retain the same port Gram data, an injective
embedding, and every actual cross contact after gluing. U1 rejects one
candidate already at that stage. U2 passes that metric stage but not the
color-relation stage. U3 shows that extending the simple diamond wire
cannot repair the latter at six colors. These three conclusions narrow
the construction problem without claiming it is solved.

| Map | Preserved | Lost / required sidecar |
|---|---|---|
| Actual coloring -> selected frame chart | only if an explicit extraction is proved | other proper colorings; the wheel hostile |
| Full boundary table -> pair shadows | pairwise feasibility | joint obstruction and extension multiplicities |
| Geometric seam -> palette permutation | complete metric under that isometry | arbitrary color names are not geometric points |
| Separate metrics -> one common interface | nothing without equality/gluing check | U1's rank-three obstruction |
| Unit graph -> abstract automorphism | adjacency | nonedge distances; unit-star hostile |
| Diamond chain -> matrix power | all endpoint-color extension counts | prescribed three-color equality does not survive q=6 |

## 8. Reproduction and declared universe

Run from the repository root:

    python 05-knowledge/results/rank_geometry_metagraph_20261026/unit_distance/check.py
    python -O 05-knowledge/results/rank_geometry_metagraph_20261026/unit_distance/check.py

The standard-library checker imports no repository code. Explicit raising
conditions survive optimization. It checks 576 full four-state table
compositions; both local hexagon Grams and their joint linear constraints;
the sharp three-dimensional model; all 48 palette maps and all fifteen
pairings of the fixed hexagon; every six-color boundary assignment for
the wheel and rational star; 46,182 complete proper colorings for strips
of lengths one and two at q=3,...,6; and exact strip contacts and transfer
powers through length twenty. The arbitrary-length conclusions rest on
the proofs, not extrapolation. Normal and optimized receipts are saved
beside this note. No large SAT census or external formal build was used.
