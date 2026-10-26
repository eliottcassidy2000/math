# Four endpoint states require a retained palette

2026-10-26. **PROVED in the stated scopes; independent
[audit PASS](../information/SIX_COLOR_AUDIT.md).** The
statements below have elementary proofs and exact replay controls. They
provide scoped obstructions to a proposed gadget compiler, not a new
Hadwiger–Nelson bound or a six-colour equality gadget. The finite results
are **FINITE-EXACT** in their explicitly enumerated universes.

## Inheritance and connection contract

The starting objects are the full boundary relations in
[Session62 palettes, sections 1–3](../../session62_agent_theme_bridges_20261025/palettes/RESULT.md),
the four-state group `E ⋊ S3 = S4` in
[the paired-frame packet, sections 1–3](../../hadwiger_six_seven_eight_20261025/fano_eight/RESULT.md),
and the exact prescribed six-palette gadgets in
[the cube-boundary packet, section 5](../../hadwiger_six_seven_eight_20261025/cube_boundary/RESULT.md).
Their full relation and local-coefficient results are inherited. The new
work asks which of their states can actually be distinguished by an
**unprescribed** graph fragment.

The canonical hostile is the induced triangle with three pendant ports:
it enforces NAE at palette three, but it does not enforce it at palette six.
The corrected near miss is retaining lists after the vertices prescribing
them have become ordinary uncoloured vertices. The least-used sidecar is
the exact action of global `S6` on the proposed state representation.
The existing equality-test gate is retained: a nonvacuous six-colour
equality gadget whose omitted edge has unit-separated endpoints is
equivalent to a finite seven-chromatic unit-distance graph. See
[Session61 palette-deck, section 5](../../session61_hadwiger_six_seven_20261025/palette-deck/RESULT.md)
and [THM-4661, odd-degree carriers and equality testers](../../../../01-canon/theorems/THM-4661-odd-degree-carriers-a-failed-six-colour-exterior-has-at-least-twenty-vertices-and-equality-testers-are-seven-colour-certificates.md).

Anchor: unprescribed six-colour realizability. Niche: chromatic thresholds
for three-port relations. Wildcard: a small exact geometric realization of
the conditional holonomy hostile. The live board is **four endpoint states /
global colour symmetry / full boundary partition / fresh colours / exact
unit edges / existential pair assignment**. The symmetry probe makes a
single-frame table constant; the recolouring probe locates the first
possible nonconstant unprescribed relations; the metric probe then shows
why a genuine drawing still does not remove the existential pair assignment.

| Source → target | Map and preserved predicate | Lost data and required sidecar | Decisive test |
|---|---|---|---|
| Proper six-colourings → three-port boundary table | Restrict a colouring; support means extension | Interior forgotten; retain every joint equality partition and its count if counting | Compare all five three-port partitions |
| Endpoint frame `x in F2^3` → three colour values | `x_i` chooses one endpoint of the fixed disjoint pair `P_i` | All eight tuples become one global `S6` orbit; a palette reference must remain accessible to distinguish them | Enumerate all 720 colour permutations |
| Two paired frames → difference kernel | Retain the cross-frame equality partition, equivalently `x+y` | Individual endpoint origins disappear; keep actual pair assignments and shared-port identifications | Test all 24 four-state permutations against diagonal translations |
| A fixed pair assignment → six-colouring | `c(v)=2q(v)+b(v)` | One `q` is only one sector; keep or quantify over every `q` | Exact seven-point example has an inconsistent sector but 1,728 successful sectors |

No literature priority is claimed. The arguments use ordinary colouring
symmetry, recolouring, and the elementary centralizer of a regular abelian
group. The contribution is their explicit application and exact metric
controls at this particular six-colour interface.

## 1. A single unprescribed frame has no four-state constraint

Let `G` be any finite graph with three distinct ordered ports `u1,u2,u3`.
There are no prescribed colours or named lists anywhere in `G`. Write
`w_G(a1,a2,a3)` for the number of proper six-colourings extending that
specified boundary assignment. Every global permutation of the six colours
acts bijectively on extensions, so `w_G` depends only on the equality
partition of the three values.

Fix disjoint colour pairs `P_i={a_i^0,a_i^1}`, whose union is the palette,
and put `a(x)=(a_1^{x1},a_2^{x2},a_3^{x3})` for `x in A=F2^3`.
Every `a(x)` is an injective triple. All injective triples lie in the same
`S6` orbit, of size `6*5*4=120`. Therefore

    w_G(a(x)) is independent of x in A.                         (1)

In particular its support on the eight frames is either empty or all
eight. On either four-point completion `p(x)=x1+x2+x3=eta`, it is empty or
all four, with constant weight. This holds regardless of `chi(G)` or any
Euclidean realization. It is not just a failure of parity compression.

**Scope.** The theorem concerns a fragment accessed through these three
vertices only. A graph with additional retained palette ports has a larger
joint relation and is not covered by (1). A variable describing a pairing,
a Fano completion, or a palette reference is additional state; it does not
become an enforceable condition because a drawing or an abstract CSP names
it. A graph may also have four equality-pattern states for unrelated
reasons. Four NAE equality patterns are not the four endpoint states here:
the endpoint tuples are all in the single all-distinct pattern.

## 2. Two unprescribed paired frames permit only translation kernels

Now retain two ordered triples of ports with the same assigned pairs
`P_i` in the corresponding positions. Restrict the six-port boundary table
to `a(x),a(y)`. Swapping the two colours in pair `P_i` is a global palette
permutation, so it flips coordinate `i` in both frames. Consequently

    w(x+b,y+b)=w(x,y) for every b in A,
    w(x,y)=f(x+y), where f(d)=w(0,d).                           (2)

The eight values of `d` are exactly the eight equality partitions appearing
in this restricted slice: in coordinate `i` the two port colours are equal
iff `d_i=0`, and cross-coordinate colours never agree. Thus (2) follows
both from the group action and from the complete partition table.

Restrict both frames to the even completion
`E={000,110,101,011}`. The same proof for `b in E` gives a difference
kernel on the four-point additive group `E`. There are exactly `2^4=16`
possible translation-invariant Boolean relations on `E x E`, all of form

    R_D = {(x,y): x+y in D}, D subset E.                        (3)

This enumerates invariant tables, **not** realizable graph gadgets. If the
table is the graph of a total deterministic map `T:E→E`, (2) forces
`T(x+b)=T(x)+b`, hence `T(x)=x+T(0)`. Exactly four of the 24 permutations
of the four states pass this necessary test: the translations. The other
20 elements of the abstract `S4` cannot be enforced by such an
unprescribed fragment with these fixed corresponding pair assignments.
For example a nonidentity order-three map fixing `0` fails: any equivariant
map fixing `0` must fix every `x`.

This does not forbid physical reordering of ports or different local pair
assignments. Those change the chart and must be carried explicitly. With
arbitrary fixed pairs and three ports on each side, relabeling coordinates
does not by itself prove that a nontranslation table is the boundary
relation of a graph. Nor does a table passing (3) automatically supply the
pair-domain constraints. In particular, this is a stricter typing test for
the earlier abstract `S4` compatibility system, not a global impossibility
theorem for larger unit-distance gadgets.

## 3. A sharp chromatic staircase for three independent ports

Here no pairs or lists are prescribed. Let `G` be `k`-colourable and let a
desired proper boundary equality partition have `t` nonsingleton blocks
and `s` singleton blocks. Starting from a proper `k`-colouring:

1. Recolour each nonsingleton block with its own fresh colour. Its vertices
   are independent because the desired boundary assignment is proper.
2. Keep the singleton vertices in their old colours except that repeated
   singleton colours are separated with distinct fresh colours.

At most `t+max(0,s-1)` fresh colours are used. Nonsingleton blocks have
different new colours; singleton blocks are all different; old interior
colours are unchanged. Every edge stays proper. A final global permutation
matches any particular named boundary assignment with that partition.
Thus every such assignment extends using at most

    k+t+max(0,s-1) colours.                                   (4)

For three independent ports there are five partitions, written
`000,001,010,011,012`. Formula (4) has the following exact consequences:

* If `chi(G)<=4`, every six-colour boundary assignment extends.
* If `chi(G)<=5`, the four nondiscrete partitions always extend. The only
  possible omitted orbit is `012` (all three distinct).
* If a six-colourable graph has exactly the NAE relation (it omits `000`
  and admits the other four partitions), then `chi(G)=6`. More strongly,
  `chi(G-B)>=6`: if the exterior were five-colourable, colouring all three
  independent ports with a fresh sixth colour would extend `000`.

The first conclusion is about support; extension counts can still vary
between the five equality partitions. The hypotheses refer to `G` itself,
except in the explicitly stronger last necessary condition. They are not
the fixed-eight-port field theorem
[THM-4614, three-colourable contexts](../../../../01-canon/theorems/THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md),
whose different boundary geometry remains essential.

These bounds are sharp for abstract graphs. Joining `K4` to three
independent ports gives a five-chromatic graph with six-colour weights

    (120,24,24,24,0)                                         (5)

on the five partitions. Joining a new `K3` to all six vertices of the
triangle-with-three-pendants gadget gives a six-chromatic graph with weights

    (0,48,48,48,12).                                         (6)

For (5), the clique needs four unused colours. For (6), the three added
clique colours leave exactly the original three-colour NAE problem, with
`(6-r)_3` choices for the added clique when the ports use `r` colours.
Both examples contain `K4`, so neither is a planar unit-distance graph.
The Gram matrix based at one of four mutually unit-separated points would
have diagonal `1`, off-diagonal `1/2`, and determinant `1/2`; it has rank
three and cannot arise from three vectors in the plane. This supplies a
metric hostile, not a missing geometric realization.

**Gluing corollary.** Suppose a graph is a literal union of fragments with
private interiors, each four-colourable, and each meeting the rest in at
most three ports. Include every seam edge in a separately retained port
graph `H`; there are no undecomposed edges between private interiors.
If `H` is six-colourable, the union is six-colourable: choose its colouring
and independently extend across each fragment using (4). Boundary edges
inside a fragment are included in `H`, so their restrictions are proper.
Repeated use, cyclic incidence, or higher-order summaries of these
fragments cannot introduce an extra six-colour obstruction. A first
nontrivial three-port fragment must leave this four-colourable class, or
have more retained ports. Arbitrary cross-interior unit edges invalidate
this decomposition and must be included before applying the statement.

## 4. An exact unit-distance graph realizes the conditional hostile

Let `alpha=(1+i sqrt(3))/2`. Put

    (A0,A1,A2)=(0,1,alpha),
    (B0,B1,B2)=(10,11,10+alpha), Z=100.                       (7)

All seven points are distinct. Their complete induced unit graph is two
disjoint equilateral triangles and one isolate. This follows exactly from
`N(a+b alpha)=a^2+ab+b^2`; the checker tests all 21 pairs.

Conditionally prescribe pair `P1` at the A vertices, pair `P2` at the B
vertices, and pair `P3` at `Z`. At cell `j` write its three endpoint bits as
`x_j=(a_j,b_j,z)`. Each seam from `j` to `j+1 mod3` consists of the unit
edges `A_j A_(j+1)` and `B_j B_(j+1)` and the literal common vertex `Z`.
Its exact conditional relation is

    x_(j+1)=x_j+(1,1,0).                                    (8)

The gain has even parity, so each seam preserves either four-state
completion. Yet the three seam gains sum to `(1,1,0)`, and the resulting
frame system has no section. In contrast, the actual unprescribed graph
has the explicit three-colouring

    (c(A0),c(A1),c(A2),c(B0),c(B1),c(B2),c(Z))
       =(0,1,2,0,1,2,0).                                   (9)

There are exactly `86400=(6*5*4)^2*6` proper six-colourings.
The graph therefore realizes the abstract conditional obstruction with
actual complete unit-distance geometry, while disproving the implication
from that obstruction to non-six-colourability.

The missing quantifier can be counted, not merely named. Fix the palette
pairs `{0,1},{2,3},{4,5}` and let `q:V→{0,1,2}` assign a pair to each
vertex. The inherited decomposition `c=2q+b` says the compatible bit
assignments are proper two-colourings of the subgraph `H_q` retaining
same-`q` edges. For each triangle exactly `3^3-3=24` assignments of `q`
avoid a monochromatic odd cycle. There are therefore

    24^2*3=1728 successful q among 3^7=2187.                 (10)

The complete census has 972, 648, and 108 successful assignments with
respectively 5, 6, and 7 components of `H_q`. Their bit-extension counts
sum to `972*2^5+648*2^6+108*2^7=86400`, agreeing with an independent direct
six-colour backtracking count. The fixed `q` used in (8) is among the 459
unsuccessful choices. Keeping `q` existential is indispensable.

This is the ordered gain hostile already present abstractly in Session61,
now realized by an exact seven-point graph and audited across *all* pair
assignments. It is not a new gain-cohomology theorem. No equality between
distinct vertices is forced: the third coordinate is the same vertex `Z`.

## 5. Removing the prescriptions from the previous exact gadgets

The [previous coordinate receipt](../../hadwiger_six_seven_eight_20261025/cube_boundary/receipt.json)
contains two exact unit graphs `F` and `G`. The current checker reads those
coordinates as input, independently recomputes every pair distance, verifies
every claimed edge and absence of an edge, and checks that every formerly
fixed-colour vertex is a leaf. No pin is kept in this calculation.

For the six-vertex triangle-with-three-pendants core, the ordinary
six-colour weights on `000,001,010,011,012` are

    (60,68,68,68,71).                                        (11)

This follows by inclusion–exclusion on the three spoke inequalities, or
by direct exact backtracking. Every freely coloured leaf contributes a
factor five. The additional separate centre contributes six. Consequently

    w_F = 6*5^25 * (60,68,68,68,71),
    w_G = 6*5^15 * (1,1,1,1,1).                             (12)

All five boundary orbits survive. On the disjoint-pair endpoint slice,
`F` has the constant weight `6*5^25*71`, and `G` is also constant, exactly
as (1) requires. The previous prescribed binary NAE table remains correct
for its original prescribed problem; it becomes a different table when
the prescriptions are removed. This is an explicit conversion receipt,
not a retraction of the earlier scoped result. The graphs retain their
exact 32/19 vertices and 31/15 edges.

## 6. What remains necessary for a genuine six-colour obstruction

The actual advance is a sharper interface target. A single endpoint frame
with only three ports is invisible to unprescribed colour symmetry. Two
fixed paired frames permit only the difference kernels (2), with just four
deterministic maps on a completion; the full abstract `S4` is not freely
available. A palette reference or a larger jointly retained boundary must
carry whatever distinctions a construction uses. An ordinary four-colourable
three-port fragment is transparent at palette six, and a six-colour NAE
fragment requires a six-chromatic exterior. These are necessary gates,
not sufficient existence claims.

The exact geometric example shows that even a legitimate unit-distance
implementation of a conditional frame loop does not remove the need to
quantify over all pair assignments. A genuine obstruction still requires
the complete unprescribed boundary relations, consistent joint palette
alignment, all actual cross-unit edges, and an empty final composition.
Neither this packet nor the inherited prescribed gadgets supply it.

## Exact replay

Run in this directory:

```text
python check.py > receipt.json
python -O check.py > optimized.json
```

The two JSON outputs agree byte-for-byte. The receipt includes the checker
hash, the prior coordinate-receipt hash, all seven new coordinates and all
their unit edges, an actual three-colouring, all exact extension counts,
and the consequence objects above. There are 6,728 always-active checks.

The declared finite universes are: all 720 colour permutations; all 64
pairs of eight-point frames; all 65,536 Boolean relations on two four-state
frames; all 24 permutations of four states; all 2,187 ternary assignments
on the seven-point graph; all five three-port partitions on both abstract
sharp controls and on the unit examples; and 2,205 constructive recolouring
controls (every `k=1,...,6` initial three-port colouring and every target
partition, in the densest compatible graph with one interior representative
of each old colour). The latter is a test of the construction, not an
exhaustive census of finite graphs. The previous unit gadgets are tested
on every coordinate pair, with all leaves unpinned. The proofs establish
the unrestricted graph statements; the finite universes serve as positive
and hostile controls. No search over all unit-distance graphs is claimed.
