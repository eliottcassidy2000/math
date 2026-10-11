# Context can forget a leaf, but counting must keep its weight

2026-10-11. **PROVED:** the elementary abstract interface, weighted extension,
and unit-separation equivalence below. **FINITE-EXACT:** the declared finite
controls. **OPEN:** a six-colour unit-separation equality tester, equivalently
a finite unit-distance graph that is not six-colourable. No plane decision
or universal Collatz consequence is claimed.

## Inheritance and the new target

The [complete boundary relation](boundary_palette_compression_20261011.md)
already remembers all extendible equality partitions of named ports, with
positive lifts. Its signed specialization and
[adaptive integration](plane_cycle_boundary_integration_20261011b.md) retain
actual odd-cycle witnesses. Those are inherited, not re-proved as new work.
The useful extension here is to distinguish three tasks: recover the graph,
decide extendibility under permitted contexts, or count all extensions.

Closest proved mechanism: complete boundary relations. Canonical hostile:
an inaccessible leaf is later exposed as a new port. Corrected near miss:
forgetting colour names without keeping their multiplicities when counting.
Least-used sidecar: a conditional extension weight. The board is **marked
leaf pair / future port access / equality partition / conditional weight /
geometric distance / equality tester**. The general continuation principle
is inherited from [THM-2961, ordered continuation and context tax](../../01-canon/theorems/THM-2961-ordered-continuation-poset-and-finite-group-context-tax.md),
whose linked THM-2176 gives the equality version. We prove the graph instances
directly and claim no priority for standard colouring-interface algebra.

## 1. Two sides of the star/path coin under restricted observation

Let B be the three named vertices of a path P3, with centre c and two leaves.
For r>=0 form S_r by adding r more leaves at c, and I_r by adding r isolated
vertices. Only B may be touched by a future context H; H's interior is disjoint
from the added vertices. All graphs and palettes here are finite, and k>=2.

For each proper k-colouring of H union P3, the extra vertices have exactly
(k-1)^r extensions on S_r and k^r extensions on I_r. Therefore

\[
 \Sigma_k(S_r,B)=\Sigma_k(I_r,B),\qquad
 k^r\,P_{S_r\cup_B H}(k)=(k-1)^r\,P_{I_r\cup_B H}(k).       \tag{1}
\]

The second equality remains meaningful when both chromatic-polynomial values
are zero. Hence the two fragments are interchangeable for Boolean colouring
extension in every permitted context, but not for counting. With six colours
their relative multiplicity is (5/6)^r whenever a ratio is defined.

The [converse codec packet](tournament_converse_context_20261011d.md) proves
that these exact types are paired by transported tournament converse, with
the retained leaf pair marked for a reversible labelled graph operation.
That marker solves a reconstruction problem. If the other leaves are forever
inaccessible, the smaller Boolean interface does not need to reconstruct them.
If their ports may later reappear, the quotient is no longer valid. Already
at r=1, a boundary assignment making c and the extra vertex equal extends
to I_1 and fails on S_1. At k=1 even the unexposed core has no proper colouring;
the stated k>=2 avoids describing that vacuous case as a positive extension.

## 2. Multiplication conditioned on context; addition over alternatives

For a fragment F with ordered boundary B and palette k, let pi be a canonical
partition of B into j<=k blocks. Define w_F(pi) as the number of proper
colourings of F extending the particular boundary assignment pi, using the
fixed names 0,...,j-1. Colour permutations show that the number depends only
on the partition. It can be zero; its positive support is exactly Sigma_k.

There are (k)_j=k!/(k-j)! assignments of distinct actual colour names to these
j blocks. Consequently

\[
 P_F(k)=\sum_\pi (k)_{|\pi|}w_F(\pi).                       \tag{2}
\]

If F and H intersect exactly in B and their union has precisely the declared
edges, fixing the boundary makes their interior choices independent:

\[
 w_{F\cup_B H}(\pi)=w_F(\pi)w_H(\pi).                     \tag{3}
\]

This is a bijection between extension sets: restrict a colouring to the two
interiors, or glue two assignments agreeing on B. It is not a probabilistic
independence assumption. Geometry and any omitted cross-edges remain separate
obligations. Boolean support turns (3) into intersection; keeping integer
weights retains multiplicities as well.

Forgetting one boundary vertex v has a different operation. Fix a partition
pi of the remaining ports with j blocks. Let pi_i be its extension assigning
v to existing block i, and pi_new its extension assigning v a new colour.
Recanonicalize after insertion in the original port order. Then

\[
 w_{F,B\setminus\{v\}}(\pi)
 =\sum_{i=0}^{j-1}w_{F,B}(\pi_i)+(k-j)w_{F,B}(\pi_{new}).    \tag{4}
\]

There are j already named choices and k-j unused colour names. For j=k the
new-block term is zero. For empty remaining boundary, (4) is k times the
one-port weight. A single isolated vertex thus has one orbit but k colourings.
The coefficient k-j is information lost by a bare set projection.

Equations (2)-(4), boundary-edge filtering, and positive witness lifts are a
compositional way to combine additive alternatives with multiplicative
extensions. They do not say that all integer or prime events are independent;
the multiplication proof explicitly uses disjoint interiors after conditioning.

## 3. Abstract minimality has a constructive separating context

For fixed B and k, two fragments have the same Boolean outcome after gluing
every abstract graph context if and only if their complete Sigma_k relations
are equal. Sufficiency is inherited gluing. Here is the converse, filling the
abstract minimality gap without asserting it for geometric contexts.

For any partition pi of B with j<=k blocks, create k fresh anchor vertices
A_0,...,A_(k-1) forming K_k. Connect each port b to every anchor except
A_(pi(b)). The anchors exhaust the palette, so each port is forced to take
the colour of its one non-neighbour. Conversely those assignments are proper.
Thus this context has exactly the singleton relation {pi}. Its conditional
weight is (k-j)!, the permutations of the unused anchors. If two relations
differ at pi, this context makes one union colourable and the other not.

The quantifier is all abstract finite contexts with those named ports and
disjoint interiors. A smaller context bank may identify more fragments;
success on that bank gives no guarantee for a later context outside it.
This construction contains K_k. In the unit-distance plane a clique has at
most three vertices: two unit circles about distinct unit-separated points
have two intersections, whose mutual distance is sqrt(3). Thus the k=6
selector is not an admissible geometric replacement.

## 4. The missing six-colour context is exactly prize-bearing

Fix k>=2. A **unit-separation equality tester** consists of a finite plane
unit-distance subgraph F with distinct nonadjacent named ports u,v at Euclidean
distance one, such that F has a k-colouring and every k-colouring gives u,v
the same colour. These are subgraphs: the geometrically available edge uv is
deliberately omitted. Induced unit-distance graphs cannot have such ports.

There is such a tester if and only if there is a finite plane unit-distance
graph that is not k-colourable.

* From a tester, add its available unit edge uv. Its equality relation and
  the edge's inequality relation have empty intersection, so the union is
  not k-colourable. It is (k+1)-colourable: give u a fresh colour in any
  k-colouring of F.
* Conversely take an edge-minimal non-k-colourable unit-distance subgraph G
  on a finite point set. For any edge uv, F=G-uv is k-colourable by minimality.
  Every such colouring must give u,v equal colours or it would colour G.
  Coordinates and unit separation are inherited unchanged.

This finite equivalence is elementary. Graph compactness additionally
identifies the finite obstruction with failure of a global k-colouring; it
is not used to assert that such an obstruction exists. At k=2 a two-edge
equilateral-triangle path is a positive control. At k=3 every edge deletion
of the exact Moser spindle is a positive control, checked below.

At k=6, importing a geometric equality test at unit separation would already
settle the desired non-six-colourability direction. The abstract anchor
construction cannot be silently substituted. A productive search can now
target a complete two-port relation {equal}, exact distance one, and a
six-colouring lift; it needs no guessed universal boundary-degree bound.
This is equivalent to the finite obstruction problem, not a claimed reduction
in complexity or an easier existence theorem.

The opposite direction also has an exact target: every k-colourable finite
unit-distance subgraph must admit a k-colouring separating each proposed
unit-separated nonedge. If this holds universally, begin with an empty graph
on any finite point set and insert its unit edges one at a time, recolouring
at every step. All finite unit graphs are then k-colourable. This is an edge
release formulation of the universal case, not a proof of its premise.

## 5. A geometric Kempe guard for the two-port search

Let F be an equality tester and fix any proper k-colouring. Its ports u,v
have the same colour a. For every other colour b, they must lie in the same
connected component of the subgraph using only a and b: otherwise swapping
a and b in the component of u gives a proper colouring separating the ports.
In particular all k colours must occur. This supplies k-1 simple alternating
u-to-v paths, each of even length.

At most two channels can have a path of length two. Such a path has a common
unit neighbour of u,v as its middle vertex; two unit circles whose centres
are distance one apart meet in exactly two points. Different channels require
different middle colours. Consequently, for k>=4, at least k-3 channels have
shortest length at least four. The short channels contribute one non-a vertex
each, and each long channel contributes at least two. These contributions
are disjoint across channels because their non-a colours differ. Thus there
are at least 2k-4 non-a vertices. A long simple path has an internal a vertex
in addition to u,v, proving

\[
                    |V(F)|\ge 2k-1 \qquad(k\ge4).           \tag{5}
\]

For a six-colour tester, every colouring therefore has five Kempe channels,
at least three of them long, and the fragment has at least eleven vertices.
This is an elementary necessary filter, not a claimed competitive record
for unit-distance chromatic lower bounds or a sufficient forcing criterion.

The accompanying `kempe_separate` routine searches each two-colour component.
It either returns an actual separated colouring, or the actual paths blocking
all single swaps from the supplied colouring. A path result does not certify
forced equality over other colourings. The distinction has a geometric
four-vertex hostile: put u=(0,0), v=(1,0), and w_+=(1/2,sqrt(3)/2),
w_-=(1/2,-sqrt(3)/2), retaining just the four edges from u,v to w_+,w_-.
The colouring (0,0,1,2) blocks both three-colour channels. Nevertheless
(0,1,2,2) is a proper colouring separating u,v. Recolouring contexts matters
even when every channel of one colouring is present.

The abstract graph K_(k+1) with uv deleted realizes all k-1 channels with
length two. For k>=4 it violates the geometric common-neighbour cap, which
is another explicit reason an abstract equality gadget cannot be imported
as a unit-separation gadget.

## Reproduction and exact controls

[Script](../../04-computation/experiments/boundary_context_weights_20261011d.py),
[saved output](boundary_context_weights_20261011d.out).

```text
python -B -X utf8 04-computation/experiments/boundary_context_weights_20261011d.py
python -B -O -X utf8 04-computation/experiments/boundary_context_weights_20261011d.py
```

The checker covers all 76 labelled graphs of orders0..4, every boundary
subset and k=1..4: 4,396 complete weighted tables, compared with independent
literal colour enumeration and the inherited Boolean compiler. Every
single-port projection in this universe is checked. It also checks116
abstract singleton selectors (boundary size0..4, k=1..6),256 two-port glued
pairs, and320 hidden-leaf contexts (16 context graphs, k=2..6, r=0..3).
Geometric controls retain the inherited exact algebraic spindle incidence;
all11 edge deletions have precisely the equal-port three-colour relation.
Exposed ports, a six-anchor abstract gadget, empty supports, the single-vertex
orbit/count discrepancy, malformed partitions and Boolean integer aliases
are explicit hostiles. An additional32,438 single-swap inputs cover every
proper colouring with k=2..4 of every labelled graph of orders2..4 and every
nonedge. Returned colourings and all returned path witnesses are checked
against the source. The rhombus blocks the false sufficiency implication;
abstract complete-minus-edge graphs at k=4..7 block geometric overreach.
All **51,658** gates remain active under optimization. Independent peers
audited all five proofs, the weighted API, and the later Kempe API. The
optimized replay matches the saved output; eight additional peer controls
cover noncontiguous labels, reversed ports and malformed types. These extra
controls are not included in the production count above. LF SHA256 of the
script is `be7bb71e9a3e8dcacd6ac8b7a612a9c822d8254e9023c32c02c1055c8e98954b`;
the output hash is `a6b1ccc36490a0c68876328ac1abf061cf7cd89f6c7b5e3d3eef576ec049e39e`.
