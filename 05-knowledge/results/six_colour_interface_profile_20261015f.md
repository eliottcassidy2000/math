# Six-colour extension depends on how many classes touch the ports

2026-10-15. **PROVED:** the interface-colour extension lemma and its
conditional application to the existing eight ports in `Q(zeta55)`.
**FINITE-EXACT:** palette options and the cross-intersecting pair census.
**OPEN:** an authenticated six-colour forcing context and Hadwiger–Nelson.
Self-audited; no novelty claim for the product-of-palettes method.

**Continuation:** [safe singleton palettes](compression_contexts_20261015g.md#3-a-safe-singleton-extends-457-of-the-marked-boundary-patterns)
show that a four-colouring with at most three interface classes preserves
all 457 boundary patterns using at most five colours. This is another
necessary forcing screen, not a constructed four-chromatic field graph.

[Component palette matching](six_colour_component_pairs_20261016a.md)
further retains independent two-class flips and the exact nine-point
incidences. Two admissible pairs suffice, even when they overlap; all
patterns using at most four colours extend across any four-colourable
exterior, without the port-free-class assumption.

## 1. Inheritance and the new coordinate

[THM-4614, three-colourable contexts](../../01-canon/theorems/THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md)
already proves that every six-colouring of the marked pentagon/rhombus
ports extends across any three-colourable exterior in `Q(zeta55)`.
Its exact, context-independent census has nine possible double-port
vertices and no exterior vertex with three port neighbours.

The closest mechanism is that theorem's disjoint-pair extension. The
canonical hostile is a large exterior that leaves every boundary pattern
unchanged. The corrected near miss is to confuse four-chromaticity of the
exterior with sufficient boundary forcing. The least-used coordinate is
**which exterior colour classes actually meet a port neighbour**.

We also recovered the stronger full-relation lesson from
[the triangle-product repair](unit_triangle_product_boundary_20261013c.md):
one retained chord preserves the complete three-colouring relation, but
the assertion fails with a fourth colour. The
[hinge-rotation study](hn_type_m_hinge_rotations_and_forced_radii_20261014.md)
likewise finds that centre-monochromatic supports alone miss necessary
cross-colour constraints. The correct consumer remains the full labelled
boundary relation at the specified palette.

## 2. A general extension lemma

Let `G` be a finite graph with marked ports `B`, and let the exterior
`H=G-B` have a supplied proper `k`-colouring. Suppose at most `ell` of
its colour classes meet `N_G(B)`. Every exterior vertex has at most two
port neighbours. Let `d` be the number of exterior vertices with two port
neighbours of **distinct prescribed colours** in a fixed proper
`m`-colouring of `G[B]`.

**Proposition H1 (PROVED).** If

\[
                 k+\ell\le m,\qquad d<\binom m2,         \tag{1}
\]

that prescribed port colouring extends to `G`.

**Proof.** Give each of the `ell` interface classes a disjoint pair of
colours, and every other exterior class one disjoint colour. A uniformly
random injection of these labelled slots into the `m` colours makes each
particular interface pair uniform among all `binom(m,2)` pairs. A vertex
fails precisely when its two differently coloured port neighbours exhaust
its class's pair. Its failure probability is `1/binom(m,2)`. The union
bound in (1) leaves an unblocked assignment. Choose an available colour
at every vertex. Exterior edges join different classes and hence disjoint
palettes; port edges are proper by construction. QED.

Classes with a singleton palette have **no** port neighbours. This is the
extra hypothesis that prevents the incorrect inference from ordinary
four-colourability alone. Repeated class/pair blockers need only be counted
once, so counting vertices is a convenient sufficient bound.

## 3. Apply it without changing the eight geometric ports

Let `B` be exactly the existing pentagon plus rhombus, sharing 0, from
[the audited census](six_colour_core_doubling_20261012.md). Let `H` be any
finite unit-distance exterior in `Q(zeta55)`. Its port-neighbour bound is
two and its entire double-port bank has size nine, independent of how the
exterior was generated.

For **every** proper six-colouring of these same ports, extension follows
from any of the following supplied exterior colourings:

| total exterior classes | classes allowed to touch the ports | palette assignments | assignments certainly surviving nine blockers |
|---:|---:|---:|---:|
| 3 | 3 | 90 | 36 |
| 4 | 2 | 180 | 72 |
| 5 | 1 | 360 | 144 |

The first row is inherited THM-4614. The last two are the new screens.
These counts retain labelled class slots and arbitrary prescribed port
colours; they are not counts of actual exterior colourings.

Define the anchored chromatic profile

\[
 \kappa_\ell(H;B)=\min\{k:\ H\text{ has a proper }k\text{-colouring
 with }N_G(B)\text{ using at most }\ell\text{ classes}\},
\]

with value infinity if the set is empty. A context that eliminates even
one of the bare six-colour boundary patterns must therefore satisfy

\[
       \kappa_3\ge4,\qquad\kappa_2\ge5,\qquad\kappa_1\ge6. \tag{2}
\]

The first condition is ordinary non-three-colourability. The second says
more: if the exterior has a four-colouring, **every** such colouring must
use at least three classes on its port-neighbour set. The third excludes
a five-colouring with that set monochromatic. Neither added condition is
asserted to be sufficient for forcing.

This gives a cheap next search predicate before computing all 570 boundary
patterns. Test the actual exterior with the extra restriction on its
port-neighbour classes. A witness satisfying one of the table rows gives
an extension compiler for *every* boundary pattern, so more expensive
forcing tests on that context are unnecessary. A failed restricted
colouring test merely passes this necessary screen; it is not a unit
distance or six-colour forcing certificate.

## 4. Retain blocker support, not just its total

In the `(k,ell,m)=(4,2,6)` case, let `F_1,F_2` be the forbidden pairs
for the two interface classes. A palette assignment exists exactly when
some pair outside `F_1` is disjoint from some pair outside `F_2`; the two
remaining colours fill the singleton classes. Thus failure of this palette
method means the two families of allowed pairs are cross-intersecting.

The complete finite universe of nonempty subfamilies of the fifteen pairs
of six colours has 32,767 members. For each first family, intersecting its
pair-neighbourhoods gives the maximal possible second family. Exact census:
if both families are nonempty and cross-intersecting, their combined size
is at most 10 (36 first families attain that maximum). Consequently a
failure with both classes retaining at least one allowed pair needs at
least 20 distinct class/pair blockers. If a class has no allowed pair,
15 blockers in that class already obstruct this palette method.

These are **FINITE-EXACT** support statements, stronger than a scalar
union bound in this small universe. They do not say that a graph with
15 or 20 such blockers has no six-colouring: a colouring need not respect
the supplied class palettes. The one-class bank of all fifteen pairs is
an explicit hostile to confusing failure of a sufficient method with
failure of the underlying colouring problem.

## 5. Connection and reproduction

Source: a literal unit graph, its marked ports and a proper exterior
colouring. Map: the colouring to labelled class palettes and forbidden
pairs. Preserved predicate: existence of a lifted colouring for each
specified boundary assignment. Lost information: exterior colourings
outside those palettes and their multiplicities. Needed sidecar: the
actual port incidences and the witness lift. No geometric vertices or
edges are inferred from an abstract palette obstruction.

The common move with [native Collatz packets](collatz_reset_debt_20261015f.md)
is to retain the fibres where a repair is actually allowed. Native signs
may cancel only at a common authenticated terminal; colour blockers matter
in their assigned interface classes. Equal scalar totals do not determine
either outcome. This is a precise analogy of certificate design, not a
map from colourings to Collatz or a transfer of either open conjecture.

[The control script](../../04-computation/experiments/marked_debt_controls_20261015f.py)
enumerates all 90/180/360 labelled assignments and all 32,767 first pair
families. It checks exact blocker counts, nine-blocker witness guarantees,
and the palette-method hostile. Normal and `-O` agree on
[the output](marked_debt_controls_20261015f.json).
