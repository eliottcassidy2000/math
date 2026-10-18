# Exceptional phases and the thirteen pair boundary threshold

2026-10-18, `clock_duality_oct18`. **PROVED:** the eight-anchor projection obstruction and the algebraic event sieve below. **FINITE-EXACT:** the exhaustive supergraph census. The elementary graph facts used in the proof are not claimed as new general graph theory. The new application is a quantitative filter for this repository's marked cubic construction. **OPEN:** a unit-distance configuration retaining an equal-endpoint six-colouring and eliminating all separating six-colourings.

## Inheritance and the research target

The anchor is the Hadwiger–Nelson boundary-support construction; the niche is marked snark transport; the wildcard is the primitive cyclic support recovered from PR93. The concept board is: full boundary support, fibre projection, exceptional phase, shade transport, primitive cyclic address.

The closest proved mechanism is [THM-4116, boundary-state gluing](../../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md): compatibility requires the ordered boundary table, not merely a nonzero count or a minor. The corrected near miss is the family of reflection-linked copies which still retracts onto its source. The hostile is the three-edge-colourable Desargues graph containing a Petersen subdivision, authenticated in [the marked snark packet](snark_boundary_transfer_20261018a.md). The least-used relevant sidecar is the projective chart for a unit direction: its missing point must survive any Pythagorean or modular parametrization.

[THM-4661, odd-degree carriers](../../01-canon/theorems/THM-4661-odd-degree-carriers-a-failed-six-colour-exterior-has-at-least-twenty-vertices-and-equality-testers-are-seven-colour-certificates.md) proves an exterior floor of20 and gives an exhaustive, not proof-checked SAT floor of22 in its stated triangle-free carrier. These do **not** automatically apply to the different field
`F=Q(zeta15,sqrt(-11),w)` of [THM-4700, unit-path ring](../../01-canon/theorems/THM-4700-the-unit-path-ring-of-the-six-colour-field-short-reciprocal-walks-and-the-heron-census-of-coefficient-triangles.md), which contains `sqrt(-3)`. The tests below retain that distinction.

The imported PR's four commits and nine final files are preserved on `main`; [the provenance record](pr93_preservation_20261018a.json) retains the original commit hashes and file hashes. The verified [Git bundle](pr93_original_history_20261018.bundle) also preserves the four original commit objects and their exact intermediate versions after the PR branch is deleted. Its base prerequisite is already in main history. Captured local cherry-pick IDs refer to the pre-rebase preservation step, not permanent mainline commit IDs.

## A projection with twelve new pairs cannot force equality

Let `B` be the eight-vertex graph with edges

```
01,12,23,34,40,05,56,67,70.
```

It is `C5` and `C4` sharing vertex0. It has9 edges, degree4 at0 and degree2 at each other vertex. Its570 proper six-colour partitions are the current complete marked port support.

Let `G` be any finite graph containing these eight marked vertices. Suppose a map

```
pi : V(G) -> V(B)
```

fixes all eight marked vertices and never identifies the endpoints of an edge. Define `Q` on the eight anchors by retaining every edge of `B` and every distinct projected edge `{pi(u),pi(v)}` of `G`.

**Theorem.** If `|E(Q) \ E(B)| <= 12`, then for every original nonedge `xy` of `B`, `G` has a proper six-colouring in which `x` and `y` have different colours. The number of vertices of `G` is unrestricted.

**Proof.** A graph `H` on eight vertices with chromatic number at least7 contains a `K7`. Indeed, two disjoint missing edges permit two paired colour classes and four singleton classes. Otherwise the edges of the complement are pairwise intersecting, hence form a star or a triangle (or a subset of one). A missing triangle permits a three-vertex colour class and five singleton classes. In the remaining star case, the other seven vertices form a clique.

If `Q` had no separating six-colouring at `xy`, then `Q+xy` would require at least7 colours and contain a `K7`, say on `W`. The old graph has at most7 edges in `W`, because deleting any vertex removes at least2 of its9 edges. At most one edge of this clique is the newly inserted `xy`. Thus `Q[W]` has at least20 edges and contributes at least `20-7=13` pairs outside `B`, a contradiction. Pull back a separating colouring of `Q` by `pi`. Every edge of `G` projects to an edge of `Q`, and the marked endpoints are fixed. This is the desired colouring. QED.

This proof also handles a noncolourable quotient: such a quotient cannot occur under the twelve-pair bound. If the `K7` omits a terminal or `xy` was already a projected edge, the lower bound in the argument is at least14 instead of13.

**Sharp abstract boundary.** Omit any degree-two vertex `v`. On the other seven vertices add every pair except a chosen old nonedge `xy`. Keep only the two old edges at `v`. This adds exactly13 pairs and gives `K7-xy` plus the eighth vertex. The terminals are forced equal, a proper six-colouring exists, and the eighth vertex can always be coloured. There are `7*14=98` such minimum augmented graphs. They are abstract projected graphs; their dense cliques are not planar unit-distance realizations.

The bound counts **distinct projected pairs**, not physical cross-edges. It is conditional on the eight-anchor projection and the absence of collapsed edges. Extra anchors, fibre edges producing loops, or coincident geometric points with incompatible anchor labels require separate analysis. Passing this necessary filter does not prove that the geometric lift forces equality: the lift can colour different members of a fibre differently even when the quotient cannot.

## From arbitrary phase search to finitely many algebraic events

Let `F` be a fixed conjugation-stable field containing all anchors, coefficients and map entries. Let `g_i` be specified Möbius maps preserving the unit circle. In a chosen embedding consider

```
z_(c,i)(q) = c + k (g_i(q)-g_0(q)),    q*conj(q)=1.
```

The zeroth copy is the literal marked source. Clearing denominators in

```
N(c'-c+k(g_j(q)-g_i(q)))-1
```

gives a polynomial of degree at most4. Poles and introduced zero factors must be excluded, and collisions checked independently. [The exact cubic packet](cubic_boundary_search_20261018a.md) implements this operation over the48-dimensional coefficient field and proves the generic retraction for its `b,c` coefficients.

Fix an irreducible event factor `f` and a root satisfying the chosen norm-one embedding and all pole/collision guards. Every nonzero event polynomial vanishes there exactly when `f` divides it. Form `Q_f` from all projected pairs whose events vanish, including identically present edges. Any separating six-colouring of `Q_f` pulls back to the complete geometric graph at this phase. Therefore a candidate with a loopless eight-anchor projection can be rejected before the full colouring search if its event factor supplies at most12 distinct new projected pairs.

For the fixed finite maps and anchors, the nonzero event factors form a finite list. Away from that list only the generic edges remain. This is a complete event reduction for that specified one-parameter family, not for arbitrary unit-distance configurations. New fields, new maps and new anchors enlarge the family and its event list.

For rotations `g_i(q)=alpha_i q` the sharper equation is quadratic. Set

```
d=c'-c, v=k(alpha_j-alpha_i).
P(q)=d*conj(v)+(N(d)+N(v)-1)q+v*conj(d)q^2.
```

At a valid unit phase, `P(q)=0` is exactly the unit-distance equation. If all its coefficients lie in `F`, an element of degree3 over `F` cannot satisfy a nonzero such polynomial. One must not call a phase cubic over a smaller field while silently placing its own cubic coefficients in `F`.

For the designed exceptional edge `d=u-v`, with `N(u)=1`, this factors as

```
P(q)=(q-1)(Cq-A),   A=d*conj(v), C=v*conj(d).
```

When `A,C` are nonzero, the other root is `q*=A/C`, also of norm1. These are the two circle-intersection branches, with possible tangency when they coincide. This gives exact exceptional phases instead of a floating-point phase mesh.

## What the larger cubic tests actually establish

For each `k=b,c` the exact geometric packet tests four phase maps `q,u q,u^2 q,u^3 q` and the eleven-anchor base consisting of the old eight ports and

```
d_j=u^j-k(u^j-1), j=1,2,3.
```

At `q=1` there are44 distinct points,39 unit edges, and three designed cross-edges. This is33 exterior vertices relative to the full eleven-anchor base, or36 beyond the original eight marked ports. All570 original boundary partitions extend. The three new edges genuinely break the obvious layer projection, but a repaired retraction sends each new spoke tip to the old neighbour1 of vertex0 and preserves every source colour.

All six other event roots, three for each coefficient, were also tested. Each gives44 distinct points and37 edges: four copies of the old nine-edge graph and one new spoke. The same boundary-preserving retraction works. This closes these particular designed quadratic events; it does not enumerate every event factor of every larger family.

The precise next target is a **simultaneous event factor** whose projected constraints remove separating support, while the full graph still has an equal-endpoint six-colouring. A new edge or a large vertex count alone does not meet that target. The thirteen-pair theorem applies when there is a loopless projection onto the original eight anchors; the eleven-anchor quotient has to be tested on its own full support.

## Marked reductions linking the three lanes

The connection with snarks is an exact interface prescription. [The snark packet](snark_boundary_transfer_20261018a.md) retains pair charge, shade and wire parity when suppressing paths; a bare Petersen minor does not retain those data. Its two-defect counts can be contracted with a proposed connector's full six-colour boundary table. The missing object is a geometric connector whose table removes all incompatible separating states while keeping an equal state.

The modular connection has the same explicitly limited role. [The PR93 extension](pr93_modular_bridge_20261018a.md) retains the primitive cyclic address in `P^1(Z/p^r)` and the differential frame. In the stated inert CM setting, the primitive trace is divisible by `p^(2r)`; restoring a scalar backtracking branch at depth2 contributes `p^2` and destroys precision4. Neither a modular trace nor a finite norm-one circle is a colouring certificate.

The Pythagorean connection is the projective norm-one parametrization. Over `Q(i)`, `(1+it)/(1-it)` gives rational unit directions, with `t=infinity` supplying the omitted point `-1`. The modular two-chart map `[z] -> z/conj(z)` is its finite-ring analogue. The useful preserved structure is the full projective address and involution, not a numerical coincidence among primes. In all three lanes the cheapest decisive test is to restore the supposedly dispensable branch or label and check whether the target predicate survives.

## Reproduction and exact universe

Run:

```
python -B -X utf8 04-computation/experiments/exceptional_factor_support_20261018a.py --out 05-knowledge/results/exceptional_factor_support_20261018a.json
```

The script exhausts all `2^19=524288` simple supergraphs of the fixed nine-edge source and all19 original nonedge terminal pairs. It computes support from the570 proper six-colour partitions and independently compares it with the eight-vertex `K7` characterization: **9,961,472 exact comparisons**. There are206 empty-support quotients and2808 graph/terminal equality cases with nonempty support. The minimum is13 added pairs at every terminal pair. The98 minimum examples retain four partitions in91 cases and five in7 cases; witnesses are stored in [the JSON receipt](exceptional_factor_support_20261018a.json).

Positive and hostile controls are the unaugmented full-support graph, `K7` minus one edge, and closing that missing edge. An independent proof audit confirmed the projection conditions and the sharp13-pair count. No geometric realizability or universal Hadwiger–Nelson conclusion is inferred from this finite universe.
