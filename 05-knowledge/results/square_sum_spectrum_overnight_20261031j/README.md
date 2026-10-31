# Mixed square-class trees and the exact boundary state beyond trees

**Status: PROVED scoped composition and boundary formulas; FINITE-EXACT independent controls.** Universal raw tournament-H coverage and arbitrary interval-ownership extensions remain OPEN. This packet extends [packet i](../square_sum_spectrum_overnight_20261031i/README.md), without changing its frozen bytes.

## Inheritance: what is reused, and what is added

The general named-boundary transfer principle, including singleton paths, weights and a connectivity test, is already proved in [Session32, PORT_TRANSFER.md, section3](../session32_square_monster_20261023/ports/PORT_TRANSFER.md). Section2 below specializes that stronger transfer to cycles with prescribed removed edges; it is not claimed as a new general dynamic-programming theorem. The new arithmetic content is that packet i realizes every required tree edge with one source copy, and different tree edges can use different square classes. The local cyclic-matching count and the marked atlas restriction below make the specialized state space explicit.

The tournament [Session40 forest rule, sections4--5](../session40_marked_ports_square_sum_20261023/marked-ports/PATH_TRANSPORT_AND_FOREST_RULE.md) already distinguishes removable tree phases from retained joint cycle holonomies. Its states are affine-F3 attachment phases, whereas ours are terminal path matchings. We inherit the tree/circuit distinction without equating those state spaces. The least-used sidecar recovered here is the [105 matching/atlas dictionary](../hamming_matching_atlases_20261018e.md), whose marked eight-slot map gives the exact restriction in section3.

## 1. Pairwise arithmetic types suffice on a tree

Let T be a finite tree. At each vertex v put a finite positive-integer square-sum carrier G_v and designate a different edge for each incident edge of T. Condition the source language jointly on all those edges being present. For each edge vw of T, assume the two designated edges have endpoints in one common square class d_vw. Different edges of T may use different classes. Ports may share vertices within a carrier; they may not reuse the same edge.

Then square injections of the individual sources and two seams per tree edge produce a square-sum carrier with:

- exactly the sum of the source vertex counts;
- Hamiltonian cycles in bijection with the product of the jointly port-conditioned source cycle languages;
- unchanged degree at every cut endpoint, hence preserved subcubic degree;
- preserved unused ports, source labels through explicit injections, individual source weights, and additive branch clocks.

**Proof.** Join distinct tree components one edge at a time. The pending port in either component is the square-scaled image of its original edge and retains its square class. It has not been removed, because each original edge is assigned only once. Apply packet i's one-copy theorem. All other pending ports survive; this also holds when two ports share a vertex, since each removed incidence is replaced by one seam incidence. The new two-edge cut forces the inverse cycle factorization. Induct over the tree edges. At the end, every global cycle decodes uniquely to the original port-conditioned cycles. In particular no consistency condition is silently deferred to a closing circuit: a tree has none. This is an exact bijection, not a lower bound from independently chosen examples.

Joint conditioning matters. If two requested ports require contradictory source bits, the local language is empty and so is the product. Separate nonempty unary port lists do not certify the joint source language. A port's name is not a consumable resource distinct from its edge: registering the same physical edge twice is rejected.

The weighted formula uses each source cycle's internal-edge weight after deleting all its consumed ports. The output sum is the product of these source sums times the product of all seam weights. This is valid over any commutative coefficient ring, including signed weights. Owner, scale and source-bit addresses are retained for inversion; the resulting graph is a selected sparse carrier, not a full induced graph or an initial interval.

### Native consequence and independent controls

`tree_joins.py` joins four actual source carriers through classes 1, 5 and 10, with different primitive shapes at the first join. The result has **11297 vertices and exactly 2^97 Hamiltonian cycles**. Three 2812-vertex sources each pay their bit-22 condition and contribute 24 free bits; the 2861-vertex source contributes 25. Eight complete branch words are walked and decoded back to all four source edge sets. The count for every word follows from the input closure theorem and the tree bijection, not from these eight walks.

The Euler-brick roots 44,117,240 provide a three-vertex square-sum carrier because their pairwise hypotenuses are 125,244,267. Its three ports share vertices. A three-leaf star of four such triangles is joined using different port shapes; the final arithmetic cycle has twelve distinct vertices. This explicitly checks that vertex sharing is allowed while edge reuse is not.

`audit_tree.py` is an independent graph-only consumer. Subdividing one K4 edge three times gives a seven-vertex subcubic source with two cycles and four forced ports. Literal Hamiltonian-cycle enumeration checks both path and star attachment trees under all eight endpoint-pairing choices: sixteen contexts, 256 output cycles, with every source recovered. It does not use the arithmetic producer's cycle-walk routine.

## 2. What a cyclic attachment must retain

Suppose a source cycle has r distinct designated edges. Give the two endpoints of each edge separate **incidence slots**, even if two slots refer to the same physical vertex. Delete the designated edges. The remaining cycle breaks into r paths, some of length zero when designated edges are adjacent. Their endpoints define a perfect matching P of the 2r slots. Let M_0 pair the two slots belonging to each original deleted edge. Then P union M_0 is one alternating cycle.

Conversely a matching P is a possible abstract single-cycle boundary type exactly when P union M_0 is connected. The number of these abstract types is

    2^(r-1) (r-1)!.

To count, order the r labelled deleted edges cyclically and choose their endpoint orientations. Fix the first edge's orientation to quotient reversal. There are (r-1)! cyclic orders and 2^(r-1) remaining orientations. This counts abstract boundary types; a particular carrier need not realize all of them. At r=1,2,3,4 the counts are 1,2,8,48, compared with 1,3,15,105 unrestricted perfect matchings.

For several sources and a fixed seam pairing M, take the disjoint union of their path matchings P_v. Assume the source vertex sets are disjoint and the proposed seam edges are distinct, loop-free and authenticated in the intended graph; otherwise the slot model alone does not certify a simple-graph realization. The sewn graph is one Hamiltonian cycle **if and only if**

    kappa((union_v P_v) union M)=1,

where kappa counts connected components in the matching multigraph. Each alternating component expands into exactly one output cycle. Zero-length paths must be retained as incidences; erasing them loses this equivalence when ports share a vertex.

Aggregate each source's internal path weights by its exact matching, obtaining W_v(P_v). The exact cyclic composition formula is

    (product of seam weights) *
    sum_(P_v) [kappa((union_v P_v) union M)=1] product_v W_v(P_v).

On trees the connectedness indicator is identically one on admissible source types, by the preceding induction. On interaction graphs with circuits it is an additional constraint. This identifies the missing boundary coordinate; it does not assert arithmetic seam feasibility around arbitrary circuits.

### Three-source cyclic incidence hostile

Take three abstract triangles. Each offers two adjacent ports; the remaining edge and the shared vertex give one nontrivial path and one zero-length path. Join them around an interaction triangle. There are eight endpoint-pairing words. Exactly the four even-parity words give a nine-cycle. The other words give components of sizes 4+5 or 3+6. All local cycles exist and every final vertex still has degree two. The first failed implication is therefore "locally valid two-seam joins imply global connectedness" once the interaction graph ceases to be a tree.

The geometric/arithmetic status of those cyclic seam systems is not asserted. They are a graph-level hostile for an unconditional product argument. `boundary_matchings.py` independently reproduces their cycle counts from the twelve incidence slots and exhaustively checks the abstract type formula through r=6. For a three-leaf tree all eight admissible center matchings pass, as the theorem predicts.

## Board

* **Anchor:** interval absorption must preserve the local matching or prove a cheaper complete boundary relation; numeric port values alone do not suffice.
* **Niche:** compile cyclic arithmetic attachments with exact matching tables and retained zero-length paths. First separate the arithmetic compatibility problem from the connectedness selector above.
* **Wildcard:** tournament composition already needs source state and path data. This packet gives another exact instance where local counts multiply on a tree but a circuit forces a retained boundary state. It supplies no universal tournament or Collatz convergence implication.

This is a proof by explicit maps and inverse maps, with a precise tree/circuit boundary. No novelty claim is made for classical matching-transfer methods or Euler bricks.

## 3. An exact 48-of-105 bridge to the earlier Hamming atlases

For r=4, label the eight incidence slots by F2^3 so that the original port pairs are M_0={{X,X+a}:X in F2^3}, with a fixed nonzero a, each pair counted once. Retain this labelling, its origin and basis. The earlier matching/atlas dictionary sends every perfect matching P on these eight **dual slopes** to an optimal four-chart partition of the 32-word code H8+I2. Thus it also sends exactly the 48 cyclic source-boundary types to a distinguished subset of those 105 atlases. The preserved extra predicate is connectedness of P union M_0; it is not a comparison of two unmarked lists of the same size.

The prior decomposition 105=7+42+56 is by one, two or four different pair differences. The new connectedness restriction has the exact split

    48 = 0 + 24 + 24.

No parallel atlas works: its two constant translations with M_0 generate at most a four-element subgroup, so the eight slots do not form one alternating cycle. For a two-direction atlas with directions d,e, let H=span(d,e). Its two H-cosets use d and e respectively. If a lies in H, both matchings preserve the cosets and the union is disconnected. If a is outside H, following the alternating edges visits all eight slots: the four-step displacement is d+e, and the first half visits distinct H coordinates and both cosets. There are four planes not containing a, three unordered direction-bases per plane, and two assignments to the cosets, giving 4*3*2=24. The remaining 48-24=24 connected matchings have four directions.

`boundary_matchings.py` enumerates all105 and checks this restriction and its full component table. This gives a concrete test for a four-port cycle compiler inside the previously authenticated atlas domain. It does not geometrically realize a six-colour selector. Changing or forgetting the marked slot-to-slope identification changes the dictionary, just as the prior atlas paper's physical-translation hostile warns.
