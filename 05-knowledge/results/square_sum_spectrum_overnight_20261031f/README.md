# Square-sum carriers closed under edge queries; exact tournament renewal prices

**Status: PROVED scoped all-size theorems with exact finite certificates and same-agent proof review. Independent computation paths are supplied; no independent subagent audit or priority claim. Universal tournament raw-H completeness remains OPEN.**

This packet inherits [packet e](../square_sum_spectrum_overnight_20261031e/README.md) at `9a518109cd`: a uniform two-choice weak-pair supply in one phase fibre, exponential square-sum families, native ports, and an exact midpoint-fan deficit grammar. Its finite proof base is112..2811; outside that interval the two choices repeat the same A refinement.

## Inheritance and board

The closest proved mechanism is independent signed-offset composition with fixed label phase. The canonical hostile is still duplicate295/missing293 for different phases. The corrected near miss here is subtler: a synchronized code can be smaller than the Hamiltonian family of its own union graph. The least-used sidecar is the distinction between a primitive local switch and several local switches tied to one bit.

| Lane | Object / operation | Retained information | New result or boundary |
|---|---|---|---|
| Anchor | Square-sum cycles and local edge requirements | Phase, switch owner, atom cuts | Exact required/forbidden-edge compiler |
| Anchor | Union graph of a recursive cycle family | Primitive switches rather than synchronized bits | Fully refined subcubic graph has exactly2^b Hamiltonian cycles |
| Niche | Fixed left/center/right apex decorations | Dyadic leading coefficient and actual positions | Fixed boundary decorations cannot enlarge the pure-power deficit grammar |
| Wildcard | Repeated mixed singleton word | Cyclic correlations and integer period | Strict positive drift, with a sharp all-length minimum |

No fresh incoming origin/main mathematics appeared at the mid-session fetch. The current routes were checked against the existing bounded-run obstruction, THM-4748's fan, and prior necklace/periodic-edit work; their different predicates are kept separate below.

## 1. An exact edge atlas

Fix N and a nonempty A/C schedule w from packet e. Put b=R(w), the product of its25/49 radices. Each b-bit string constructs one cycle. All choices have the same label phase at every interface. Let C_0 be the all-zero reference and E_0 its undirected edge set.

**Theorem 1 (disjoint edge switches).** There are nonempty edge sets R_i,A_i, i=1,...,b, such that

    E(C_z) = (E_0 minus union_(z_i=1) R_i)
             union union_(z_i=1) A_i.                 (1)

Here R_i is contained in E_0, A_i is disjoint from E_0, |R_i|=|A_i|>=2, and the sets of vertices incident to R_i union A_i are disjoint for different i. Thus the edge symmetric-difference distance is exactly

    |E(C_z) triangle E(C_y)|
      =sum_(z_i != y_i) (|R_i|+|A_i|).               (2)

**Proof.** At a finite leaf, an odd segment reversal changes only two boundary edges, and an exchange of two even blocks changes only three boundary edges. It preserves first two labels, endpoint and every label's position parity. Above a leaf, lifting a label x at offset c is a fixed vertex map qx+(-1)^j c, because its parity j does not depend on the choice. Each source has a disjoint label band/offset ownership class. A block transports all its edges; a split tail removes only the invariant edge1--3. Therefore the removed/added edge sets transport by fixed injections, and different leaf supports remain disjoint. Connectors are constant. Induction proves(1).

An unrefined leaf above2811 repeats its one bit in several descendant copies. Its removed/added sets are the disjoint unions of those local changes; different top-level bits still have disjoint support. This explains, rather than invalidates, the larger weights in that case. At N70312, schedule A has25 bits of edge weight100; schedule AA has625 bits of weight4. Equation(2) follows from disjoint supports.

The native two-layer carrier obeys the same rule. Its2--7 cut is between fixed atoms, never inside a varying local switch, and it is constant across choices. Reassemble the two fragments or transport their unchanged internal edges. The prescribed native path, seed and rays remain constant.

`edge_atlas.py` computes the atlas by transporting local differences, without constructing every one-bit flip at the full size. An independent consumer applies the edge replacements, checks square sums and degree two, and walks the resulting Hamiltonian cycle. It agrees with literal path construction on complete64-choice controls, mixed schedules, reversal/exchange leaves and tied larger leaves. `native_edge_atlas.py` does the same through the protected cuts at270112 and270113.

## 2. Exact edge constraints and the cost of retaining ports

A query specifies whether some unordered edge must be present or absent. By(1), each edge has one of four types: always present, always absent, present exactly at z_i=0, or present exactly at z_i=1. Therefore any finite collection of queries reduces to constants and single-bit assignments.

**Theorem 2.** Inside this authenticated family, a query is impossible exactly when a constant is contradicted or one bit is forced to both values. Otherwise, if r distinct bits are fixed, it has exactly2^(b-r) cycle extensions. A witness and every free bit are returned by the compiler.

This is an iff inside the stated family. It is stronger than checking separate local requirements and assuming they glue: contradictions on shared bits are retained. The compiler rejects the requirement that an edge removed by one switch and an incompatible edge added by that switch both be present.

**Corollary (reference-port budget).** Choose any reference member of the family and any k of its edges. At least2^max(b-k,0) members retain them. The same lower bound holds if k arbitrary vertices have their complete reference neighbourhoods prescribed: a changed vertex belongs to at most one bit support, so its two incident-edge requirements cost at most one bit. The selected reference supplies consistency.

This is a uniform statement over the selected ports, not just a few specially chosen vertices. It does not assert that arbitrary edges from the entire square-sum graph are jointly realizable. For the protected family at every N>=270112, b>=2401. The explicit four changed-vertex controls leave2397 free bits while retaining the native path and seed. For growing N, packet e's choice of schedules gives b>=(2/225-o(1))N; retaining k=o(N) arbitrary reference neighbourhoods preserves that exponential rate.

## 3. Refinement closes the actual union graph

Let G be the union of all edges occurring in the family. In general the coarse code need not describe all Hamiltonian cycles of G. At N70312 the A code has25 synchronized bits, while AA unties them into625 bits. Their reference cycles and union graphs are identical. Requiring an old edge in one of these descendant switches and a new edge in another contradicts the coarse bit but has an actual Hamiltonian realization in the same graph. `closed_edge_carrier.json` records the explicit two-edge witness.

This is an obstruction to erasing internal ties, not to the existence of the cycle. The repair is finite and native: append A to the schedule until its leaf index lies in112..2811. The leaf index strictly decreases; its lower bound112 is retained because an A step is appended only when the current index is at least2812. The union graph is unchanged, since each formerly synchronized alternative merely toggled all of these disjoint descendant switches together.

**Theorem 3 (closed subcubic carrier).** For a fully refined schedule with b primitive leaf switches, its union graph G has maximum degree3 and has **exactly2^b undirected Hamiltonian cycles**. These are precisely the decoded cycles. The edge compiler of Theorem2 is therefore an exact Hamiltonian-extension solver for G, not just for a selected subfamily. All fixed native ports remain fixed in the protected version.

The proof uses the local gadgets, not an extrapolation from finite cycle counts:

* **Odd reversal.** Its four changed vertices have variable edges forming a4-cycle. Two opposite ports are joined by the unchanged interior of the reversed segment. That segment has at least three vertices, so its internal degree-two vertices force a common edge at both ports. Each changed vertex has just one common edge available and must use at least one variable edge. The two forced ports must use exactly one. The only possible variable choices are the two perfect matchings of the4-cycle. Thus the original ordering and its reversal are the only Hamiltonian choices through this gadget.
* **Even-block exchange.** The six changed vertices have variable edges forming a6-cycle. The two unchanged block interiors force common edges at four ports. **Every one of the33 exchange certificates has both block lengths at least4**; this stronger condition is checked in the complete bank, not assumed from evenness. There are exactly three local degree patterns. Two are the alternating perfect matchings, giving the two block orders. The third closes a cycle using the two blocks and their neighbouring ports, disconnecting it from the remainder of G. It cannot occur in a Hamiltonian cycle: the global anchored vertex1 is outside that local cycle. Hence only the two intended choices remain.

The changed supports are disjoint. Each unchanged block interior lies within its own primitive leaf copy, so it contains no changed vertex from another gadget; its internal vertices really have degree two in the full union graph. Once the gadget choices are fixed, every common edge is forced. Every independent choice is already known to yield a Hamiltonian cycle by the recursive constructor. This proves both directions and the exact count. A changed vertex has two variable neighbours and one common neighbour; every other vertex has degree two, proving the degree bound.

The local variable-edge enumerations have2 and3 degree-feasible cases respectively. An independent undirected subset-DP count on the contracted gadgets with their forced internal paths gives exactly two Hamiltonian cycles in each. All5400 finite leaf moves were checked for the required interior lengths. These checks certify the finite gadget types; the disjoint-support argument proves the all-size statement.

For every N>=2812, choose the largest admissible A-only depth. It is fully refined, so this produces a concrete subcubic spanning subgraph of the square-sum graph with exactly2^b cycles and b>(2N+1)/5625. The protected version exists at every N>=270112. Mixed radix schedules retain the asymptotic rate2/225. If r of the b primitive moves are reversals and s are exchanges, then

    |E(G)|=N+2r+3s,     r+s=b.

Deleting the invariant anchor edge gives the corresponding exact count2^b of Hamiltonian paths between its two named endpoints. No assertion is made about all Hamiltonian paths with arbitrary endpoints, or all Hamiltonian cycles of the full square-sum graph.

## 4. Tournament progress

[TOURNAMENT-PRICES.md](TOURNAMENT-PRICES.md) contains two new scoped mechanisms:

1. Fixed left/center/right one-apex decorations have an exact quadratic in2^a. If their count is a pure power of two minus a fixed deficit for infinitely many a, dyadic-unit rigidity eliminates the boundary decorations and forces the previous symmetric midpoint cluster. In particular that attempted repair cannot supply deficit9.
2. Every repeated mixed singleton word has strictly positive linear renewal drift. At period length L, its sharp minimum is

       (2^floor(L/2)+2^ceil(L/2)-2)/(2^L-1),

   attained exactly by a single zero-run and a single one-run whose lengths differ by at most one. The proof uses cyclic mismatch correlations. Fixed native prefix/suffix contexts cannot remove the linear term. A changing period may make the price arbitrarily small, so this does not imply a global forbidden value.

The affine classification extends to any fixed finite number of affine apex clusters: positivity forces difference-closed slope support, and reflection then rules out all but the midpoint cluster. [BOUNDED-FAN-COST.md](BOUNDED-FAN-COST.md) strengthens this further. For every fixed deficit outside that midpoint grammar, including9, boundedly many apex losses permit only bounded exponents k in H=2^k-c. An explicit coarse cutoff follows by decomposing the power equation into minimal vanishing sums and retaining their integer polyhedra. Along unbounded one-apex realizations, the number of losses must grow at least on the order of log k/log log k. This does not restrict general two-feedback constructions.

The current exact raw-H coverage through99,999,999 except7,21 is inherited unchanged. The universal spectrum remains OPEN; H/Aut is not substituted for H anywhere in these tests.

## Reproduction and audit boundary

Run from this directory:

```sh
python3 edge_atlas.py
python3 native_edge_atlas.py
python3 closed_edge_carrier.py
python3 periodic_drift.py
python3 three_region_fans.py
python3 affine_cluster_fans.py
python3 vanishing_power_clusters.py
python3 deficit_nine_queries.py
python3 verify_packet.py
```

Normal and optimized Python receipts agree. `VERIFICATION.json` lists universes, positive and hostile controls; `MANIFEST.json` pins inherited proof inputs and packet bytes. Same-agent proof review is supplemented by independent edge-union walking, local undirected subset DP, literal tournament adjacency counts, and exact rational coefficient checks. No second-person or subagent review is claimed.
