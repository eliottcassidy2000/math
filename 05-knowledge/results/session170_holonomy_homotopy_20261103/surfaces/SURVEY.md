# Topological threads: source, connection, next direction, limitation

This is a focused metanalysis of named sources, not a repository-wide census or a literature survey. Navigation entry points verified at the S169 baseline were [S163][S163], [S166][S166], [S167][S167], [S168][S168] and [S169][S169]. Earlier exact files were located with searches restricted to `00-navigation` or `05-knowledge/results`, then read individually. The historical sources are evidence for their stated constructions, not blanket endorsements of every historical research claim.

The new mathematical deliverable is [NOTE.md](NOTE.md): the sharp genus-n realization of full signed-unit transport and the resulting exact capacity of S17's singular finite-field complexes. The computations in `check.py` validate explicit permutations, not just a thematic correspondence. References below to the companion monodromy packet mean [Flatness has a price](../monodromy/NOTE.md), which owns the stronger compatible-quotient classification and physical torus construction.

## 1. Homotopy: singular incidence loops are usable memory

**Source.** [S17 homotopy supplement][homotopy], sections “The incidence graph” through “Integral homology and examples.” It derives K_q ≃ (wedge of c genus-g_p surfaces) ∨ (wedge of b circles), with b=(q/p−1)(q²−1)/2−c+1. At q=9 it obtains thirty spheres and fifty-one circles.

**Connection.** S168 prices full G_n transport in 2n independent generators. The S17 decomposition says exactly which generators survive normalization and which arise from gluing normalization components together. Neither a drawing nor Euler characteristic alone provides this information.

**New useful direction, executed.** Theorem 1 and Corollary 2.1 of NOTE.md pair S168's generators into commuting handle pairs and prove full G_n capacity iff 2c g_p+b≥2n. At q=9, n≤25 comes entirely from incidence loops; at q=11, n≤26 comes entirely from handles. This is a real composition of the two source theorems.

**Limitation.** The homotopy decomposition is not a homeomorphism and does not turn the singular K_q into one surface. The resulting covers are formal transport covers, not an embedding in a tournament-class graph. No computational enumeration of 2^560 sheets is implied by the q=25 capacity bound.

## 2. Holonomy: a scalar section is weaker than a full frame

**Source.** [S62 twisted pair lifts][palette], §§1–3. An unordered local three-pair palette has affine F₂³⋊S₃ transport. Parity observes its two Fano completions. The even subspace E has twisted obstruction dimension 2(β−1)+dim E^H, giving 2β, 2β−1 or 2β−2 according to the permutation holonomy subgroup.

**Connection.** This earlier exact loss computation is the linear-cohomological counterpart of the companion monodromy packet's normal-orbit quotient. It also gives a useful warning: a global section may exist even when the local system is not trivial. In the new setting a quotient cover may extend even though the original sheet cover does not.

**Next direction.** On an actual marked graph, calculate both the normal-orbit loss and a chosen linear observation's H¹ loss. These measure different resources and need not agree. The cheap hostile is two independent cycles with identical three-cycle palette holonomy: the scalar completion obstruction vanishes, but two even-frame obstruction dimensions remain.

**Limitation.** S62 does not realize arbitrary affine palette transports by planar unit-distance exteriors. Its three-cycle sufficiency statement is for one graph cycle, not all graphs. Conflating a section with a trivialization would erase the principal distinction.

## 3. Cover monodromy: topology does not cure one-step erasure

**Source.** [S166 local monodromy][local], §§1–4, and [S168 full transport structure][tower], §§3–5. S166 separates rank loss, hidden sheets and source-language exclusion, and proves rigid base-identity colors can hide arbitrarily many sheets. S168 identifies full monodromy generation, connectedness and deck group as different questions.

**Connection.** A graph has no two-cell relations. Filling cells changes the permitted transport, and the companion monodromy packet computes the finest remaining sheet observation exactly. The action quotient can be substantially smaller than the abstract group quotient because point stabilizers matter.

**Next direction, partly executed.** Use a normal-closure partition algorithm on selected physical loops, rather than inferring global invertibility from rigid colors. For a new physical model, first authenticate local edge-lifting bijections; only then apply this calculation.

**Limitation.** The guard move D is singular and is not a covering transport. A quotient of unmarked tournament observations may also fail local bijectivity. Neither case is repaired by assigning a monodromy group after the fact.

## 4. Braided triple information: actual pair marks exceed support counts

**Source.** [The three-cusp closure theorem][braid], §§2–6. Its specified D₄ Artin marking has three actual cusp pairs and three actual node pairs. The proof transports those pairs by verified simultaneous conjugations/Hurwitz moves, retains the difference between actual retained labels and all inertia-fixed labels, and uses commuting-support incidence to eliminate a formal degree-fifteen survivor.

**Connection.** This is a nonlinear instance of the same “relations cost information” issue. A numerical support configuration can satisfy every Euler and intersection inequality while no permutations realize the required commuting relations. The formal survivor has pairwise singleton intersections; commuting permutations cannot have exactly one jointly moved point because their moved-support intersection is a union of nontrivial cycles.

**Next direction.** When importing a braid-inspired clock, record full local permutation relations and access marks before reducing to orbit or support sizes. The smallest decisive test is to run actual commutators on the candidate support configuration, not to add another scalar moment. Conversely, the surface construction succeeds because its handle commutators are explicitly identity permutations.

**Limitation.** This source is a whole-support exclusion under explicit geometric and marking hypotheses, not a proof of the full planar Jacobian conjecture. It is not evidence that three tournament states inherit Artin braid relations. A braid drawing or a triple count alone supplies no map.

## 5. Associators: an exact algebraic analogue of filling cells

**Source.** [S34 winner algebra and associator bridge][associator], §§2–5. The commutative tournament winner product is associative exactly for transitive tournaments. Its universal associative quotient identifies precisely vertices in the same strongly connected component. Separately, the Z₂ fusion rule admits distinct coherent associators, with ω(a,b,c)=(−1)^(abc) and a compatible braiding r(a,b)=i^(ab).

**Connection.** Both that quotient and the companion monodromy packet's cell quotient are universal repairs: imposing a new law has a computable kernel, which may erase the phenomenon one wanted to study. The pentagon and hexagon require explicit higher data beyond a triple's incidence graph, just as a cover extension requires relator holonomy beyond a graph's edge list.

**Next direction.** Compare the associative SCC quotient and a relator-induced orbit quotient on a common *explicitly constructed* finite carrier. A useful positive target would enlarge the state so a cyclic associator becomes a boundary without identifying all SCC vertices. This requires actual chain data and coherence equations, not a name change.

**Limitation.** The cited source constructs a winner-algebra toy and a pointed fusion example, not a vertex algebra for arbitrary tournaments or a physics equivalence. Its two distinct constructions should not themselves be identified.

## 6. Trivalent quotient: edge multiplicity stores topological capacity

**Source.** [Blue quotient at order nine][blue], §§1 and 4. After black contraction the verified quotient has one hub with 32,760 loops, five leaves, and one trivalent pure-blue vertex attached to the hub by three parallel edges. The same packet's seven-port measurement theorem supplies seven additional physical probes after 105 closures leave a seven-dimensional kernel.

**Connection.** A pair of vertices joined by three individual edges has cycle rank two; its simple-graph shadow has cycle rank zero. Dropping parallel edges therefore deletes the entire local monodromy capacity. This is the same kind of retained-mark obligation as the actual local pair marks in the braid source, but the structures are different.

**New finite deduction.** The hub–trivalent subgraph has β=2. By S166's rigid-ambiguity construction, it supports connected d-sheet formal permutation covers with trivial deck group for every d≥3: assign two independent cycle transports (1 2) and (1 2 … d). Thus that tiny local quotient can carry unbounded formal hidden-sheet ambiguity, even though it cannot carry full G_n transport for n≥3 (which requires β≥2n). The two statements price different requirements. One does not obtain either from vertex degrees alone.

**Limitation.** This deduction concerns formal covers of the retained multigraph, not additional tournament classes. The order-nine census is a finite certificate. It does not prove all-order black-core connectivity. Likewise 112 closure measurements determine the specified 91-coordinate seven-port interface, not every graph invariant.

## 7. Cellular automata: reset maps and invertible holonomy need separate carriers

**Sources.** [Orientation, monodromy and Rule 30][rule30], §§2, 4, 6; [S168 literal cellular covers][cellular], §§1–7.

The first source distinguishes four binary return maps: identity, flip, reset-to-zero and reset-to-one. Only identity and flip belong to invertible covering holonomy. The second realizes an actual C₆→K₃ double cover per tournament cell and supplies a reversible-gate compiler with explicit external cell order and schedule.

**Connection.** A deck bit can record a logical loop but cannot by itself make a clean encoding of an arbitrary irreversible map reversible. A relator-orbit calculation can further determine when those genuine physical deck bits are erased by changing the base topology.

**Next direction.** Before compiling a proposed dynamical rule, classify its return maps as permutations or rank-dropping maps, then measure the history register required by the latter. For reversible rules, test logical return and deck return separately. The inherited off-family widths 2 and 4 double physical periods, whereas the special powers-of-three clock has an anti-period cancellation.

**Limitation.** The Rule 30 random-input trace law does not control the fixed single-seed center. The tournament construction is scheduled, marked dynamics, not autonomous evolution determined by an unmarked tournament or its Hamilton count.

## 8. Jacobian collision tests: differential rank and loop relations answer different questions

**Source.** [S169 collision metanalysis][collision], especially its first three ledger rows; [S166 local monodromy][local], §§1, 5–6. The former supplies precise collision-complete physical test banks and distinguishes ambient branch velocities from legal tournament deformations. The latter gives an explicit real-analytic constant-Jacobian noninjective map, and an invertible three-count guard observer whose allowed source image is still restricted.

**Connection.** A complete first-order test cannot establish global source membership or trivial covering holonomy. Conversely, a monodromy calculation is not a test for rank loss. The new surface construction deliberately addresses only the hidden-sheet case.

**Next direction.** For any candidate new tournament operation, test in order: whether the observation is well-defined; whether its differential/kernel is controlled; whether local transitions are bijective; whether relator holonomies are paid; and whether the desired output belongs to the actual construction language. A failure at one stage names the missing information rather than demanding a larger undifferentiated state vector.

**Limitation.** These examples are not polynomial Keller-map counterexamples or a reduction of the planar Jacobian problem. “Locally invertible” and “constructible by this positive operation language” remain independent predicates.

## 9. Unit distance: a marked common graph is stronger than an analogy

**Source.** [S166 shared-root diamonds][diamonds], §§1–7. Its D(H) has injective planar unit realizations exactly for linear forests H, while allowing intermodule collisions enlarges the class to all bipartite H. The Moser spindle occurs as an explicit induced seven-class subgraph in the strong order-six tournament single-arc-flip graph. Its longer diamond chain is rigid but eventually nontraceable.

**Connection.** This source gives an unusually strong cross-domain transfer: the same fully certified finite graph appears in two settings. It also isolates a global admissibility coordinate, injectivity, that neither local edge constraints nor infinitesimal rigidity captures. That is analogous in role, not identical in mathematics, to source-language membership and covering-sheet data.

**Next direction.** Take one actual flip-graph embedding, retain its representative/edge witnesses, and compute the relator holonomy induced by literal relabelings around its cycles. This could connect the current formal cell-descent theorem to a second physical example beyond the hexagon. No such lift is claimed in this packet.

**Limitation.** The common graph is the single-arc-flip class graph, not the blue/black complement-line metagraph. Its embedding proves no plane-color bound and no unit-spine theorem for extremal configurations.

## Priority conclusion

The most productive new bridge is homotopy decomposition → commuting generator pairs → exact cover capacity; it is completed in NOTE.md. The next useful bridge is an authenticated physical cycle → explicit relator permutation → exact observation quotient. Broad transfers from homology rank, symmetry, braid vocabulary, or scalar response counts should wait until this map and its loss are written down.

[S163]: https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/00-navigation/S163.md
[S166]: https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/00-navigation/S166.md
[S167]: https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/00-navigation/S167.md
[S168]: https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/00-navigation/S168.md
[S169]: https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/00-navigation/S169.md
[homotopy]: https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/05-knowledge/results/session17_prime_surfaces_20261020/independent_audit/HOMOTOPY_SUPPLEMENT.md
[palette]: https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/05-knowledge/results/session62_agent_theme_bridges_20261025/holonomy/TWISTED_PAIR_LIFTS.md
[local]: https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/05-knowledge/results/session166_planar_metagraph_20261103/local_monodromy/NOTE.md
[tower]: https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/05-knowledge/results/session168_computational_structures_20261103/tower/STRUCTURE.md
[braid]: https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/05-knowledge/results/planar_jc48_sep08_three_cusp_closure.md
[associator]: https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/05-knowledge/results/session34_gauge_representations_20261023/vertex/vertex_simplices_associator_bridge.md
[blue]: https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/05-knowledge/results/blue_quotient_20261020b/README.md
[rule30]: https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/05-knowledge/results/orientation_monodromy_rule30_20261026/README.md
[cellular]: https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/05-knowledge/results/session168_computational_structures_20261103/cellular/NOTE.md
[collision]: https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/05-knowledge/results/session169_runner_jacobian_20261103/collision/META_LEDGER.md
[diamonds]: https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/05-knowledge/results/session166_planar_metagraph_20261103/unit_diamonds/NOTE.md
