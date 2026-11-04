# Constructive synthesis of three uploaded graph-theory papers

All source claims below are taken from the supplied PDFs, not independently reproduced proofs. PDF page numbers are one-indexed. The 4CT paper's printed page is PDF page minus one; the hypohamiltonian paper uses matching printed and PDF pages; Thomassen's journal pages begin at 90. Separate source-issue notes contain narrowly scoped corrections and a proposed implementation repair.

## 1. A boundary-state view of reductions

The most useful common idea is to separate a large graph into a small patch, its boundary, and an untouched exterior. One need not remember the patch's internal construction: remember the finite information the exterior can observe.

- In the 4CT paper, the interface is a ring coloring, together with its Kempe reducibility level. The local objective is eventually to make that coloring extendible.
- In the hypohamiltonian paper, the interface is a spanning path-forest's endpoint decomposition. The objective is to preserve non-Hamiltonicity while supplying a compatible cycle after each vertex deletion.
- In Thomassen's paper, equal-length cycle changes may fail to improve length; the extra information is the cycle's bridge structure and an auxiliary weight. This suggests that a search should preserve relevant interface state while using a secondary potential to orient otherwise neutral moves.

This is a methodological bridge, not a theorem identifying these different state spaces.

## 2. Precise source results: linear reducibility in triangulations

**Inoue, Kawarabayashi, Miyashita, Mohar, Thomassen, Thorup.** *The Four Color Theorem with Linearly Many Reducible Configurations and Near-Linear Time Coloring*. arXiv:2603.24880v3. Uploaded title page dated October 5, 2026. [Source](https://arxiv.org/abs/2603.24880).

- Theorem 1.1, PDF p. 2: an O(n log n) planar 4-coloring algorithm.
- Section 2, PDF p. 4: an obstructing cycle has length 3 or 4 with at least one vertex on each side, or length 5 with at least two on each side. “Internally 6-connected” means there are no such cycles. A separating 5-cycle around a single degree-5 vertex is therefore not obstructing.
- Lemma 3.2, PDF p. 9: their set contains exactly 8202 D-reducible configurations. Each has at most 19 vertices, ring length at most 18, diameter at most 4, and all ring colorings are 25-extendible. Except for one radius-3 configuration, each has radius at most 2. Only its designated center can have degree greater than 8, and then its degree is at most 12.
- Theorem 3.3, PDF p. 9: in linear time one finds linearly many pairwise non-touching objects, each a separating 3- or 4-cycle or a configuration from the set. Every involved vertex has degree at most 12.
- Corollary 3.4, same page: one obtains either linearly many non-touching **induced** D-reducible configurations or linearly many non-touching obstructing cycles. Here length-5 obstructing cycles enter when non-induced configurations are converted to short cycle obstructions. Do not conflate this corollary with the 3/4-cycle statement of Theorem 3.3.
- Lemma 3.1, PDF p. 8: the diameter bound alone is not sufficient to turn a non-induced occurrence into an obstructing cycle. The lemma also needs a specified condition on length-4 boundary paths and external degree.
- Theorem 5.1, PDF p. 14: in an internally 6-connected triangulation, if every vertex in B12(v) has zero final charge, B14(v) contains an induced reducible configuration. Thus zero curvature is not an obstruction to finding reductions. “Flat” refers to the paper's **final discharged charge**, not simply degree 6.
- Theorem 6.5, PDF p. 16: positive-charge vertices and zero-charge vertices of degree at least 9 yield reductions or short separating cycles in degree-bounded radius-2 neighborhoods; an all-zero-charge, degree-at-most-8 radius-12 neighborhood also does so. Degree-bounded balls avoid high-degree hubs destroying packing estimates.
- Theorem 14.2, PDF p. 55: improving the color-extension subroutine improves overall coloring through a constant-factor recursion. Section 15 announces linear-time coloring and a bounded deletion set on fixed surfaces as forthcoming work; they are **not established results of this uploaded paper**.

**Constructive direction:** For the user's triangulation work, distinguish three layers in a certificate: (1) a local occurrence with degree/embedding conditions, (2) inducedness or a concrete short obstruction, and (3) boundary-state extendibility. These separate checks make local replacements composable. The paper's flat-neighborhood idea argues for sampling and classifying large zero-final-charge patches instead of only chasing positive-charge defects.

## 3. Precise source results: equal-length cycle moves with a potential

**Carsten Thomassen.** *Chords in longest cycles in 3-connected graphs*. J. Combinatorial Theory B 179 (2026), pp. 90–117. [DOI](https://doi.org/10.1016/j.jctb.2026.03.001).

- Main result, journal p. 90 / PDF p. 1: **some** longest cycle in every 3-connected graph has a chord. This is weaker than the older conjecture that **every** longest cycle has a chord.
- Theorem 2, journal p. 92 / PDF p. 3: in a connected multigraph with no adjacent even-degree vertices, if all vertices outside a cycle C have even degree, there is a distinct cycle containing every odd-degree vertex.
- Theorem 3 and Corollary 1, same page: the Fleischner–Stiebitz 3-coloring theorem for a Hamiltonian cycle plus disjoint triangles supplies an independent transversal containing one vertex from each triangle.
- Corollary 3, journal pp. 93–94 / PDF pp. 4–5: a chordless cycle whose vertices have degree at least 3, with a bridge having three private attachments in every bridge-component, can be replaced by a cycle at least as long. Only selected private attachment vertices are lost. In the equal-length case, the exchanged portions have very restricted lengths.
- Theorem 4, journal p. 95 / PDF p. 6: the stronger technical statement assumes 2-connectivity, cycle-vertex degrees at least 3, a 3-fan to the cycle from every outside vertex of degree at least 3, and no degree-2 vertex adjacent to two vertices of degree at least 3. Then either a longer cycle exists or an equally long chorded one exists. It also preserves specified attachments of even-star bridges.
- Journal p. 91 / PDF p. 2 explicitly says polynomial running time would depend on polynomial algorithms for the lollipop step and the relevant cycle-plus-triangles coloring problem. The paper does not supply an unconditional polynomial algorithm for finding a chorded cycle at least as long as an arbitrary input chordless cycle.

**Constructive direction:** A Hamiltonian or near-Hamiltonian search can treat equal-length reroutings as useful rather than failures, while recording bridge signatures and looking for a strictly monotone secondary objective. This is a proposal for search design, not a guarantee that an arbitrary simple bridge-size score works; the paper explicitly explains why greedily maximizing the largest bridge can fail (journal p. 95).

**Direct consequence with limited novelty:** A hypohamiltonian graph has longest-cycle length n−1 and is 3-connected, so this theorem ensures a chorded (n−1)-cycle. This fact is a cheap certificate target, not a route to constructing a hypohamiltonian graph. In fact the stronger chord observation for hypohamiltonian graphs is elementary from minimum degree ≥3: each vertex on an (n−1)-cycle has at most one neighbor off the cycle, and a chordless such cycle would force every cycle vertex to be adjacent to the omitted vertex, creating a Hamiltonian cycle by an insertion. Therefore this particular application does not need the deep theorem.

## 4. Precise source results: hypohamiltonian patch replacement

**Durias, Guerrero, Jaiyen, Kaplan, Lapine, Peralta, Shrivastav.** *Local modifications in planar hypohamiltonian graphs*. arXiv:2610.21210v1, October 2026. [Source](https://arxiv.org/abs/2610.21210).

- Theorem 1.2, p. 2: planar hypohamiltonian graphs exist at every order n≥36.
- Proposition 2.2, pp. 5–6: each of the three displayed patch replacements H→Hi preserves Hamiltonicity in both directions. “Contains the patch” includes the prescribed external edges and excludes additional edges to patch vertices (Remark 2.1); ordinary subgraph containment is insufficient.
- Theorem 2.4, pp. 6–7: the replacements preserve hypohamiltonicity. H→H1 is the +2 replacement; H→H3 is an alternative iterable +4 replacement. H2 demonstrates further local flexibility.
- The interface consists of exactly 14 endpoint decompositions. All four boundary terminals lie in one bipartition class; that class has two more vertices than the other. Every spanning cycle fragment therefore has exactly two path components, including possible isolated terminal components. There are 12 single-path-plus-singleton patterns and 2 noncrossing two-path patterns. A crossing pairing is excluded by planarity.
- Corollary 2.10, p. 8: the imperfect +2 replacement with one edge removed is valid if two specified exterior Hamiltonian paths exist, both from terminal 1 to 3, one retaining all four terminals and one retaining only terminals 1 and 3. This is an important distinction between universal local replacement and replacement requiring an exterior witness.
- Theorem 2.13, p. 10: the +3 modification transfers Hamiltonian cycles/paths **forward**. It is not a general hypohamiltonicity-preserving operation. Its internal deletion cases and non-Hamiltonicity need extra information in the instances used.
- Proposition 2.17, pp. 13–14: for the generalized Thomassen +4 operation, two adjacent cubic vertices on the quadrilateral suffice for forward preservation of hypohamiltonicity. One cubic vertex suffices for the stated converse. Proposition 2.19 adds two exterior Hamiltonian-path conditions to obtain a forward result with only one cubic vertex. Do not freely apply +4 to every facial quadrilateral.
- Counts established by the constructions: at least 12 planar hypohamiltonian graphs of order 36, 8 of order 37, 66 of order 38, and 48 of order 39 (Corollaries 2.12, 2.14, 2.5, 2.15).
- Theorem 2.22 and Corollary 2.23, pp. 15–16: four planar hypohamiltonian graphs of orders n1,...,n4 give a planar hypotraceable graph of order Σni−6. Consequently planar hypotraceable graphs exist at order 130 and every order ≥132.
- Corollary 2.24, p. 16: almost hypohamiltonian planar graphs exist at order 31 and every order ≥33.
- Proposition 4.1, p. 21, is restricted to pentagulations of the particular five-quadrilateral patch H with its prescribed attachment behavior. It does not prove the general conjecture that order 34 is minimal. That conjecture and its three supporting heuristics remain conjectural (pp. 22–23).

## 5. A derived, reusable patch-transfer lemma

This is an elementary reformulation developed here from the papers' method, not a new theorem claimed by the sources.

Let P and Q be finite graph patches with the same labeled terminal set T. Glue either to any exterior X by identifying T, with disjoint edge sets and no other shared vertices. Assume X has at least two vertices outside T, so that at least one remains after an exterior vertex deletion. A feasible fragment is a spanning linear forest of the patch in which every nonterminal has degree 2 and terminals have degree 0,1,2. Exclude closed cycle components. Its signature consists of the terminal degree vector and the matching pairing the degree-1 terminals along nontrivial path components. Degree-0 terminals record isolated components. Write S(P) for the signature set, and use the same definition after deleting a vertex.

**Lemma.** Suppose:

1. S(P)=S(Q).
2. For every terminal t, S(P−t)⊆S(Q−t).
3. For every nonterminal v of Q, some nonterminal w of P satisfies S(P−w)⊆S(Q−v).

Then if X∪P is hypohamiltonian, X∪Q is hypohamiltonian.

**Proof.** A hypothetical Hamiltonian cycle in X∪Q restricts to a fragment in Q. Replace it with an equal-signature fragment in P. Terminal degrees remain 2 after gluing, and the induced pairing of exterior path components is unchanged, so the result is one Hamiltonian cycle, contradiction. For deletion of an exterior vertex use S(P)⊆S(Q) and the old deletion cycle. For deletion of terminal t use condition 2. For deletion of an internal vertex v of Q choose w from condition 3 and transfer the cycle of (X∪P)−w. Each transfer preserves the exterior and yields one cycle for the same degree/pairing reason. These cover all vertex deletions.

No planarity is needed for the lemma. Planarity reduces the state space to noncrossing pairings and permits embedding-preserving gluing. The two-exterior-vertex assumption excludes a Hamiltonian cycle living entirely inside the patch both before and after any one exterior vertex deletion.

### Concrete independently generated example

`boundary_signature_demo.py` transcribes H and H1 from Figure 2 using an outer 8-cycle with terminals 0,2,4,6. It exhaustively enumerates only the small patches' spanning linear forests, not global graphs. It uses Python's standard library, requires no installation, and outputs `boundary_signature_demo.json` with actual fragment witnesses.

The run found exactly 14 full signatures on each patch, equal full signature sets, all four terminal-deletion inclusions, and at least one old internal deletion witness covering every new internal deletion signature requirement. The JSON contains edge lists, coverage choices, and representatives. This is a reusable prototype for automatically certifying future replacement patches; it is not an independent verification of the paper's full construction or minimality conjecture.

## 6. Connections to the other research themes

**Square-sum:** The path-signature language applies to Hamiltonian-path searches in any graph, including a square-sum graph. A direct gadget transplant is unjustified: new arithmetic labels create edges to many existing labels and generally do not preserve the small interface or planarity. The sound next question is whether naturally occurring small separators in a concrete square-sum instance admit reusable endpoint-state certificates. This is a research direction, not an asserted structural property.

**Tournament versus even-automorphism graph equinumerosity:** None of these papers proves or directly treats that identity. Interface equivalence preserves Hamiltonian behavior, not automorphism groups, orbit sizes, or enumeration. A potentially useful but unproved direction is an equivariant signature calculus: retain the boundary group's action on states and require compatible patch maps before attempting unlabelled transfer. This must be matched to the exact definition of “even automorphism” in the user's project. No bijection follows from the 14-state calculation alone.

**Best next constructive step:** Search bounded planar bipartite four-terminal patches with imbalance 2, classify them by full and vertex-deletion signature sets, and retain inclusion-minimal certificates. The target is new universally safe +k replacements or a clear list of minimal exterior witnesses for imperfect replacements. This converts the paper's hand-drawn modifications into a finite, auditable discovery problem and supplies a concrete connection to the triangulation project's reducible-configuration libraries.
